
#include "sensor_filter.h"
#include "sensor_filter_mat.h"
#include "UTM.h"

namespace mobili::ins {

// coordinate values
static int utm_zone_;
static double origin_x_;
static double origin_y_;

// filter values
static int64_t current_timestamp_;
static float state_[9];  // x, y, z, rx, ry, rz, vx, vy, vz
static float state_rotmat_[9];
static float covariance_[81];
static float gravity_[3] = {0, 0, 9.81};

void InitFilter() {
  utm_zone_ = 0;
  origin_x_ = 0;
  origin_y_ = 0;
  current_timestamp_ = 0;

  // initialize the initial state
  std::fill_n(&state_[0], 9, 0.0f);
  std::fill_n(&state_rotmat_[0], 9, 0.0f);
  std::fill_n(&covariance_[0], 9 * 9, 0.0f);
  state_rotmat_[0] = 1.0f;
  state_rotmat_[4] = 1.0f;
  state_rotmat_[8] = 1.0f;

  // initialize position covariance
  covariance_[0] = 1e2;
  covariance_[10] = 1e2;
  covariance_[20] = 1e-5;
  // initialize rotation covariance
  covariance_[30] = 1e-2;
  covariance_[40] = 1e-2;
  covariance_[50] = 1e4;
  // initialize velocity covariance
  covariance_[60] = 1e-2;
  covariance_[70] = 1e-2;
  covariance_[80] = 1e-5;

  InitBuffers();
}

// x, y, z, weight
constexpr float kAccStableThreshold = 0.02;
static float recent_acc_average_[4];
static int64_t last_static_timestamp_ = 0;
void AddStaticMeasurement(int64_t timestamp) {
  if (timestamp - last_static_timestamp_ < 1e8) {
    return;
  }
  last_static_timestamp_ = timestamp;

  InsertVelocityMeasure(0, 0, 0, 1e-2, 1e-2, 1e-2);
}

bool CheckStatic(float* acc_input) {
  float acc_norm = std::sqrt(acc_input[0] * acc_input[0] + acc_input[1] * acc_input[1] +
                             acc_input[2] * acc_input[2]);
  // there might be error with acc
  if (acc_norm > 100.0f) return false;

  // check if the acc is local stable
  bool flag_stable = false;
  if (recent_acc_average_[3] > 0.0f) {
    flag_stable = std::abs(acc_input[0] - recent_acc_average_[0]) < kAccStableThreshold &&
                  std::abs(acc_input[1] - recent_acc_average_[1]) < kAccStableThreshold &&
                  std::abs(acc_input[2] - recent_acc_average_[2]) < kAccStableThreshold;
  }

  // compute the variance of acc to tell if car static
  float weight_inv = 1.0f / (1.0f + recent_acc_average_[3]);
  recent_acc_average_[0] =
      weight_inv * (recent_acc_average_[3] * recent_acc_average_[0] + acc_input[0]);
  recent_acc_average_[1] =
      weight_inv * (recent_acc_average_[3] * recent_acc_average_[1] + acc_input[1]);
  recent_acc_average_[2] =
      weight_inv * (recent_acc_average_[3] * recent_acc_average_[2] + acc_input[2]);
  if (recent_acc_average_[3] < 10) {
    recent_acc_average_[3] += 1;
  }

  if (!flag_stable) return false;
  // ACC might has large error
  // if (std::abs(acc_norm - 9.81) < 0.1) return true;

  return true;
}

// buffer some static values used in the calculation
static float acc_input_[3];
static float acc_world_[3];

static float gyr_input_[3];
static float acc_input_hat_[9];
static float gyr_input_mat_[9];

static float temp_matrix_1_[9];
static float temp_matrix_2_[9];
static float temp_matrix_F_[81];
static float temp_matrix_F_T_[81];
static float temp_matrix_l_1_[81];

void InitBuffers() {
  std::fill_n(recent_acc_average_, 4, 0.0f);
  std::fill_n(temp_matrix_F_, 9 * 9, 0.0f);

  temp_matrix_F_[0] = 1;
  temp_matrix_F_[10] = 1;
  temp_matrix_F_[20] = 1;

  temp_matrix_F_[60] = 1;
  temp_matrix_F_[70] = 1;
  temp_matrix_F_[80] = 1;
}

int InsertImu(int64_t timestamp, float ax, float ay, float az, float gx, float gy, float gz) {
  if (current_timestamp_ == 0) {
    current_timestamp_ = timestamp;
    return 0;
  }
  // update the state
  float dt_sec = 1e-9 * static_cast<float>(timestamp - current_timestamp_);
  if (dt_sec < 0.0f) return -1;
  current_timestamp_ = timestamp;

  /////////////////////// state update ///////////////////////
  // update position
  VectorAdd(&state_[0], &state_[6], dt_sec, &state_[0], 3);

  // compute acc in world coordinate
  VectorAssign(acc_input_, ax, ay, az);
  MatrixVectorMultiple(state_rotmat_, acc_input_, acc_world_, 3, 3);
  VectorSub(acc_world_, gravity_, acc_world_, 3);

  // update velocity
  VectorAdd(&state_[6], acc_world_, dt_sec, &state_[6], 3);

  // update rotation
  VectorAssign(gyr_input_, dt_sec * gx, dt_sec * gy, dt_sec * gz);
  ComputeRotationMatrix(gyr_input_, gyr_input_mat_);
  // HatPlusIdentity(gyr_input_mat_, gyr_input_);
  MatrixMultiple(state_rotmat_, gyr_input_mat_, temp_matrix_1_, 3, 3, 3);

  Hat(acc_input_hat_, acc_input_);
  MatrixMultiple(state_rotmat_, acc_input_hat_, temp_matrix_2_, 3, 3, 3);
  std::memcpy(state_rotmat_, temp_matrix_1_, sizeof(state_rotmat_));

  // get the rotation vector, and normalize
  // AxisAngleFromRotationMatrix(state_rotmat_, &state_[3]);

  /////////////////////// covariance update ///////////////////////

  // assign F matrix
  temp_matrix_F_[6] = dt_sec;
  temp_matrix_F_[16] = dt_sec;
  temp_matrix_F_[26] = dt_sec;

  // TODO(yeliu) : review the jacobians
  // InverseRightJacobian(&state_[3], temp_matrix_1_);
  // float* RwT = acc_input_hat_;
  // Transpose(gyr_input_mat_, RwT, 3, 3);
  // MatrixMultiple(temp_matrix_1_, RwT, gyr_input_mat_, 3, 3, 3);
  //
  // Vector3Assign(&gyr_input_mat_[0], &temp_matrix_F_[30]);
  // Vector3Assign(&gyr_input_mat_[3], &temp_matrix_F_[39]);
  // Vector3Assign(&gyr_input_mat_[6], &temp_matrix_F_[48]);

  temp_matrix_F_[30] = gyr_input_mat_[0];
  temp_matrix_F_[39] = gyr_input_mat_[1];
  temp_matrix_F_[48] = gyr_input_mat_[2];
  temp_matrix_F_[31] = gyr_input_mat_[3];
  temp_matrix_F_[40] = gyr_input_mat_[4];
  temp_matrix_F_[49] = gyr_input_mat_[5];
  temp_matrix_F_[32] = gyr_input_mat_[6];
  temp_matrix_F_[41] = gyr_input_mat_[7];
  temp_matrix_F_[50] = gyr_input_mat_[8];

  temp_matrix_F_[57] = -dt_sec * temp_matrix_2_[0];
  temp_matrix_F_[58] = -dt_sec * temp_matrix_2_[1];
  temp_matrix_F_[59] = -dt_sec * temp_matrix_2_[2];
  temp_matrix_F_[66] = -dt_sec * temp_matrix_2_[3];
  temp_matrix_F_[67] = -dt_sec * temp_matrix_2_[4];
  temp_matrix_F_[68] = -dt_sec * temp_matrix_2_[5];
  temp_matrix_F_[75] = -dt_sec * temp_matrix_2_[6];
  temp_matrix_F_[76] = -dt_sec * temp_matrix_2_[7];
  temp_matrix_F_[77] = -dt_sec * temp_matrix_2_[8];

  MatrixMultiple(temp_matrix_F_, covariance_, temp_matrix_l_1_, 9, 9, 9);
  Transpose(temp_matrix_F_, temp_matrix_F_T_, 9, 9);
  MatrixMultiple(temp_matrix_l_1_, temp_matrix_F_T_, covariance_, 9, 9, 9);

  covariance_[0] += 1e-1;
  covariance_[10] += 1e-1;
  covariance_[20] += 1e-1;
  // initialize rotation covariance
  covariance_[30] += 1e-4;
  covariance_[40] += 1e-4;
  covariance_[50] += 1e-4;
  // initialize velocity covariance
  covariance_[60] += 1e-2;
  covariance_[70] += 1e-2;
  covariance_[80] += 1e-2;

  // LOG_EVERY_N(INFO, 10) << dt_sec << " : " << state_[0] << " " << state_[1] << " " << state_[2]
  //                       << " " << state_[3] << " " << state_[4] << " " << state_[5] << " "
  //                       << state_[6] << " " << state_[7] << " " << state_[8];
  if (CheckStatic(acc_input_)) {
    AddStaticMeasurement(timestamp);
    return 2;
  }
  return 1;
}

static float temp_vec_err_[3];
static float temp_vec_delta_x_[9];
static float temp_matrix_S_[9];
static float temp_matrix_S1_[9];
static float temp_matrix_S2_[9];
static float temp_matrix_PHT_[27];
static float temp_matrix_K_[27];

bool InsertGpsMeasure(int64_t /*timestamp*/, double lat, double lon, size_t num_sta) {
  // convert the gps to utm coordinate
  // transform wgs to utm
  double utm_x;
  double utm_y;
  int utm_zone = utm::LatLonToUTMXY(lat, lon, 0, utm_x, utm_y);
  if (utm_zone_ == 0) {
    utm_zone_ = utm_zone;
    origin_x_ = utm_x;
    origin_y_ = utm_y;

    // reset the position to be zeros
    std::fill_n(state_, 3, 0.0f);
  }

  temp_vec_err_[0] = (utm_x - origin_x_) - state_[0];
  temp_vec_err_[1] = (utm_y - origin_y_) - state_[1];
  temp_vec_err_[2] = -state_[2];

  // insert gps to filter
  temp_matrix_S_[0] = covariance_[0] + 1e-2;
  temp_matrix_S_[1] = covariance_[1];
  temp_matrix_S_[2] = covariance_[2];
  temp_matrix_S_[3] = covariance_[9];
  temp_matrix_S_[4] = covariance_[10] + 1e-2;
  temp_matrix_S_[5] = covariance_[11];
  temp_matrix_S_[6] = covariance_[18];
  temp_matrix_S_[7] = covariance_[19];
  temp_matrix_S_[8] = covariance_[20] + 1e-4;

  // if (cholsl(temp_matrix_S_, temp_matrix_S1_, temp_matrix_S2_, 3)) {
  //   return false;
  // }
  if (!invert3x3Matrix(temp_matrix_S_, temp_matrix_S1_)) {
    InitFilter();
    return false;
  }

  Vector3Assign(&covariance_[0], &temp_matrix_PHT_[0]);
  Vector3Assign(&covariance_[9], &temp_matrix_PHT_[3]);
  Vector3Assign(&covariance_[18], &temp_matrix_PHT_[6]);

  Vector3Assign(&covariance_[27], &temp_matrix_PHT_[9]);
  Vector3Assign(&covariance_[36], &temp_matrix_PHT_[12]);
  Vector3Assign(&covariance_[45], &temp_matrix_PHT_[15]);

  Vector3Assign(&covariance_[54], &temp_matrix_PHT_[18]);
  Vector3Assign(&covariance_[63], &temp_matrix_PHT_[21]);
  Vector3Assign(&covariance_[72], &temp_matrix_PHT_[24]);

  // int arows, int acols, int bcols
  MatrixMultiple(temp_matrix_PHT_, temp_matrix_S1_, temp_matrix_K_, 9, 3, 3);
  MatrixMultiple(temp_matrix_K_, temp_vec_err_, temp_vec_delta_x_, 9, 3, 1);

  if (std::isnan(temp_vec_delta_x_[0])) {
    InitFilter();
    return false;
  }

  // LOG(INFO) << "update : " << temp_vec_delta_x_[0] << " " << temp_vec_delta_x_[1] << " "
  //           << temp_vec_delta_x_[2] << " " << temp_vec_delta_x_[3] << " " << temp_vec_delta_x_[4]
  //           << " " << temp_vec_delta_x_[5] << " " << temp_vec_delta_x_[6] << " "
  //           << temp_vec_delta_x_[7] << " " << temp_vec_delta_x_[8] << " ";
  // if (abs(temp_vec_delta_x_[3]) + abs(temp_vec_delta_x_[4]) + abs(temp_vec_delta_x_[5]) > 2.0) {
  //   LOG(WARNING) << temp_vec_delta_x_[3] << " " << temp_vec_delta_x_[4] << " "
  //                << temp_vec_delta_x_[5] << " ";
  // }

  // constrain the rotation update
  // ConstrainValue(&temp_vec_delta_x_[3], 1.0);
  // ConstrainValue(&temp_vec_delta_x_[4], 1.0);
  // ConstrainValue(&temp_vec_delta_x_[5], 1.0);

  // update state
  VectorAdd(state_, temp_vec_delta_x_, 1.0, state_, 9);

  // update rotation
  ComputeRotationMatrix(&temp_vec_delta_x_[3], temp_matrix_S_);
  MatrixMultiple(state_rotmat_, temp_matrix_S_, temp_matrix_S1_, 3, 3, 3);
  std::memcpy(state_rotmat_, temp_matrix_S1_, sizeof(state_rotmat_));

  // update covariance
  // set (I - KH)
  std::fill_n(temp_matrix_l_1_, 9 * 9, 0.0f);
  SetDiagonalOnes(temp_matrix_l_1_, 9);
  VectorSub(&temp_matrix_l_1_[0], &temp_matrix_K_[0], &temp_matrix_l_1_[0], 3);
  VectorSub(&temp_matrix_l_1_[9], &temp_matrix_K_[3], &temp_matrix_l_1_[9], 3);
  VectorSub(&temp_matrix_l_1_[18], &temp_matrix_K_[6], &temp_matrix_l_1_[18], 3);
  VectorSub(&temp_matrix_l_1_[27], &temp_matrix_K_[9], &temp_matrix_l_1_[27], 3);
  VectorSub(&temp_matrix_l_1_[36], &temp_matrix_K_[12], &temp_matrix_l_1_[36], 3);
  VectorSub(&temp_matrix_l_1_[45], &temp_matrix_K_[15], &temp_matrix_l_1_[45], 3);
  VectorSub(&temp_matrix_l_1_[54], &temp_matrix_K_[18], &temp_matrix_l_1_[54], 3);
  VectorSub(&temp_matrix_l_1_[63], &temp_matrix_K_[21], &temp_matrix_l_1_[63], 3);
  VectorSub(&temp_matrix_l_1_[72], &temp_matrix_K_[24], &temp_matrix_l_1_[72], 3);

  MatrixMultiple(temp_matrix_l_1_, covariance_, temp_matrix_F_T_, 9, 9, 9);
  std::memcpy(covariance_, temp_matrix_F_T_, sizeof(covariance_));

  return true;
}

static float temp_vel_loc_[3];
static float temp_matrix_RT_[9];
static float temp_matrix_H_[27];
static float temp_matrix_HT_[27];
static float temp_matrix_KH_[81];

bool InsertVelocityMeasure(float vx, float vy, float vz, float cov_x, float cov_y, float cov_z) {
  // compute error
  Transpose(state_rotmat_, temp_matrix_RT_, 3, 3);
  MatrixVectorMultiple(temp_matrix_RT_, &state_[6], temp_vel_loc_, 3, 3);
  temp_vec_err_[0] = vx - temp_vel_loc_[0];
  temp_vec_err_[1] = vy - temp_vel_loc_[1];
  temp_vec_err_[2] = vz - temp_vel_loc_[2];

  /////////// compute H matrix  /////////////
  std::fill_n(&temp_matrix_H_[0], 3, 0.0f);
  std::fill_n(&temp_matrix_H_[9], 3, 0.0f);
  std::fill_n(&temp_matrix_H_[18], 3, 0.0f);

  // temp_matrix_S_ -> hat(-RTV)
  float* hat_neg_RTV = temp_matrix_S_;
  // Vector3Negative(temp_vel_loc_);
  Hat(hat_neg_RTV, temp_vel_loc_);
  Vector3Assign(&hat_neg_RTV[0], &temp_matrix_H_[3]);
  Vector3Assign(&hat_neg_RTV[3], &temp_matrix_H_[12]);
  Vector3Assign(&hat_neg_RTV[6], &temp_matrix_H_[21]);

  Vector3Assign(&temp_matrix_RT_[0], &temp_matrix_H_[6]);
  Vector3Assign(&temp_matrix_RT_[3], &temp_matrix_H_[15]);
  Vector3Assign(&temp_matrix_RT_[6], &temp_matrix_H_[24]);

  // compute PHT
  Transpose(temp_matrix_H_, temp_matrix_HT_, 3, 9);
  MatrixMultiple(covariance_, temp_matrix_HT_, temp_matrix_PHT_, 9, 9, 3);

  // compute S
  MatrixMultiple(temp_matrix_H_, temp_matrix_PHT_, temp_matrix_S_, 3, 9, 3);
  temp_matrix_S_[0] += cov_x;
  temp_matrix_S_[4] += cov_y;
  temp_matrix_S_[8] += cov_z;

  // compute K
  if (!invert3x3Matrix(temp_matrix_S_, temp_matrix_S1_)) {
    InitFilter();
    return false;
  }
  MatrixMultiple(temp_matrix_PHT_, temp_matrix_S1_, temp_matrix_K_, 9, 3, 3);

  // compute update
  MatrixMultiple(temp_matrix_K_, temp_vec_err_, temp_vec_delta_x_, 9, 3, 1);

  // update state
  VectorAdd(state_, temp_vec_delta_x_, 1.0, state_, 9);

  // update rotation
  ComputeRotationMatrix(&temp_vec_delta_x_[3], temp_matrix_S_);
  MatrixMultiple(state_rotmat_, temp_matrix_S_, temp_matrix_S1_, 3, 3, 3);
  std::memcpy(state_rotmat_, temp_matrix_S1_, sizeof(state_rotmat_));

  // update covariance
  std::fill_n(temp_matrix_l_1_, 9 * 9, 0.0f);
  SetDiagonalOnes(temp_matrix_l_1_, 9);
  MatrixMultiple(temp_matrix_K_, temp_matrix_H_, temp_matrix_KH_, 9, 3, 9);
  MatrixSub(temp_matrix_l_1_, temp_matrix_KH_, temp_matrix_l_1_, 9, 9);

  MatrixMultiple(temp_matrix_l_1_, covariance_, temp_matrix_F_T_, 9, 9, 9);
  std::memcpy(covariance_, temp_matrix_F_T_, sizeof(covariance_));

  return true;
}

void GetLocalXY(double lat, double lon, float* x, float* y) {
  double utm_x;
  double utm_y;
  int utm_zone = utm::LatLonToUTMXY(lat, lon, 0, utm_x, utm_y);
  *x = utm_x - origin_x_;
  *y = utm_y - origin_y_;
}

void GetPosition(float* result) {
  result[0] = state_[0];
  result[1] = state_[1];
  result[2] = state_[2];
}

void GetVelocity(float* result) {
  result[0] = state_[6];
  result[1] = state_[7];
  result[2] = state_[8];
}

void GetRotationMatrix(float* result) { std::memcpy(result, state_rotmat_, sizeof(state_rotmat_)); }

}  // namespace mobili::ins
