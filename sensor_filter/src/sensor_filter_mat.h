
#pragma once

#include <math.h>
#include <algorithm>  // std::clamp

namespace mobili::ins {

void MatrixVectorMultiple(float* mat, float* vec, float* res, size_t row, size_t col) {
  size_t i, j;
  for (i = 0; i < row; i++) {
    res[i] = 0.0f;
    for (j = 0; j < col; j++) {
      res[i] += mat[i * col + j] * vec[j];
    }
  }
}

void VectorAdd(float* vec_a, float* vec_b, float factor, float* res, size_t dof) {
  size_t i;
  for (i = 0; i < dof; i++) {
    res[i] = vec_a[i] + factor * vec_b[i];
  }
}

void SetDiagonalOnes(float* mat, size_t dof) {
  size_t i = 0;
  size_t idx = 0;
  for (i = 0; i < dof; i++) {
    mat[idx] = 1.0f;
    idx = idx + dof + 1;
  }
}

void VectorSub(float* vec_a, float* vec_b, float* res, size_t dof) {
  size_t i;
  for (i = 0; i < dof; i++) {
    res[i] = vec_a[i] - vec_b[i];
  }
}

void MatrixSub(float* mat_a, float* mat_b, float* mat_res, size_t row, size_t col) {
  size_t i, j;
  for (i = 0; i < row; i++) {
    for (j = 0; j < col; j++) {
      size_t idx = i * col + j;
      mat_res[idx] = mat_a[idx] - mat_b[idx];
    }
  }
}

void Vector3Negative(float* vec) {
  vec[0] = -vec[0];
  vec[1] = -vec[1];
  vec[2] = -vec[2];
}

void VectorAssign(float* vec, float x, float y, float z) {
  vec[0] = x;
  vec[1] = y;
  vec[2] = z;
}

void VectorAssign(float* res, float* vec, size_t dof) {
  size_t i;
  for (i = 0; i < dof; i++) {
    res[i] = vec[i];
  }
}

void Vector3Assign(float* vec, float* res) {
  res[0] = vec[0];
  res[1] = vec[1];
  res[2] = vec[2];
}

void Hat(float* mat, float* vec) {
  mat[0] = 0;
  mat[1] = -vec[2];
  mat[2] = vec[1];

  mat[3] = vec[2];
  mat[4] = 0;
  mat[5] = -vec[0];

  mat[6] = -vec[1];
  mat[7] = vec[0];
  mat[8] = 0;
}

void HatPlusIdentity(float* mat, float* vec) {
  mat[0] = 1;
  mat[1] = -vec[2];
  mat[2] = vec[1];

  mat[3] = vec[2];
  mat[4] = 1;
  mat[5] = -vec[0];

  mat[6] = -vec[1];
  mat[7] = vec[0];
  mat[8] = 1;
}

// C <- A * B
void MatrixMultiple(float* a, float* b, float* c, int arows, int acols, int bcols) {
  int i, j, l;
  for (i = 0; i < arows; ++i) {
    for (j = 0; j < bcols; ++j) {
      c[i * bcols + j] = 0;
      for (l = 0; l < acols; ++l) {
        c[i * bcols + j] += a[i * acols + l] * b[l * bcols + j];
      }
    }
  }
}

void Transpose(float* a, float* at, int m, int n) {
  int i, j;
  for (i = 0; i < m; ++i)
    for (j = 0; j < n; ++j) {
      at[j * m + i] = a[i * n + j];
    }
}

int choldc1(float* a, float* p, int n) {
  int i, j, k;
  float sum;

  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      sum = a[i * n + j];
      for (k = i - 1; k >= 0; k--) {
        sum -= a[i * n + k] * a[j * n + k];
      }
      if (i == j) {
        if (sum <= 0) {
          return 1; /* error */
        }
        p[i] = sqrt(sum);
      } else {
        a[j * n + i] = sum / p[i];
      }
    }
  }

  return 0; /* success */
}

int choldcsl(float* A, float* a, float* p, int n) {
  int i, j, k;
  float sum;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) a[i * n + j] = A[i * n + j];
  if (choldc1(a, p, n)) return 1;
  for (i = 0; i < n; i++) {
    a[i * n + i] = 1 / p[i];
    for (j = i + 1; j < n; j++) {
      sum = 0;
      for (k = i; k < j; k++) {
        sum -= a[j * n + k] * a[k * n + i];
      }
      a[j * n + i] = sum / p[j];
    }
  }

  return 0; /* success */
}

int cholsl(float* A, float* a, float* p, int n) {
  int i, j, k;
  if (choldcsl(A, a, p, n)) return 1;
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      a[i * n + j] = 0.0;
    }
  }
  for (i = 0; i < n; i++) {
    a[i * n + i] *= a[i * n + i];
    for (k = i + 1; k < n; k++) {
      a[i * n + i] += a[k * n + i] * a[k * n + i];
    }
    for (j = i + 1; j < n; j++) {
      for (k = j; k < n; k++) {
        a[i * n + j] += a[k * n + i] * a[k * n + j];
      }
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
      a[i * n + j] = a[j * n + i];
    }
  }

  return 0; /* success */
}

bool invert3x3Matrix(const float m[9], float invOut[9]) {
  // Row-major input
  float det = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) +
              m[2] * (m[3] * m[7] - m[4] * m[6]);

  if (std::fabs(det) < 1e-6f) {
    // Matrix is singular, can't invert
    return false;
  }

  float invDet = 1.0f / det;

  invOut[0] = (m[4] * m[8] - m[5] * m[7]) * invDet;
  invOut[1] = -(m[1] * m[8] - m[2] * m[7]) * invDet;
  invOut[2] = (m[1] * m[5] - m[2] * m[4]) * invDet;

  invOut[3] = -(m[3] * m[8] - m[5] * m[6]) * invDet;
  invOut[4] = (m[0] * m[8] - m[2] * m[6]) * invDet;
  invOut[5] = -(m[0] * m[5] - m[2] * m[3]) * invDet;

  invOut[6] = (m[3] * m[7] - m[4] * m[6]) * invDet;
  invOut[7] = -(m[0] * m[7] - m[1] * m[6]) * invDet;
  invOut[8] = (m[0] * m[4] - m[1] * m[3]) * invDet;

  return true;
}

void ComputeRotationMatrix(float* theta, float* matrix) {
  float angle = std::sqrt(theta[0] * theta[0] + theta[1] * theta[1] + theta[2] * theta[2]);
  if (angle < 1e-6f) {
    HatPlusIdentity(matrix, theta);
  } else {
    // Rodrigues formula
    float axis[3] = {theta[0] / angle, theta[1] / angle, theta[2] / angle};
    float K[9];
    Hat(K, axis);

    float K2[9];
    MatrixMultiple(K, K, K2, 3, 3, 3);
    float sin_a = std::sin(angle);
    float cos_a = std::cos(angle);

    for (int i = 0; i < 9; ++i) {
      float I = (i % 4 == 0 ? 1.0f : 0.0f);
      matrix[i] = I + sin_a * K[i] + (1.0f - cos_a) * K2[i];
    }
  }
}

// Extracts axis-angle from rotation matrix R (row-major)
// Output: axis[3] (unit vector), angle (in radians)
void AxisAngleFromRotationMatrix(const float R[9], float axis[3]) {
  float trace = R[0] + R[4] + R[8];
  float cos_theta = (trace - 1.0f) * 0.5f;
  cos_theta = std::clamp(cos_theta, -1.0f, 1.0f);
  float angle = std::acos(cos_theta);

  // Near zero angle -> identity matrix
  if (angle < 1e-6f) {
    axis[0] = 0.0f;
    axis[1] = 0.0f;
    axis[2] = 0.0f;
    return;
  }

  // Near 180 degrees: special handling
  if (M_PI - angle < 1e-3f) {
    // R is nearly symmetric and trace ~ -1
    float xx = (R[0] + 1.0f) * 0.5f;
    float yy = (R[4] + 1.0f) * 0.5f;
    float zz = (R[8] + 1.0f) * 0.5f;
    float xy = (R[1] + R[3]) * 0.25f;
    float xz = (R[2] + R[6]) * 0.25f;
    float yz = (R[5] + R[7]) * 0.25f;

    if (xx > yy && xx > zz) {
      axis[0] = std::sqrt(xx);
      axis[1] = xy / axis[0];
      axis[2] = xz / axis[0];
    } else if (yy > zz) {
      axis[1] = std::sqrt(yy);
      axis[0] = xy / axis[1];
      axis[2] = yz / axis[1];
    } else {
      axis[2] = std::sqrt(zz);
      axis[0] = xz / axis[2];
      axis[1] = yz / axis[2];
    }
    return;
  }

  // General case
  float sin_theta = std::sin(angle);
  axis[0] = angle * (R[7] - R[5]) / (2.0f * sin_theta);  // (R32 - R23)
  axis[1] = angle * (R[2] - R[6]) / (2.0f * sin_theta);  // (R13 - R31)
  axis[2] = angle * (R[3] - R[1]) / (2.0f * sin_theta);  // (R21 - R12)
  return;
}

void InverseRightJacobian(const float phi[3], float mat[9]) {
  float theta2 = phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2];
  float phi_half[3] = {phi[0] * 0.5f, phi[1] * 0.5f, phi[2] * 0.5f};

  // mat = I + 0.5 * [phi]x
  HatPlusIdentity(mat, phi_half);
  if (theta2 < 1e-8) {
    return;
  }
  float theta = std::sqrt(theta2);
  float phi_scaled[3] = {phi[0] / theta, phi[1] / theta, phi[2] / theta};
  float factor = 1.0f - 0.5f * theta * (1.0f + std::cos(theta)) / std::sin(theta);

  mat[0] += factor * (phi_scaled[0] * phi_scaled[0] - 1);
  mat[1] += factor * phi_scaled[0] * phi_scaled[1];
  mat[2] += factor * phi_scaled[0] * phi_scaled[2];

  mat[3] += factor * phi_scaled[1] * phi_scaled[0];
  mat[4] += factor * (phi_scaled[1] * phi_scaled[1] - 1);
  mat[5] += factor * phi_scaled[1] * phi_scaled[2];

  mat[6] += factor * phi_scaled[2] * phi_scaled[0];
  mat[7] += factor * phi_scaled[2] * phi_scaled[1];
  mat[8] += factor * (phi_scaled[2] * phi_scaled[2] - 1);
}

void ConstrainValue(float* value, float threshold) {
  if (*value < -threshold) {
    *value = -threshold;
  } else if (*value > -threshold) {
    *value = threshold;
  }
}

}  // namespace mobili::ins
