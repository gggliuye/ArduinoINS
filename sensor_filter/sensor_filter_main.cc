// Copyright 2025 DeepMirror Inc. All rights reserved.

#include <filesystem>
#include "common/base/glog.h"
#include "ikfom/common/record/record_reader.h"
#include "ikfom/common/visualizer.h"
#include "ikfom/ins/ins_phone.h"
#include "map/mapping/dataset/make_dataset.h"
#include "map/session/session_factory.h"
#include "map/utils/arcgis_projector.h"

#include "arduino/sensor_filter/sensor_filter.h"

DEFINE_string(session_name, "", "session_name");
DEFINE_string(imu_record, "mobili.ins.imu", "imu name");
DEFINE_string(gnss_record, "mobili.ins.fix", "gnss name");

int main(int argc, char** argv) {
  DM_InitGoogleLogging(argc, argv);

  dm::DmViewer viewer;
  mobili::ins::InitFilter();

  LOG(INFO) << "Load data from : " << FLAGS_session_name;

  auto session = dm::map::GetDefaultSession(FLAGS_session_name);
  dm::ikfom::RecordReader reader(session,
                                 std::vector<std::string>{FLAGS_imu_record, FLAGS_gnss_record});

  dm::ikfom::RecordReader::Message message;
  Eigen::Quaterniond rot_gyr = Eigen::Quaterniond::Identity();
  std::map<int64_t, Sophus::SE3d> base_trajectory;

  bool is_static = false;
  while (!reader.ReachEnd()) {
    double ratio = 100.0 * reader.current_idx() / reader.total_size();
    if (!reader.ReadMessage(&message)) continue;

    if (message.record_name == FLAGS_imu_record) {
      dm::internal::proto::data::ImuData raw_proto;
      raw_proto.ParseFromString(message.content);

      auto acc = raw_proto.linear_acceleration();
      auto gyr = raw_proto.angular_velocity();

      int ret = mobili::ins::InsertImu(message.timestamp, acc.x(), acc.y(), acc.z(), gyr.x(),
                                       gyr.y(), gyr.z());
      // mobili::ins::InsertImu(message.timestamp, acc.x(), acc.y(), acc.z(), 0, 0, gyr.z());
      is_static = (ret == 2);
      {
        static int64_t ts = message.timestamp;
        double dt = static_cast<double>(message.timestamp - ts) * 1e-9;
        ts = message.timestamp;
        rot_gyr =
            rot_gyr *
            Sophus::SO3d::exp(dt * Eigen::Vector3d(gyr.x(), gyr.y(), gyr.z())).unit_quaternion();
      }
    } else if (message.record_name == FLAGS_gnss_record) {
      dm::internal::proto::data::GNSSData raw_proto;
      raw_proto.ParseFromString(message.content);

      auto pos = raw_proto.wgs_position();
      int num_sta = std::stoi(raw_proto.nmea_sentence());
      if (!mobili::ins::InsertGpsMeasure(message.timestamp, pos.y(), pos.x(), num_sta)) break;
      // if (!mobili::ins::InsertGpsMeasure(message.timestamp, pos.y(), pos.x(), num_sta)) {
      //   LOG(ERROR) << "InsertGpsMeasure failed!";
      //   continue;
      // }
      if (!mobili::ins::InsertVelocityMeasure(0, 0, 0, 1e-1, 1e4, 1e-4)) {
        LOG(ERROR) << "InsertVelocityMeasure failed!";
      }

      if (raw_proto.solve_status() !=
          dm::internal::proto::data::GNSSData_SolveStatus_STATUS_NO_FIX) {
        float x, y;
        mobili::ins::GetLocalXY(pos.y(), pos.x(), &x, &y);
        viewer.PublishGps(Eigen::Vector3d(x, y, 0), raw_proto.solve_status());

        //
        // float rotation_matrix[9];
        // mobili::ins::GetRotationMatrix(rotation_matrix);
        // // we are using row major, eigen is column major
        // Eigen::Matrix3f rot = Eigen::Map<Eigen::Matrix3f>(rotation_matrix).transpose();
        viewer.PublishPose(Sophus::SE3d(rot_gyr, Eigen::Vector3d(x, y, 0)), "gyro_rotation");
      }
    } else {
      LOG_EVERY_N(ERROR, 100) << message.record_name << " no call back.";
    }

    if (reader.current_idx() % 10 != 0) continue;

    // Get current pose
    {
      float position[3];
      float velocity[3];
      float rotation_matrix[9];
      mobili::ins::GetPosition(position);
      mobili::ins::GetVelocity(velocity);
      mobili::ins::GetRotationMatrix(rotation_matrix);

      // auto result_map = Eigen::Map<Eigen::MatrixXf>(result.mesh2camera, 4, 4);
      Eigen::Vector3f pos = Eigen::Map<Eigen::Vector3f>(position);
      Eigen::Vector3f vel = Eigen::Map<Eigen::Vector3f>(velocity);
      // we are using row major, eigen is column major
      Eigen::Matrix3f rot = Eigen::Map<Eigen::Matrix3f>(rotation_matrix).transpose();

      Sophus::SE3d pose(Eigen::Quaternionf(rot).cast<double>(), pos.cast<double>());
      viewer.PublishPoseVelocity(pose, vel.cast<double>(), "internal");

      if (is_static) {
        viewer.PublishPoseVelocity(pose, Eigen::Vector3d(0, 0, 4), "static");
      }
    }
    // std::optional<dm::ikfom::PoseData> pose = filter->GetCurrentPose();
    // // add to trajectory
    // if (pose) {
    //   base_trajectory[message.timestamp] = pose->pose_se3();
    // }
    // if (pose) {
    //   pose->position -= filter->origin().value();
    //   // Eigen::Vector3d velocity = Eigen::Vector3d(0, 0, pose->velocity(0));
    //   if (FLAGS_show_zero_velocity && filter->LastZeroVelocityCheck()) {
    //     viewer.PublishPoseVelocity(pose->pose_se3(), Eigen::Vector3d(0, 0, 4), "internal");
    //     continue;
    //   }
    //   // Eigen::Vector3d(0, 0, -pose->velocity.norm() * 3)
    //   viewer.PublishPoseVelocity(pose->pose_se3(),
    //                              Eigen::Vector3d(0, 0, -pose->velocity.norm() * 3), "internal");
    // }
  }

  LOG(INFO) << "Done.";

  viewer.Run();

  LOG(INFO) << "Finished.";
  return 0;
}
