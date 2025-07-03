## run with esp32

```
bazel run -c opt //map/tools:imu_gps_esp32_to_record_main -- \
-session_path=/LidarMapping/data/esp32/mobili_7


SESSION_NAME=mobili_7
bazel run -c opt //ikfom/ins/processor:phone_fusion_main -- \
-map_storage_input_directories=/LidarMapping/data/esp32 \
-session_name=${SESSION_NAME} \
-show_connection=true \
-imu_record=mobili.imu \
-gnss_record=mobili.fix


bazel run -c opt //arduino/sensor_filter:sensor_filter_main -- \
-map_storage_input_directories=/LidarMapping/data/esp32 \
-session_name=${SESSION_NAME} \
-imu_record=mobili.imu \
-gnss_record=mobili.fix
```
