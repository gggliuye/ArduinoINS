package(default_visibility = ["//visibility:public"])

cc_library(
    name = "sensor_filter",
    hdrs = [
        "sensor_filter_mat.h",
        "sensor_filter.h",
    ],
    srcs = ["sensor_filter.cc"],
    deps = [
        "@dm_core_cpp//common/base:glog",
        "@dm_core_cpp//math:sophus_utils",
        "@dm_core_cpp//third_party/utm",
    ],
)

cc_binary(
    name = "sensor_filter_main",
    srcs = ["sensor_filter_main.cc",],
    deps = [
        ":sensor_filter",
        "//ikfom/common:visualizer_opencv",
        "//ikfom/common:visualizer",
        "//ikfom/common/record:record_reader",
        "//ikfom/ins:ins_phone",
        "//map/mapping/dataset:make_dataset",
        "@dm_core_map//map/session:sensor_utils",
        "@dm_core_map//map/session:session_factory",
        "@dm_core_map//map/utils:arcgis_projector",
        "@opencv",
    ],
)
