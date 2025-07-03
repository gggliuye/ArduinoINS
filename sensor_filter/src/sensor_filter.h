
#pragma once
#include <Arduino.h>
#include <math.h>
#include <cstring>
// #include <iostream>
// #include "common/base/glog.h"

namespace mobili::ins {

void InitFilter();
void InitBuffers();

// -1 : error, 0: init, 1 : normal, 2: static
int InsertImu(int64_t timestamp, float ax, float ay, float az, float gx, float gy, float gz);
bool InsertGpsMeasure(int64_t timestamp, double lat, double lon, size_t num_sta);
bool InsertVelocityMeasure(float vx, float vy, float vz, float cov_x, float cov_y, float cov_z);

void GetLocalXY(double lat, double lon, float* x, float* y);
void GetPosition(float* result);
void GetVelocity(float* result);
void GetRotationMatrix(float* result);

}  // namespace mobili::ins
