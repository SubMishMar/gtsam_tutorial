/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file HandEyeCalibration AX = XB
 * @brief Example of application of ISAM2 for HandEyeCalibration AX = XB
 */

#include <boost/program_options.hpp>

// GTSAM related includes.


#include <gtsam/navigation/ImuFactor.h>


#include <cstring>
#include <fstream>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <ceres/autodiff_cost_function.h>
#include <ceres/numeric_diff_cost_function.h>

class rotationError {
private:
    // All the given info
    const Eigen::Vector3d axis_cam_;
    const Eigen::Vector3d axis_lidar_;

public:
    rotationError(const Eigen::Vector3d axis_cam,
                  const Eigen::Vector3d axis_lidar):
            axis_cam_(axis_cam), axis_lidar_(axis_lidar)
    {}

    template  <typename  T>
    bool operator() (const T* const calib_angle,
                     T* residual) const {
        Eigen::Matrix<T, 3, 1> _axis_cam;
        _axis_cam(0) = T(axis_cam_.x());
        _axis_cam(1) = T(axis_cam_.y());
        _axis_cam(2) = T(axis_cam_.z());
        Eigen::Matrix<T, 3, 1> _axis_lidar;
        _axis_lidar(0) = T(axis_lidar_.x());
        _axis_lidar(1) = T(axis_lidar_.y());
        _axis_lidar(2) = T(axis_lidar_.z());

        Eigen::Matrix<T, 3, 1> _axis_lidar_transformed;
        ceres::AngleAxisRotatePoint(calib_angle, _axis_lidar.data(), _axis_lidar_transformed.data());

        residual[0] = _axis_lidar_transformed(0) - _axis_cam(0);
        residual[1] = _axis_lidar_transformed(1) - _axis_cam(1);
        residual[2] = _axis_lidar_transformed(2) - _axis_cam(2);
        return true;
    }
};

using namespace std;
using namespace gtsam;

namespace po = boost::program_options;

struct IMUCalibration {
    double accelerometer_sigma;
    double gyroscope_sigma;
    double integration_sigma;
    double accelerometer_bias_sigma;
    double gyroscope_bias_sigma;
    double average_delta_t;
};

struct ImuMeasurement {
    uint32_t time_sec;
    uint32_t time_nsec;
    double timestamp; // time_sec + time_nsec/1e9
    double dt;
    Vector3 accelerometer;
    Vector3 gyroscope;  // omega
};

struct LOPoseMeasurement {
    uint32_t time_sec;
    uint32_t time_nsec;
    double timestamp; // time_sec + time_nsec/1e9
    Vector3 position;  // x,y,z
    Quaternion quat; // qw, qx, qy, qz
};

const string output_filename = "IMULOPoseResults.csv";

string imu_meta_data_filename;
string imu_data_filename;
string lopose_data_filename;

void loadIMULOPoseData(IMUCalibration& imu_calibration,
                       vector<ImuMeasurement>& imu_measurements,
                       vector<LOPoseMeasurement>& lopose_measurements) {
    string line;
    // Read IMU metadata and compute relative sensor pose transforms
    // AccelerometerBiasSigma GyroscopeBiasSigma AverageDeltaT
    ifstream imu_metadata(imu_meta_data_filename.c_str());
    std::cout <<"-- Reading sensor metadata" << std::endl;
    getline(imu_metadata, line, '\n');  // ignore the first line
//    // Load IMU calibration
    getline(imu_metadata, line, '\n');
    sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf",
            &imu_calibration.accelerometer_sigma, &imu_calibration.gyroscope_sigma, &imu_calibration.integration_sigma,
            &imu_calibration.accelerometer_bias_sigma, &imu_calibration.gyroscope_bias_sigma,
            &imu_calibration.average_delta_t);
    printf("%lf %lf %lf %lf %lf %lf\n",
           imu_calibration.accelerometer_sigma, imu_calibration.gyroscope_sigma, imu_calibration.integration_sigma,
           imu_calibration.accelerometer_bias_sigma, imu_calibration.gyroscope_bias_sigma,
           imu_calibration.average_delta_t);
    // Read IMU data
    // Time dt accelX accelY accelZ omegaX omegaY omegaZ
    printf("-- Reading IMU measurements from file\n");{
        ifstream imu_data(imu_data_filename.c_str());
        uint32_t time_sec = 0, time_nsec = 0;
        double acc_x = 0, acc_y = 0, acc_z = 0, gyro_x = 0, gyro_y = 0, gyro_z = 0;
        int count = 0;
        double timestamp_prev = 0;
        std::cout.precision(20);
        while (!imu_data.eof()) {
            getline(imu_data, line, '\n');
            sscanf(line.c_str(), "%u, %u, %lf, %lf, %lf, %lf, %lf, %lf",
                   &time_sec, &time_nsec, &acc_x, &acc_y, &acc_z, &gyro_x, &gyro_y, &gyro_z);
            ImuMeasurement measurement;
            measurement.time_sec = time_sec;
            measurement.time_nsec = time_nsec;
            measurement.timestamp = double(time_sec) + double(time_nsec)/1e9;
            if(count == 0) {
                measurement.dt = measurement.timestamp;
            }
            else {
                measurement.dt = measurement.timestamp - timestamp_prev;
            }
            measurement.accelerometer = Vector3(acc_x, acc_y, acc_z);
            measurement.gyroscope = Vector3(gyro_x, gyro_y, gyro_z);
            imu_measurements.push_back(measurement);
            count++;
            timestamp_prev = measurement.timestamp;
        }
    }

    // Read Lopose data
    // TimeSec,TimeNSec, X,Y,Z, qx, qy, qz, qw
    printf("-- Reading LOPose measurements from file\n"); {
        ifstream lopose_data(lopose_data_filename.c_str());
//        getline(lopose_data, line, '\n');  // ignore the first line

        uint32_t time_sec = 0, time_nsec = 0;
        double lopose_x = 0, lopose_y = 0, lopose_z = 0;
        double loquat_x = 0, loquat_y = 0, loquat_z = 0, loquat_w = 0;
        while (!lopose_data.eof()) {
            getline(lopose_data, line, '\n');
            sscanf(line.c_str(), "%u, %u, "
                                 "%lf,%lf,%lf, "
                                 "%lf,%lf,%lf,%lf",
                    &time_sec, &time_nsec,
                    &lopose_x, &lopose_y, &lopose_z,
                    &loquat_x, &loquat_y, &loquat_z, &loquat_w);
            LOPoseMeasurement measurement;
            measurement.time_sec = time_sec;
            measurement.time_nsec = time_nsec;
            measurement.timestamp = double(time_sec) + double(time_nsec)/1e9;
            measurement.position = Vector3(lopose_x, lopose_y, lopose_z);
            measurement.quat = Quaternion(loquat_w, loquat_x, loquat_y, loquat_z);
            lopose_measurements.push_back(measurement);
        }
        return;
    }
}

int main(int argc, char* argv[]) {
    imu_meta_data_filename = argv[1];
    imu_data_filename = argv[2];
    lopose_data_filename = argv[3];

    std::cout << "imu_meta_data_filename: " << imu_meta_data_filename << std::endl;
    std::cout << "imu_data_filename: " << imu_data_filename << std::endl;
    std::cout << "gps_data_filename: " << lopose_data_filename << std::endl;

    IMUCalibration imu_calibration;
    vector<ImuMeasurement> imu_measurements;
    vector<LOPoseMeasurement> lopose_measurements;
    loadIMULOPoseData(imu_calibration, imu_measurements, lopose_measurements);

    // Configure different variables
    double g = 9.8;
    auto w_coriolis = Vector3::Zero();  // zero vector

//    // Configure noise models
    auto noise_model_lopose = noiseModel::Diagonal::Precisions((Vector6()
                                                                       << Vector3::Constant(1),
                                                                          Vector3::Constant(1)).finished());
//
    // Set initial conditions for the estimated trajectory
    // initial pose is the reference frame (navigation frame)
    // the system is stationary at the beginning
    auto current_pose_global = Pose3(Rot3(), Vector3::Zero());
    Vector3 current_velocity_global = Vector3::Zero();
    auto current_bias = imuBias::ConstantBias();  // init with zero bias
//
    auto sigma_init_x = noiseModel::Diagonal::Precisions((Vector6() << Vector3::Constant(0),
            Vector3::Constant(1.0)).finished());
    auto sigma_init_v = noiseModel::Diagonal::Sigmas(Vector3::Constant(1000.0));
    auto sigma_init_b = noiseModel::Diagonal::Sigmas((Vector6() << Vector3::Constant(0.100),
                                                                   Vector3::Constant(5.00e-05))
                                                     .finished());
//
    // Set IMU preintegration parameters
    Matrix33 measured_acc_cov = I_3x3 * pow(imu_calibration.accelerometer_sigma, 2);
    Matrix33 measured_omega_cov = I_3x3 * pow(imu_calibration.gyroscope_sigma, 2);
    // error committed in integrating position from velocities
    Matrix33 integration_error_cov = I_3x3 * pow(imu_calibration.integration_sigma, 2);
//
    auto imu_params = PreintegratedImuMeasurements::Params::MakeSharedU(g);
    imu_params->accelerometerCovariance = measured_acc_cov;     // acc white noise in continuous
    imu_params->integrationCovariance = integration_error_cov;  // integration uncertainty continuous
    imu_params->gyroscopeCovariance = measured_omega_cov;       // gyro white noise in continuous
    imu_params->omegaCoriolis = w_coriolis;

    std::shared_ptr<PreintegratedImuMeasurements> current_summarized_measurement = nullptr;

    for (size_t i = 1; i < lopose_measurements.size(); i++) {
        double timestamp_curr = lopose_measurements[i].timestamp;
        double timestamp_prev = lopose_measurements[i-1].timestamp;
        // Summarize IMU data between the previous GPS measurement and now
        current_summarized_measurement = std::make_shared<PreintegratedImuMeasurements>(imu_params, current_bias);
        static size_t included_imu_measurement_count = 0;
        static size_t j = 0;
        while (j < imu_measurements.size() && imu_measurements[j].timestamp <= timestamp_curr) {
            if (imu_measurements[j].timestamp >= timestamp_prev) {
                current_summarized_measurement->integrateMeasurement(imu_measurements[j].accelerometer,
                                                                     imu_measurements[j].gyroscope,
                                                                     imu_measurements[j].dt);
                included_imu_measurement_count++;
            }
            j++;
        }
        auto lopose_prev = Pose3(lopose_measurements[i-1].quat, lopose_measurements[i-1].position);
        auto lopose_curr = Pose3(lopose_measurements[i].quat, lopose_measurements[i].position);
        auto delta_lopose = lopose_prev.inverse()*lopose_curr;
    }
}

