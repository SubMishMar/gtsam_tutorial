//
// Created by usl on 3/19/21.
//

#include "PointToPlaneFactor.h"

namespace lin_estimator {

    Vector1 PointToPlaneFactor::computeErrorAndJacobians(const Pose3& wT1,
                                                         const Pose3& wTm,
                                                         const Vector3& wVm,
                                                         const Vector6& Bm,
                                                         const Pose3& Tc,
                                                         OptionalJacobian<1, 6> H1,
                                                         OptionalJacobian<1, 6> H2,
                                                         OptionalJacobian<1, 3> H3,
                                                         OptionalJacobian<1, 6> H4,
                                                         OptionalJacobian<1, 6> H5) const {
        /// Pre-integrated measurements
        Rot3 deltaR = preintegrated_imu_measurements_.deltaR;
        Vector3 deltaP = preintegrated_imu_measurements_.deltaP;
        Vector3 deltaV = preintegrated_imu_measurements_.deltaV;
        double deltaT = preintegrated_imu_measurements_.deltaT;
        Vector3 gravity = preintegrated_imu_measurements_.gravity;
        Matrix93 H_bias_omega = preintegrated_imu_measurements_.H_bias_omega;
        Matrix93 H_bias_accel = preintegrated_imu_measurements_.H_bias_accel;

        Matrix33 H_delR_bias_accel = H_bias_accel.block(0, 0, 3, 3);
        Matrix33 H_delP_bias_accel = H_bias_accel.block(3, 0, 3, 3);
        Matrix33 H_delV_bias_accel = H_bias_accel.block(6, 0, 3, 3);
        Vector3 bias_accel_hat = Vector3(Bm(0), Bm(1), Bm(2));

        Matrix33 H_delR_bias_omega = H_bias_omega.block(0, 0, 3, 3);
        Matrix33 H_delP_bias_omega = H_bias_omega.block(3, 0, 3, 3);
        Matrix33 H_delV_bias_omega = H_bias_omega.block(6, 0, 3, 3);
        Vector3 bias_omega_hat = Vector3(Bm(3), Bm(4), Bm(5));

        Pose3 A = Pose3(Rot3::identity(), wVm*deltaT + 0.5*gravity*deltaT*deltaT);
        Pose3 B11 = Pose3(Rot3::identity(), H_delP_bias_omega*bias_omega_hat);
        Pose3 B12 = Pose3(Rot3::identity(), H_delP_bias_accel*bias_accel_hat);
        Pose3 B2 = Pose3(deltaR, deltaP);
        Matrix3 H_exp_Hbg;
        Pose3 B3 = Pose3(gtsam::Rot3::Expmap(H_delR_bias_omega*bias_omega_hat, H_exp_Hbg), gtsam::Vector3::Zero());
        Matrix H_1, H_2, H_3, H_4, H_5, H_6, H_7, H_8, H_9, H_10, H_11, H_12, H_13, H_14, H_15, H_16, H_17, H_18;
        Pose3 L1_T_Lmplusi = Tc.inverse(H_1).
                            compose(wT1.inverse(H_2), H_3, H_4).
                            compose(A, H_5, H_6).
                            compose(wTm, H_7, H_8).
                            compose(B11, H_9, H_10).
                            compose(B12, H_11, H_12).
                            compose(B2, H_13, H_14).
                            compose(B3, H_15, H_16).
                            compose(Tc, H_17, H_18);

        /// Jacobian of L1_T_Lmplusi wrt GTI1
        Matrix H_L1TLmplusi_GTI1 = H_17*H_15*H_13*H_11*H_9*H_7*H_5*H_4*H_2;

        /// Jacobian of L1_T_Lmplusi wrt GTIm
        Matrix H_L1TLmplusi_GTIm = H_17*H_15*H_13*H_11*H_9*H_8;

        /// Jacobian of L1_T_Lmplusi wrt A
        Matrix66 H_L1Tmplusi_A = H_17*H_15*H_13*H_11*H_9*H_7*H_6;
        Matrix33 H_L1Rmplusi_pa = H_L1Tmplusi_A.block(0, 3, 3,3);
        Matrix33 H_L1pmplusi_pa = H_L1Tmplusi_A.block(3, 3, 3,3);

        /// Jacobian of L1_T_Lmplusi wrt Tc
        Matrix6 H_L1TLmplusi_Tc;
        H_L1TLmplusi_Tc = H_18 + H_17*H_15*H_13*H_11*H_9*H_7*H_5*H_3*H_1;

        /// Jacobian of L1_T_Lmplusi wrt G_v_W
        Matrix63 H_L1TLmplusi_GvIm; ;
        H_L1TLmplusi_GvIm.block(0, 0, 3, 3) = H_L1Rmplusi_pa;
        H_L1TLmplusi_GvIm.block(3, 0, 3, 3) = H_L1pmplusi_pa*deltaT;

        /// Jacobian of L1_T_Lmplusi wrt B12
        Matrix66 H_L1TLmplusi_B12 = H_17*H_15*H_13*H_12;
        Matrix33 H_L1RLmplusi_pB12 = H_L1TLmplusi_B12.block(0, 3, 3, 3);
        Matrix33 H_L1pLmplusi_pB12 = H_L1TLmplusi_B12.block(3, 3, 3, 3);

        /// Jacobian of L1_T_Lmplusi wrt ba_hat
        Matrix63 H_L1TLmplusi_ba_hat;
        H_L1TLmplusi_ba_hat.block(0, 0, 3, 3) = H_L1RLmplusi_pB12;
        H_L1TLmplusi_ba_hat.block(3, 0, 3, 3) = H_L1pLmplusi_pB12*H_delP_bias_accel;

        /// Jacobian of L1_T_Lmplusi wrt B11
        Matrix66 H_L1TLmplusi_B11 = H_17*H_15*H_13*H_11*H_10;
        Matrix33 H_L1RLmplusi_pB11 = H_L1TLmplusi_B11.block(0, 0, 3, 3);
        Matrix33 H_L1pLmplusi_pB11 = H_L1TLmplusi_B11.block(3, 3, 3, 3);
        /// Jacobian 1 of L1_T_Lmplusi wrt bg_hat
        Matrix63 H_L1TLmplusi_bg_hat1;
        H_L1TLmplusi_bg_hat1.block(0, 0, 3, 3) = H_L1RLmplusi_pB11;
        H_L1TLmplusi_bg_hat1.block(3, 0, 3, 3) = H_L1pLmplusi_pB11*H_delP_bias_omega;

        /// Jacobian of L1_T_Lmplusi wrt B3
        Matrix66 H_L1TLmplusi_B3 = H_17*H_16;
        Matrix63 H_L1TLmplusi_RB3 = H_L1TLmplusi_B3.block(0, 0, 6, 3);
        /// Jacobian 2 of L1_T_Lmplusi wrt bg_hat
        Matrix63 H_L1TLmplusi_bg_hat2;
        H_L1TLmplusi_bg_hat2 = H_L1TLmplusi_RB3*H_exp_Hbg*H_delR_bias_omega;

        /// Jacobian of L1_T_Lmplusi wrt bg_hat
        Matrix63 H_L1TLmplusi_bg_hat = H_L1TLmplusi_bg_hat1 + H_L1TLmplusi_bg_hat2;

        /// Jacobian of L1_T_Lmplusi wrt b = [ba, bg]
        Matrix6 H_L1TLmplusi_b;
        H_L1TLmplusi_b.block(0, 0, 6, 3) = H_L1TLmplusi_ba_hat;
        H_L1TLmplusi_b.block(0, 3, 6, 3) = H_L1TLmplusi_bg_hat;

        /// Transform lidar point measurement
        Matrix H_xL1_L1TLmplusi; /// Jacobian of x_L1 wrt L1_T_Lmplusi (3 x 6)
        Point3 x_L1 = L1_T_Lmplusi.transformFrom(lidar_point_measurement_, H_xL1_L1TLmplusi);

        /// Point to plane constraint
        Point3 n_L1 = Point3(plane_param_measurement_.x(), plane_param_measurement_.y(), plane_param_measurement_.z());
        double d_L1 = plane_param_measurement_.z();

        /// Residual
        Vector1 res = Vector1(weight_*(n_L1.transpose()*x_L1 + d_L1));

        /// Jacobian of res wrt xL1 ( 1 x 3 )
        Matrix H_res_xL1 = weight_*n_L1.transpose();

        if (H1)
            (*H1) = (Matrix16() << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GTI1).finished();
        if (H2)
            (*H2) = (Matrix16() << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GTIm).finished();
        if (H3)
            (*H3) = (Matrix13() << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GvIm).finished();
        if (H4)
            (*H4) = (Matrix16() << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_b).finished();
        if (H5)
            (*H5) = (Matrix16() << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_Tc).finished();
        return res;
    }

    Vector PointToPlaneFactor::evaluateError(const Pose3& wT1,
                                             const Pose3& wTm,
                                             const Vector3& wVm,
                                             const Vector6& Bm,
                                             const Pose3& Tc,
                                             boost::optional<Matrix&> H1,
                                             boost::optional<Matrix&> H2,
                                             boost::optional<Matrix&> H3,
                                             boost::optional<Matrix&> H4,
                                             boost::optional<Matrix&> H5) const {
        Vector1 error = computeErrorAndJacobians(wT1, wTm, wVm, Bm, Tc, H1, H2, H3, H4, H5);
        return error;
    }
}