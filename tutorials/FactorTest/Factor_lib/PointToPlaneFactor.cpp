//
// Created by usl on 3/19/21.
//

#include "PointToPlaneFactor.h"

namespace gtsam {

    Vector1 PointToPlaneFactor::computeErrorAndJacobians(const Pose3& wT1, const Pose3& wTm, const Vector3& wVm, const Pose3& Tc,
                                                        OptionalJacobian<1, 6> H1, OptionalJacobian<1, 6> H2,
                                                        OptionalJacobian<1, 3> H3, OptionalJacobian<1, 6> H4) const {
        /// Pre-integrated measurements
        Rot3 deltaR = preintegrated_imu_measurements_.deltaR;
        Vector3 deltaP = preintegrated_imu_measurements_.deltaP;
        Vector3 deltaV = preintegrated_imu_measurements_.deltaV;
        double deltaT = preintegrated_imu_measurements_.deltaT;
        Vector3 gravity = preintegrated_imu_measurements_.gravity;

        Pose3 A = Pose3(Rot3::identity(), wVm*deltaT + 0.5*gravity*deltaT*deltaT);
        Pose3 B = Pose3(deltaR, deltaP);

        Matrix H_1, H_2, H_3, H_4, H_5, H_6, H_7, H_8, H_9, H_10, H_11, H_12;
        Pose3 L1_T_Lmplusi = Tc.inverse(H_1).compose(wT1.inverse(H_2), H_3, H_4).compose(A, H_5, H_6).
                compose(wTm, H_7, H_8).compose(B, H_9, H_10).compose(Tc, H_11, H_12);
        Rot3 L1_R_Lmplusi = L1_T_Lmplusi.rotation();
        Vector3 L1_p_Lmplusi = L1_T_Lmplusi.translation();

        /// Jacobian of L1_T_Lmplusi wrt GTI1
        Matrix H_L1TLmplusi_GTI1 = H_11*H_9*H_7*H_5*H_4*H_2;

        /// Jacobian of L1_T_Lmplusi wrt GTIm
        Matrix H_L1TLmplusi_GTIm = H_11*H_9*H_8;

        /// Jacobian of L1_T_Lmplusi wrt A
        Matrix66 H_L1Tmplusi_A = H_11*H_9*H_7*H_6;
        Matrix33 H_L1Rmplusi_Ra = H_L1Tmplusi_A.block(0, 3, 3,3);
        Matrix33 H_L1pmplusi_pa = H_L1Tmplusi_A.block(3, 3, 3,3);

        /// Jacobian of L1_T_Lmplusi wrt Tc
        Matrix6 H_L1TLmplusi_Tc;
        H_L1TLmplusi_Tc = H_11*H_9*H_7*H_5*H_3*H_1 + H_12;

        /// Jacobian of L1_T_Lmplusi wrt G_v_W
        Rot3 Rc = Tc.rotation();
        Rot3 G_R_I1 = wT1.rotation();
        Matrix63 H_L1TLmplusi_GvIm; ;
        H_L1TLmplusi_GvIm.block(0, 0, 3, 3) = H_L1Rmplusi_Ra;
        H_L1TLmplusi_GvIm.block(3, 0, 3, 3) = H_L1pmplusi_pa*deltaT;

        /// Transform lidar point measurement
        Matrix H_xL1_L1TLmplusi; /// Jacobian of x_L1 wrt L1_T_Lmplusi (3 x 6)
        Point3 x_L1 = L1_T_Lmplusi.transformFrom(lidar_point_measurement_, H_xL1_L1TLmplusi);

        /// Point to plane constraint
        Point3 n_L1 = Point3(plane_param_measurement_.x(), plane_param_measurement_.y(), plane_param_measurement_.z());
        double d_L1 = plane_param_measurement_.z();

        /// Residual
        Vector1 res = Vector1(n_L1.transpose()*x_L1 + d_L1);

        /// Jacobian of res wrt xL1 ( 1 x 3 )
        Matrix H_res_xL1 = n_L1.transpose();

        if (H1)
            (*H1) = (Matrix16() << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GTI1).finished();
        if (H2)
            (*H2) = (Matrix16() << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GTIm).finished();
        if (H3)
            (*H3) = (Matrix13() << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GvIm).finished();
        if (H4)
            (*H4) = (Matrix16() << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_Tc).finished();
        return res;
    }

    Vector PointToPlaneFactor::evaluateError(const Pose3 &wT1, const Pose3 &wTm, const Vector3 &wVm, const Pose3 &Tc,
                                             boost::optional<Matrix &> H1, boost::optional<Matrix &> H2,
                                             boost::optional<Matrix &> H3, boost::optional<Matrix &> H4) const {
        Vector1 error = computeErrorAndJacobians(wT1, wTm, wVm, Tc, H1, H2, H3, H4);
        return error;
    }
}