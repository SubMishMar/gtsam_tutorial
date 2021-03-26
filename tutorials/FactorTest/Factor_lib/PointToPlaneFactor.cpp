//
// Created by usl on 3/19/21.
//

#include "PointToPlaneFactor.h"

namespace gtsam {

    Vector1 PointToPlaneFactor::computeErrorAndJacobians1(const Pose3& wT1, const Pose3& wTm, const Vector3& wVm, const Pose3& Tc,
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

        /// Jacobian of L1_T_Lmplusi wrt Tc (TODO: Figure out how to use GTSAM's Jacobians)
        Matrix3 H_L1RLmplusi_Rc = Matrix3::Identity() - L1_R_Lmplusi.matrix().transpose();
        Matrix3 H_L1RLmplusi_pc = Matrix3::Zero();
        Matrix3 H_L1pLmplusi_Rc = skewSymmetric(L1_p_Lmplusi);
        Matrix3 H_L1pLmplusi_pc = L1_R_Lmplusi.matrix() - Matrix3::Identity();
        Matrix6 H_L1TLmplusi_Tc;
        H_L1TLmplusi_Tc.block(0, 0, 3, 3) = H_L1RLmplusi_Rc;
        H_L1TLmplusi_Tc.block(0, 3, 3, 3) = H_L1RLmplusi_pc;
        H_L1TLmplusi_Tc.block(3, 0, 3, 3) = H_L1pLmplusi_Rc;
        H_L1TLmplusi_Tc.block(3, 3, 3, 3) = H_L1pLmplusi_pc;

        /// Jacobian of L1_T_Lmplusi wrt G_v_W (TODO: Figure out how to use GTSAM's Jacobians)
        Rot3 Rc = Tc.rotation();
        Rot3 G_R_I1 = wT1.rotation();
        Matrix3 H_L1RLmplusi_GvIm = Matrix3::Zero();
        Matrix3 H_L1pLmplusi_GvIm = Rc.matrix().transpose()*G_R_I1.matrix().transpose()*deltaT;
        Matrix63 H_L1TLmplusi_GvIm;
        H_L1TLmplusi_GvIm.block(0, 0, 3, 3) = H_L1RLmplusi_GvIm;
        H_L1TLmplusi_GvIm.block(3, 0, 3, 3) = H_L1pLmplusi_GvIm;

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
            (*H1) = (Matrix(1, 6) << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GTI1).finished();
        if (H2)
            (*H2) = (Matrix(1, 6) << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GTIm).finished();
        if (H3)
            (*H3) = (Matrix(1, 3) << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GvIm).finished();
        if (H4)
            (*H4) = (Matrix(1, 6) << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_Tc).finished();
        return res;
    }

    Vector1 PointToPlaneFactor::computeErrorAndJacobians2(const Pose3& wT1, const Pose3& wTm, const Vector3& wVm, const Pose3& Tc,
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

        Pose3 T1 = Tc.inverse(H_1);
        Pose3 T2 = wT1.inverse(H_2);
        Pose3 T3 = T1.compose(T2, H_3, H_4);
        Pose3 T4 = T3.compose(A, H_5, H_6);
        Pose3 T5 = T4.compose(wTm, H_7, H_8);
        Pose3 T6 = T5.compose(B, H_9, H_10);
        Pose3 T7 = T6.compose(Tc, H_11, H_12);

        Pose3 L1_T_Lmplusi = T7;
        Rot3 L1_R_Lmplusi = L1_T_Lmplusi.rotation();
        Vector3 L1_p_Lmplusi = L1_T_Lmplusi.translation();

        /// Jacobian of L1_T_Lmplusi wrt GTI1
        Matrix H_L1TLmplusi_GTI1 = H_11*H_9*H_7*H_5*H_4*H_2;
//
        /// Jacobian of L1_T_Lmplusi wrt GTIm
        Matrix H_L1TLmplusi_GTIm = H_11*H_9*H_8;

        /// Jacobian of L1_T_Lmplusi wrt Tc (TODO: Figure out how to use GTSAM's Jacobians)
        Matrix3 H_L1RLmplusi_Rc = Matrix3::Identity() - L1_R_Lmplusi.matrix().transpose();
        Matrix3 H_L1RLmplusi_pc = Matrix3::Zero();
        Matrix3 H_L1pLmplusi_Rc = skewSymmetric(L1_p_Lmplusi);
        Matrix3 H_L1pLmplusi_pc = L1_R_Lmplusi.matrix() - Matrix3::Identity();
        Matrix6 H_L1TLmplusi_Tc;
        H_L1TLmplusi_Tc.block(0, 0, 3, 3) = H_L1RLmplusi_Rc;
        H_L1TLmplusi_Tc.block(0, 3, 3, 3) = H_L1RLmplusi_pc;
        H_L1TLmplusi_Tc.block(3, 0, 3, 3) = H_L1pLmplusi_Rc;
        H_L1TLmplusi_Tc.block(3, 3, 3, 3) = H_L1pLmplusi_pc;

        /// Jacobian of L1_T_Lmplusi wrt G_v_W (TODO: Figure out how to use GTSAM's Jacobians)
        Rot3 Rc = Tc.rotation();
        Rot3 G_R_I1 = wT1.rotation();
        Matrix3 H_L1RLmplusi_GvIm = Matrix3::Zero();
        Matrix3 H_L1pLmplusi_GvIm = Rc.matrix().transpose()*G_R_I1.matrix().transpose()*deltaT;
        Matrix63 H_L1TLmplusi_GvIm;
        H_L1TLmplusi_GvIm.block(0, 0, 3, 3) = H_L1RLmplusi_GvIm;
        H_L1TLmplusi_GvIm.block(3, 0, 3, 3) = H_L1pLmplusi_GvIm;

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
            (*H1) = (Matrix(1, 6) << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GTI1).finished();
        if (H2)
            (*H2) = (Matrix(1, 6) << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GTIm).finished();
        if (H3)
            (*H3) = (Matrix(1, 3) << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GvIm).finished();
        if (H4)
            (*H4) = (Matrix(1, 6) << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_Tc).finished();
        return res;
    }

    Vector1 PointToPlaneFactor::computeErrorAndJacobians3(const Pose3& wT1, const Pose3& wTm, const Vector3& wVm, const Pose3& Tc,
                                                          OptionalJacobian<1, 6> H1, OptionalJacobian<1, 6> H2,
                                                          OptionalJacobian<1, 3> H3, OptionalJacobian<1, 6> H4) const {
        /// Pre-integrated measurements
        Rot3 deltaR = preintegrated_imu_measurements_.deltaR;
        Vector3 deltaP = preintegrated_imu_measurements_.deltaP;
        Vector3 deltaV = preintegrated_imu_measurements_.deltaV;
        double deltaT = preintegrated_imu_measurements_.deltaT;
        Vector3 gravity = preintegrated_imu_measurements_.gravity;

        /// Calibration params
        Rot3 Rc = Tc.rotation();
        Vector3 pc = Tc.translation();

        /// First State
        Rot3 G_R_I1 = wT1.rotation();
        Vector3 G_p_I1 = wT1.translation();

        /// Current State
        Rot3 G_R_Im = wTm.rotation();
        Vector3 G_p_Im = wTm.translation();

        /// Current State plus
        Rot3 G_R_Implusi = G_R_Im.compose(deltaR);
        Vector3 G_p_Implusi = G_p_Im + wVm*deltaT + 0.5*gravity*deltaT*deltaT + G_R_Im*deltaP;
        Pose3 G_T_Implusi = Pose3(G_R_Implusi, G_p_Implusi);

        /// L1_T_Lmplusi
        Pose3 G_T_L1 = wT1.compose(Tc);
        Pose3 G_T_Lmplusi = G_T_Implusi.compose(Tc);
        Pose3 L1_T_Lmplusi = G_T_L1.between(G_T_Implusi);
        Rot3 L1_R_Lmplusi = L1_T_Lmplusi.rotation();
        Vector3 L1_p_Lmplusi = L1_T_Lmplusi.translation();

        /// Transform lidar point
        Matrix H_xL1_L1TLmplusi; /// Jacobian of x_L1 wrt L1_T_Lmplusi (3 x 6)
        Point3 x_L1 = L1_T_Lmplusi.transformFrom(lidar_point_measurement_, H_xL1_L1TLmplusi);

        /// Point to plane constraint
        Point3 n_L1 = Point3(plane_param_measurement_.x(), plane_param_measurement_.y(), plane_param_measurement_.z());
        double d_L1 = plane_param_measurement_.z();

        /// Residual
        Vector1 res = Vector1(n_L1.transpose()*x_L1 + d_L1);
        Matrix H_res_xL1 = n_L1.transpose(); /// Jacobian of res wrt xL1 ( 1 x 3 )

        /// Jacobians of L1_T_Lmplus wrt G_T_I1
        Matrix3 H_L1RLmplusi_GRI1 = -Rc.matrix().transpose()*G_R_Implusi.matrix().transpose()*G_R_I1.matrix();
        Matrix3 H_L1RLmplusi_GpI1 = Matrix3::Zero();
        Matrix3 H_L1pLmplusi_GRI1 =  Rc.matrix().transpose()*skewSymmetric(G_R_I1.matrix().transpose()*(G_R_Implusi.matrix().transpose()*pc + G_p_Implusi - G_p_I1));
        Matrix3 H_L1pLmplusi_GpI1 = -Rc.matrix().transpose();

        Matrix6 H_L1TLmplusi_GTI1;
        H_L1TLmplusi_GTI1.block(0, 0, 3, 3) = H_L1RLmplusi_GRI1;
        H_L1TLmplusi_GTI1.block(0, 3, 3, 3) = H_L1RLmplusi_GpI1;
        H_L1TLmplusi_GTI1.block(3, 0, 3, 3) = H_L1pLmplusi_GRI1;
        H_L1TLmplusi_GTI1.block(3, 3, 3, 3) = H_L1pLmplusi_GpI1;

        /// Jacobians of L1_T_Lmplus wrt G_T_Im
        Matrix3 H_L1RLmplusi_GRIm = Rc.matrix().transpose()*deltaR.matrix();
        Matrix3 H_L1RLmplusi_GpIm = Matrix3::Zero();
        Matrix3 H_L1pLmplusi_GRIm = -Rc.matrix().transpose()*G_R_I1.matrix().transpose()*G_R_Im.matrix()*skewSymmetric(deltaR.matrix()*pc + deltaP);
        Matrix3 H_L1pLmplusi_GpIm = Rc.matrix().transpose()*G_R_I1.matrix().transpose()*G_R_Im.transpose();

        Matrix6 H_L1TLmplusi_GTIm;
        H_L1TLmplusi_GTIm.block(0, 0, 3, 3) = H_L1RLmplusi_GRIm;
        H_L1TLmplusi_GTIm.block(0, 3, 3, 3) = H_L1RLmplusi_GpIm;
        H_L1TLmplusi_GTIm.block(3, 0, 3, 3) = H_L1pLmplusi_GRIm;
        H_L1TLmplusi_GTIm.block(3, 3, 3, 3) = H_L1pLmplusi_GpIm;

        /// Jacobians of L1_T_Lmplus wrt G_v_Im
        Matrix3 H_L1RLmplusi_GvIm = Matrix3::Zero();
        Matrix3 H_L1pLmplusi_GvIm = Rc.matrix().transpose()*G_R_I1.matrix().transpose()*deltaT;

        Matrix63 H_L1TLmplusi_GvIm;
        H_L1TLmplusi_GvIm.block(0, 0, 3, 3) = H_L1RLmplusi_GvIm;
        H_L1TLmplusi_GvIm.block(3, 0, 3, 3) = H_L1pLmplusi_GvIm;

        /// Jacobians of L1_T_Lmplusi wrt Tc
        Matrix3 H_L1RLmplusi_Rc = Matrix3::Identity() - L1_R_Lmplusi.matrix().transpose();
        Matrix3 H_L1RLmplusi_pc = Matrix3::Zero();
        Matrix3 H_L1pLmplusi_Rc = skewSymmetric(L1_p_Lmplusi);
        Matrix3 H_L1pLmplusi_pc = L1_R_Lmplusi.matrix() - Matrix3::Identity();

        Matrix6 H_L1TLmplusi_Tc;
        H_L1TLmplusi_Tc.block(0, 0, 3, 3) = H_L1RLmplusi_Rc;
        H_L1TLmplusi_Tc.block(0, 3, 3, 3) = H_L1RLmplusi_pc;
        H_L1TLmplusi_Tc.block(3, 0, 3, 3) = H_L1pLmplusi_Rc;
        H_L1TLmplusi_Tc.block(3, 3, 3, 3) = H_L1pLmplusi_pc;

        if (H1)
            (*H1) = (Matrix(1, 6) << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GTI1).finished();
        if (H2)
            (*H2) = (Matrix(1, 6) << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GTIm).finished();
        if (H3)
            (*H3) = (Matrix(1, 3) << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_GvIm).finished();
        if (H4)
            (*H4) = (Matrix(1, 6) << H_res_xL1*H_xL1_L1TLmplusi*H_L1TLmplusi_Tc).finished();

        return res;
    }

    Vector PointToPlaneFactor::evaluateError(const Pose3 &wT1, const Pose3 &wTm, const Vector3 &wVm, const Pose3 &Tc,
                                             boost::optional<Matrix &> H1, boost::optional<Matrix &> H2,
                                             boost::optional<Matrix &> H3, boost::optional<Matrix &> H4) const {
        Vector1 error = computeErrorAndJacobians1(wT1, wTm, wVm, Tc, H1, H2, H3, H4);
        return error;
    }
}