//
// Created by cc on 7/16/20.
//

#ifndef PPPLIB_INSFUNC_H
#define PPPLIB_INSFUNC_H

#include "CmnFunc.h"
#include "GnssFunc.h"

namespace PPPLib {

    static const double kConingCoeff[5][5]={
            {2/3.0,                 0,           0,           0,           0},
            {9/20.0,          27/20.0,           0,           0,           0},
            {54/105.0,       92/105.0,   214/105.0,           0,           0},
            {250/504.0,     525/504.0,   650/504.0,  1375/504.0,           0},
            {2315/4620.0, 4558/4620.0, 7296/4620.0, 7834/4620.0,15797/4620.0}
    };

    extern void Euler2Quat(const Vector3d& euler,Quaterniond& quat);
    extern void Quat2Euler(const Quaterniond& quat,Vector3d& euler);
    extern void Dcm2Euler(const Matrix3d& dcm,Vector3d& euler);
    extern void Euler2Dcm(const Vector3d& euler,Matrix3d& dcm);
    extern void Dcm2Quat(const Matrix3d& dcm,Quaterniond& quat);
    extern void Quat2Dcm(const Quaterniond& quat,Matrix3d& dcm);
    extern void Rpy2Dcm(const Vector3d& rpy, Matrix3d& dcm);
    extern void Dcm2Rpy(const Matrix3d& dcm,Vector3d& rpy);
    extern void Rpy2Quat(const Vector3d& rpy,Quaterniond& quat);
    extern void Quat2Rpy(const Quaterniond& quat,Vector3d& rpy);
    extern void Rv2Quat(const Vector3d& dtheta,Quaterniond& quat);
    extern void Rv2Dcm(const Vector3d& dtheta, Matrix3d& dcm);
    extern void Quat2Rv(const Quaterniond& quat,Vector3d& dtheta);

    Eigen::Matrix3d VectorSkew(const Eigen::Vector3d& vec);
    Eigen::Matrix3d Quaternion2RotationMatrix(const Eigen::Quaterniond& q);
    Eigen::Quaterniond RotationMatrix2Quaternion(const Eigen::Matrix3d& m);
    Eigen::Quaterniond Euler2Quaternion(const Vector3d& rpy);
    Eigen::Matrix3d Euler2RotationMatrix(const Vector3d& rpy);
    Eigen::Vector3d RotationMatrix2Euler(const Matrix3d &m);
    Eigen::Vector3d Quaternion2Euler(const Quaterniond& q);
    Eigen::Quaterniond RotationVector2Quaternion(const Vector3d& rv);
    Eigen::Vector3d RotateVec(const Vector3d& rv,const Vector3d& v);
    Eigen::Vector3d CalculateGravity(const Vector3d coord_blh,bool is_ecef);


    typedef struct {
        Vector3d rn,vn;
        double sin_lat,cos_lat,tan_lat;
        double sin_lat2,sin_lat4;
        double sq;
        double RN,cos_lat_RNh;
        double RNh,RMh;
        double g;
        Vector3d w_n_ie,w_n_en,w_n_in,w_n_ie_n;
        Vector3d w_e_ie;
        Vector3d gn,ge,gcc_n,gcc_e;
    }tEarthPar;

    typedef struct{
        cTime t_tag;
        double dt;     // time difference related to increment distance
        double dr;     // increment of distance (m)
        double vr[3];  // wheel velocity in vehicle rear frame
    }tOdoDataUnit;

    typedef struct {
        cTime t_tag;
        Vector3d gyro;      // rate rad/s
        Vector3d acce;      // rate m/s^2

        unsigned int pps;
        unsigned int imu_cnt;

        short int odo_cnt;
        tOdoDataUnit odo;
    }tImuDataUnit;

    class cImuData{
    public:
        cImuData();
        cImuData(cTime* ts,cTime* te);
        ~cImuData();

    public:
        void SetImu(tInsConf C);
        void SetImuType(IMU_TYPE type);
        void SetImuCoordType(IMU_COORD_TYPE type);
        void SetTimeSpan(cTime* ts, cTime* te);

    public:
        cTime ts_,te_;
        IMU_TYPE imu_type_;
        IMU_COORD_TYPE imu_coord_type_;
        IMU_DATA_FORMAT data_format_;
        GYRO_DATA_FORMAT gyro_format_;
        double hz_;
        vector<tImuDataUnit> data_;
    };

    typedef struct{
        cTime t_tag;
        double dt;
        cTime pt;

        ///　in body-frame
        Vector3d fb;
        Vector3d ba,bg;
        Vector3d raw_gyro_incr,raw_acce_incr;        /// increment in body-frame
        Vector3d raw_gyro_rate,raw_acce_rate;        /// raw rate in body-frame
        Vector3d cor_gyro_rate,cor_acce_rate;        /// corrected rate in body-frame

        /// in ned-frame
        Vector3d rn,vn,an,fn,rpy_n;
        Matrix3d Cbn;
        Quaterniond Qbn;

        /// in ecef-frame
        Vector3d re,ve,ae,fe,rpy_e;
        Matrix3d Cbe;
        Quaterniond Qbe;
    }tImuInfoUnit;

    typedef struct {
        Vector3d acce_bias;      /// accelerometer biases(micro-g, converted to m/s^2, body axes)
        Vector3d gyro_bias;      /// gyro biases(deg/h, converted to rad/s. body axes)
        Matrix3d Ma;             /// Accelerometer scale factor and cross coupling errors(ppm, converted to unitless, body axes)
        Matrix3d Mg;             /// gyro scale factor and cross coupling errors (ppm, converted to unitless, body axes)
        Vector3d vel_rw;         /// accelerometer noise root PSD (mg/sqrt(Hz), converted to m/s^{-1.5})
        Vector3d ang_rw;         /// gyro noise root PSD (deg/sqrt(h),converted to rad/s^{-0.5})
        Matrix3d G_err;          /// gyro g-dependent biases (deg/h/g, converted to rad-sec/m, body axes)
        double acce_quant_level; /// acce quantization level (m/s^2)
        double gyro_quant_level; /// gyro quantization level (rad/s)
        Vector3d acce_quant_res;
        Vector3d gyro_quant_res;
    }tImuErrModel;

    class cInsMech{
    public:
        cInsMech();
        ~cInsMech();

    public:
        bool InsMechanization(bool err_model,tImuInfoUnit& pre_imu_info,tImuInfoUnit& cur_imu_info,int idx);
        Eigen::MatrixXd StateTransferMat(tPPPLibConf C,tImuInfoUnit& pre_imu_info,tImuInfoUnit& cur_imu_info,int nx,double dt);

    private:
        void RotScullCorr(tImuInfoUnit& pre_imu_info,tImuInfoUnit& cur_imu_info,double dt,double *da,double *dv);
        Eigen::Quaterniond AttitudeUpdate(tImuInfoUnit& pre_imu_info,tImuInfoUnit& cur_imu_info,double dt,Vector3d da);
        Eigen::Vector3d VelocityUpdate(tImuInfoUnit& pre_imu_info,tImuInfoUnit& cur_imu_info,double dt,Vector3d dv);
        Eigen::Vector3d PositionUpdate(const tImuInfoUnit& pre_imu_info,const Vector3d& cur_vel,double dt);
        void TraceInsMechInfo(tImuInfoUnit &imu_info,bool prior,int idx);
    };

    class cInsAlign{
    public:
        cInsAlign();
        cInsAlign(cImuData imu_data,tPPPLibConf C);
        ~cInsAlign();

    public:
        bool CoarseAlign(tImuInfoUnit& imu_info);
        bool GnssSolAlign();
        bool GnssObsAlign();

    private:
        cImuData *imu_data_;
        tPPPLibConf C_;

    public:
        int imu_idx=0;
    };

    void AdjustImuData(tImuDataUnit& imu_data,IMU_COORD_TYPE coord_type,IMU_DATA_FORMAT data_format,GYRO_DATA_FORMAT gyro_val_format,double dt);

    class cInsSim {
    public:
        cInsSim();
        cInsSim(double *ep,tPPPLibConf C);
        ~cInsSim();

    private:
        void SimImuEcef(double dt,Matrix3d Cbe,Matrix3d pre_Cbe,Vector3d vel,Vector3d pre_vel,Vector3d re,tImuDataUnit& imu_data);

    public:
        void LoadSimTrj(string trj_file);
        void LoadSimImu(string imu_file);
        void SimImuObs();
        void ImuErrGrade(IMU_TYPE grade);
        void SimImuErr(Vector3d& acce, Vector3d& gyro);

    public:
        tPPPLibConf C_;
        vector<tSolInfoUnit> sim_rej_;
        cImuData sim_imu_data_;
        cTime sim_time_;

        tImuErrModel imu_err_model_;
    };


    class cIns {
    public:
        cIns();
        cIns(tPPPLibConf C);
        ~cIns();

    public:
        void InsMech(tImuInfoUnit& cur_imu_info,const tImuInfoUnit pre_imu_info,int idx);
        Eigen::MatrixXd StateTransferMat_E(tPPPLibConf C,tImuInfoUnit& pre_imu_info,tImuInfoUnit& cur_imu_info,int nx,double dt);

    private:
        void ConingScullingCompensation(int type);
        void InsUpdateEcef();

        void SolSync(tImuInfoUnit& imu_info,COORDINATE_TYPE type);
        void TraceInsMechInfo(tImuInfoUnit imu_info,bool prior,int idx);
    public:
        tPPPLibConf C_;
        vector<Vector3d> gyros_,acces_;  // multi-sample gyro and acce increment
        tImuInfoUnit *cur_imu_info_;
        tImuInfoUnit pre_imu_info_;

    private:
        double ts_=0;
        int nts_=0;                     // 默认单子样+前一周期
        double dt_=0;

        Matrix3d Mpv_;
        Vector3d dphim_,dtheta_;
        Vector3d phim_;
        Vector3d dvbm_;

        Vector3d rpy_;
        Matrix3d Cbe_,Cbn_;
        Quaterniond Qbe_,Qbn_;
    };

}

#endif //PPPLIB_INSFUNC_H
