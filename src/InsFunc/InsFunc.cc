//
// Created by cc on 7/16/20.
//

#include "InsFunc.h"

extern const double Crf[9]={0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,-1.0};

namespace PPPLib {
    extern void Euler2Quat(const Vector3d& euler,Quaterniond& quat){
        double si=sin(euler[0]/2),ci=cos(euler[0]/2);
        double sj=sin(euler[1]/2),cj=cos(euler[1]/2);
        double sk=sin(euler[2]/2),ck=cos(euler[2]/2);
        quat.w()=ci*cj*ck+si*sj*sk;
        quat.x()=si*cj*ck-ci*sj*sk;
        quat.y()=ci*sj*ck+si*cj*sk;
        quat.z()=ci*cj*sk-si*sj*ck;
        quat.conjugate();
    }

    extern void Quat2Euler(const Quaterniond& quat,Vector3d& euler){
        euler[0]=atan2(2*(-quat.w()*quat.x()+quat.y()*quat.z()),1-2*quat.x()*quat.x()-2*quat.y()*quat.y());
        euler[1]=asin(2*(-quat.w()*quat.y()-quat.x()*quat.z()));
        euler[2]=atan2(2*(-quat.w()*quat.z()+quat.x()*quat.y()),1-2*quat.y()*quat.y()-2.0*quat.z()*quat.z());
    }

    extern void Dcm2Euler(const Matrix3d& dcm,Vector3d& euler){
        euler[0]=atan2(dcm(1,2),dcm(2,2));
        euler[1]=-asin(dcm(0,2));
        euler[2]=atan2(dcm(0,1),dcm(0,0));
    }

    //321
    extern void Euler2Dcm(const Vector3d& euler,Matrix3d& dcm){
        double sin_phi=sin(euler[0]),cos_phi=cos(euler[0]);
        double sin_theta=sin(euler[1]),cos_theta=cos(euler[1]);
        double sin_psi=sin(euler[2]),cos_psi=cos(euler[2]);

        dcm(0,0)=cos_theta*cos_psi;
        dcm(0,1)=cos_theta*sin_psi;
        dcm(0,2)=-sin_theta;
        dcm(1,0)=-cos_phi*sin_psi+sin_phi*sin_theta*cos_psi;
        dcm(1,1)=cos_phi*cos_psi+sin_phi*sin_theta*sin_psi;
        dcm(1,2)=sin_phi*cos_theta;
        dcm(2,0)=sin_phi*sin_psi+cos_phi*sin_theta*cos_psi;
        dcm(2,1)=-sin_phi*cos_psi+cos_phi*sin_theta*sin_psi;
        dcm(2,2)=cos_phi*cos_theta;
    }

    extern void Dcm2Quat(const Matrix3d& dcm,Quaterniond& quat){
        double q=0.0;
        if(dcm(0,0)>=dcm(1,1)+dcm(2,2)){
            quat.x()=0.5*sqrt(1+dcm(0,0)-dcm(1,1)-dcm(2,2));
            q=4*quat.x();
            quat.w()=(dcm(2,1)-dcm(1,2))/q;
            quat.y()=(dcm(0,1)+dcm(1,0))/q;
            quat.z()=(dcm(0,2)+dcm(2,0))/q;
        }
        else if(dcm(1,1)>=dcm(0,0)+dcm(2,2)){
            quat.y()=0.5*sqrt(1.0-dcm(0,0)+dcm(1,1)-dcm(2,2));
            q=4*quat.y();
            quat.w()=(dcm(0,2)-dcm(2,0))/q;
            quat.x()=(dcm(0,1)+dcm(1,0))/q;
            quat.z()=(dcm(1,2)+dcm(2,1))/q;
        }
        else if(dcm(2,2)>=dcm(0,0)+dcm(1,1)){
            quat.z()=0.5*sqrt(1.0-dcm(0,0)-dcm(1,1)+dcm(2,2));
            q=4.0*quat.z();
            quat.w()=(dcm(1,0)-dcm(0,1))/q;
            quat.x()=(dcm(0,2)+dcm(2,0))/q;
            quat.y()=(dcm(1,2)+dcm(2,1))/q;
        }
        else{
            quat.w()=0.5*sqrt(1.0+dcm(0,0)+dcm(1,1)+dcm(2,2));
            q=4.0*quat.w();
            quat.x()=(dcm(2,1)-dcm(1,2))/q;
            quat.y()=(dcm(0,2)-dcm(2,0))/q;
            quat.z()=(dcm(1,0)-dcm(0,1))/q;
        }
        quat.normalize();
    }

    extern void Quat2Dcm(const Quaterniond& quat,Matrix3d& dcm){
        double q11=quat.w()*quat.w(),q12=quat.w()*quat.x(),
               q13=quat.w()*quat.y(),q14=quat.w()*quat.z(),
               q22=quat.x()*quat.x(),q23=quat.x()*quat.y(),
               q24=quat.x()*quat.z(),q33=quat.y()*quat.y(),
               q34=quat.y()*quat.z(),q44=quat.z()*quat.z();
        dcm(0,0)=q11+q22-q33-q44;
        dcm(0,1)=2.0*(q23-q14);
        dcm(0,2)=2.0*(q24+q13);
        dcm(1,0)=2.0*(q23+q14);
        dcm(1,1)=q11-q22+q33-q44;
        dcm(1,2)=2.0*(q34-q12);
        dcm(2,0)=2.0*(q24-q13);
        dcm(2,1)=2.0*(q34+q12);
        dcm(2,2)=q11-q22-q33+q44;
    }

    /* convert euler attitude(roll,pitch,yaw) to DCM
     * rpy: Enb(rad)
     * dcm: Cbn
     * */
    extern void Rpy2Dcm(const Vector3d& rpy, Matrix3d& dcm){
        Matrix3d Cnb;
        Euler2Dcm(rpy,Cnb);
        dcm=Cnb.transpose();
    }

    /* convert dcm to euler attitude
     * dcm: Cbn
     * rpy: Enb
     * */
    extern void Dcm2Rpy(const Matrix3d& dcm,Vector3d& rpy){
        Dcm2Euler(dcm.transpose(),rpy);
    }

    /* convert euler attitude to quaternion
     * rpy: Enb
     * quat: Qbn
     * */
    extern void Rpy2Quat(const Vector3d& rpy,Quaterniond& quat){
        Euler2Quat(rpy,quat);
        quat.conjugate();
    }

    /* convert quaternion to euler attitude
     * quat: Qbn
     * rpy:  Enb
     * */
    extern void Quat2Rpy(const Quaterniond& quat,Vector3d& rpy){
        Quat2Euler(quat.conjugate(),rpy);
    }

    /* rotation vector (angular increment) to quaternion attitude
     * dtheta: rotation vector(angular increment)
     * quat: quaternion attitude transformation
     *
     * */
    extern void Rv2Quat(const Vector3d& dtheta,Quaterniond& quat){
        const double F1=2*1;
        const double F2=F1*2*2;
        const double F3=F2*2*3;
        const double F4=F3*2*4;
        const double F5=F4*2*5;
        double n2=SQR(dtheta.norm());
        double f;
        if(n2<(PI/180.0*PI/180.0)){
            double n4=n2*n2;
            quat.w()=1.0-n2*(1.0/F2)+n4*(1.0/F4);
            f=0.5-n2*(1.0/F3)+n4*(1.0/F5);
        }
        else{
            double n_2=sqrt(n2)/2.0;
            quat.w()=cos(n_2);
            f=sin(n_2)/n_2*0.5;
        }
        quat.x()=f*dtheta[0];
        quat.y()=f*dtheta[1];
        quat.z()=f*dtheta[2];
    }

    /*
     *
     * */
    extern void Rv2Dcm(const Vector3d& dtheta, Matrix3d& dcm){
        double norm_dtheta=dtheta.norm();
        Matrix3d Ax;
        Ax=VectorSkew(dtheta);
        if(norm_dtheta>1E-8){
            Matrix3d m1=sin(norm_dtheta)/norm_dtheta*Ax;
            Matrix3d m2=(1-cos(norm_dtheta))/SQR(norm_dtheta)*Ax*Ax;
            dcm=Matrix3d::Identity()+m1+m2;
        }
        else{
            dcm=Matrix3d::Identity()+Ax;
        }
    }

    extern void Quat2Rv(const Quaterniond& quat,Vector3d& dtheta){
        double n2=acos(fabs(quat.w()));
        double k;
        if(n2>1E-40) k=2*n2/sin(n2);
        else k=2.0;
        if(quat.w()<0.0) k=-k;
        dtheta[0]=k*quat.x();
        dtheta[1]=k*quat.y();
        dtheta[2]=k*quat.z();
    }

    // =======================================================
    Eigen::Matrix3d VectorSkew(const Eigen::Vector3d& vec){
        Eigen::Matrix3d dcm=Matrix3d::Zero();

        dcm(0,1)=-vec(2);
        dcm(0,2)=vec(1);

        dcm(1,0)=vec(2);
        dcm(1,2)=-vec(0);

        dcm(2,0)=-vec(1);
        dcm(2,1)=vec(0);
        return dcm;
    }

    Eigen::Matrix3d Quaternion2RotationMatrix(const Eigen::Quaterniond& q){
        return q.toRotationMatrix();
    }

    Eigen::Quaterniond RotationMatrix2Quaternion(const Eigen::Matrix3d& m){
        Eigen::Quaterniond q(m);
        return q;
    }

    Eigen::Quaterniond Euler2Quaternion(const Vector3d& rpy){
        Eigen::AngleAxisd roll_angle(Eigen::AngleAxisd(rpy(2),Eigen::Vector3d::UnitZ()));
        Eigen::AngleAxisd pitch_angle(Eigen::AngleAxisd(rpy(1),Eigen::Vector3d::UnitY()));
        Eigen::AngleAxisd yaw_angle(Eigen::AngleAxisd(rpy(0),Eigen::Vector3d::UnitX()));

        return roll_angle*pitch_angle*yaw_angle;
    }

    Eigen::Matrix3d Euler2RotationMatrix(const Vector3d& rpy){
        Eigen::Quaterniond aa=Euler2Quaternion(rpy);
        return Quaternion2RotationMatrix(Euler2Quaternion(rpy));
    }

    template <typename Derived>
    static Eigen::Matrix<Derived,3,1> EulerAngles(const Eigen::Matrix<Derived,3,3>& m){
        Eigen::Matrix<Derived, 3, 1> res;

        const size_t i = 2;
        const size_t j = 1;
        const size_t k = 2;
        typedef Eigen::Matrix<Derived, 2, 1> Vector2;

        res[2]=atan2(m(j,k),m(k,k));
        Derived c2=Vector2(m(i,i),m(i,j)).norm();
//        if(res[2]<Derived(0))
//        {
//            res[2]+=Derived(EIGEN_PI)*2;
//        }
        res[1]=atan2(-m(i,k),c2);
        Derived s1=sin(res[2]);
        Derived c1=cos(res[2]);
        res[0]=atan2(s1*m(k,i)-c1*m(j,i),c1*m(j,j)-s1*m(k,j));

        return res;
    }

    Eigen::Vector3d RotationMatrix2Euler(const Matrix3d &m){
//        return m.eulerAngles(2,1,0);
//        return EulerAngles(m);
        Vector3d rpy;
        rpy[0]=atan2(-m(2,0),m(2,2));
        rpy[1]=asin(m(2,1));
        rpy[2]=atan2(-m(0,1),m(1,1));
        return rpy;
    }

    Eigen::Vector3d Quaternion2Euler(const Quaterniond& q){
#if 1
        return RotationMatrix2Euler(q.toRotationMatrix());
#else
        double sinr_cosp=+2.0*(q.w()*q.x()+q.y()*q.z());
        double cosr_cosp=+1.0-2.0*(q.x()*q.x()+q.y()*q.y());
        double roll=atan2(sinr_cosp,cosr_cosp);
        double pitch,yaw;

        double sinp=+2.0*(q.w()*q.y()-q.z()*q.x());
        if(fabs(sinp)>=1.0){
            pitch=copysign(M_PI/2.0,sinp);
        }
        else{
            pitch=asin(sinp);
        }

        double siny_cosp=+2.0*(q.w()*q.z()+q.x()*q.y());
        double cosy_cosp=+1.0-2.0*(q.y()*q.y()+q.z()*q.z());
        yaw=atan2(siny_cosp,cosy_cosp);

        Vector3d rpy(roll,pitch,yaw);
        return rpy;
#endif

    }

    Eigen::Quaterniond RotationVector2Quaternion(const Vector3d& rv){
        Eigen::Quaterniond qfromrv;
        Vector3d rv_2=rv*0.5;
        double norm=rv_2.norm();
        qfromrv.w()=cos(norm);
        qfromrv.vec()=norm<1E-8?rv_2:(sin(norm)/norm)*rv_2;
        return qfromrv;
    }

    Eigen::Vector3d RotateVec(const Vector3d& rv,const Vector3d& v){
        double n2=SQR(rv.norm());
        double q1,s,n,n_2;
        if(n2<1.0E-08){
            q1=1-n2*(1.0/8-n2/384);s=1.0/2-n2*(1.0/48-n2/3840);
        }
        else{
            n=sqrt(n2);n_2=n/2;
            q1=cos(n_2);s=sin(n_2)/n;
        }
        double q2=s*rv[0],q3=s*rv[1],q4=s*rv[2];
        double qo1,qo2,qo3,qo4;
        qo1=-q2*v[0]-q3*v[1]-q4*v[2];
        qo2= q1*v[0]+q3*v[2]-q4*v[1];
        qo3= q1*v[1]+q4*v[0]-q2*v[2];
        qo4= q1*v[2]+q2*v[1]-q3*v[0];

        Vector3d vo(0,0,0);
        vo[0]=-qo1*q2+qo2*q1-qo3*q4+qo4*q3;
        vo[1]=-qo1*q3+qo3*q1-qo4*q2+qo2*q4;
        vo[2]=-qo1*q4+qo4*q1-qo2*q3+qo3*q2;

        return vo;
    }

    cImuData::cImuData(){}

    cImuData::cImuData(PPPLib::cTime *ts, PPPLib::cTime *te){
        if(ts) ts_=*ts;
        if(te) te_=*te;
    }

    cImuData::~cImuData() {data_.clear();}

    void cImuData::SetImu(tInsConf C){
        imu_type_=C.imu_type;
        imu_coord_type_=C.coord_type;
        data_format_=C.data_format;
        gyro_format_=C.gyro_val_format;
        hz_=C.sample_rate;
    }

    void cImuData::SetImuType(PPPLib::IMU_TYPE type) {imu_type_=type;}

    void cImuData::SetImuCoordType(PPPLib::IMU_COORD_TYPE type) {imu_coord_type_=type;}

    void cImuData::SetTimeSpan(cTime *ts, cTime *te) {
        if(ts) ts_=*ts;
        if(te) te_=*te;
    }

    cInsMech::cInsMech() {}

    cInsMech::~cInsMech() {}

    void cInsMech::RotScullCorr(PPPLib::tImuInfoUnit &pre_imu_info, PPPLib::tImuInfoUnit &cur_imu_info, double dt,double *da,double *dv) {
        Vector3d pre_da,pre_dv,cur_da,cur_dv;
        int i;

        pre_da=pre_imu_info.cor_gyro_rate*dt;
        pre_dv=pre_imu_info.cor_acce_rate*dt;

        cur_da=cur_imu_info.cor_gyro_rate*dt;
        cur_dv=cur_imu_info.cor_acce_rate*dt;

        Vector3d t1,t2,t3,t4;
        CrossVec3(cur_da.data(),cur_dv.data(),t1.data());
        CrossVec3(pre_da.data(),cur_dv.data(),t2.data());
        CrossVec3(pre_dv.data(),cur_da.data(),t3.data());
        CrossVec3(cur_da.data(),t1.data(),t4.data());

        double a1,a2;
        double b=cur_da.norm();
        if(fabs(b)>1E-6){
            a1=(1.0-cos(b))/SQR(b);
            a2=1.0/SQR(b)*(1.0-sin(b)/b);
        }
        else{
            a1=0.5-SQR(b)/24.0+SQR(SQR(b))/720.0;
            a2=1.0/6.0-SQR(b)/120.0+SQR(SQR(b))/5040.0;
        }

        for(i=0;i<3&&dv;i++){
            dv[i]=a1*t1[i]+a2*t4[i]+1.0/12.0*(t2[i]+t3[i]);
        }

        if(da){
            CrossVec3(pre_da,cur_da,da);
            for(i=0;i<3;i++) da[i]*=1.0/12.0;
        }
    }

    Eigen::Quaterniond cInsMech::AttitudeUpdate(PPPLib::tImuInfoUnit &pre_imu_info,
                                                PPPLib::tImuInfoUnit &cur_imu_info,double dt,Vector3d da) {
        Vector3d theta_k(cur_imu_info.cor_gyro_rate*dt),theta_k_1(pre_imu_info.cor_gyro_rate*dt);

        //等效旋转矢量
        Vector3d cur_phi=theta_k+da; //单子样+前一周期
        Quaterniond quat_bb=RotationVector2Quaternion(cur_phi);

        Vector3d wiee(0,0,-OMGE_GPS);
        Vector3d zeta=wiee*dt;
        Quaterniond quat_ee=RotationVector2Quaternion(zeta);
        Quaterniond quat_k_1=RotationMatrix2Quaternion(pre_imu_info.Cbe.transpose()).conjugate();
        Quaterniond qbn_k=quat_ee*quat_k_1*quat_bb;

        return qbn_k.normalized();
    }

    Eigen::Vector3d cInsMech::VelocityUpdate(PPPLib::tImuInfoUnit &pre_imu_info,
                                             PPPLib::tImuInfoUnit &cur_imu_info,double dt,Vector3d dv) {
        Vector3d pos=pre_imu_info.re, vel=pre_imu_info.ve;
        Vector3d wiee(0,0,OMGE_GPS);
        Vector3d theta_k(cur_imu_info.cor_gyro_rate*dt),theta_k_1(pre_imu_info.cor_gyro_rate*dt);
        Vector3d vb_k(cur_imu_info.cor_acce_rate*dt+dv),vb_k_1(pre_imu_info.cor_acce_rate*dt);

        Vector3d coord_blh=Xyz2Blh(pos);
        Matrix3d Cen=CalcCen(coord_blh,COORD_NED);
        Vector3d ge=Cen.transpose()*CalculateGravity(coord_blh,false);
        Vector3d omgea_n=wiee*2.0;
        Vector3d delta_gcor=(ge-omgea_n.cross(vel))*dt;

        Matrix3d Cee=Matrix3d::Identity()-VectorSkew(wiee*0.5*dt);

//        Vector3d vrot=theta_k.cross(vb_k)*0.5;
//        Vector3d vscul=(theta_k_1.cross(vb_k)+vb_k_1.cross(theta_k))/12.0;

        Quaterniond pre_quat=RotationMatrix2Quaternion(pre_imu_info.Cbe);
        Matrix3d Cbe=pre_quat.toRotationMatrix();
        Vector3d delta_ve=Cee*Cbe*(vb_k);

        return (vel+delta_gcor+delta_ve);
    }

    Eigen::Vector3d cInsMech::PositionUpdate(const tImuInfoUnit& pre_imu_info,const Vector3d& cur_vel, double dt) {
        Eigen::Vector3d pos;
        pos=(pre_imu_info.ve+cur_vel)*0.5*dt+pre_imu_info.re;
        return pos;
    }

    bool cInsMech::InsMechanization(bool err_model,tImuInfoUnit &pre_imu_info, tImuInfoUnit &cur_imu_info,int idx) {

        for(int i=0;i<3;i++){
            if(isnan(pre_imu_info.re[i])||isnan(pre_imu_info.ve[i])||
               isinf(pre_imu_info.re[i])||isinf(pre_imu_info.ve[i])){
                LOG(ERROR)<<cur_imu_info.t_tag.GetTimeStr(4)<<" "<<idx<<" NUMERIC ERROR";
                return false;
            }
        }
        TraceInsMechInfo(pre_imu_info,true,idx);

        double dt=cur_imu_info.t_tag.TimeDiff(pre_imu_info.t_tag.t_);
        if(dt>60.0||fabs(dt)<1E-6){
            cur_imu_info.dt=dt;
            cur_imu_info.pt=pre_imu_info.pt;
            LOG(WARNING)<<"TIME DIFFERENCE TOO LARGER";
            return false;
        }

        if(err_model){

        }else{
            cur_imu_info.cor_gyro_rate=cur_imu_info.raw_gyro_incr-pre_imu_info.bg;
            cur_imu_info.cor_acce_rate=cur_imu_info.raw_acce_incr-pre_imu_info.ba;
        }

        Vector3d da,dv;
        RotScullCorr(pre_imu_info,cur_imu_info,dt,da.data(),dv.data());
        cur_imu_info.Cbe=Quaternion2RotationMatrix(AttitudeUpdate(pre_imu_info,cur_imu_info,dt,da));
        cur_imu_info.ve=VelocityUpdate(pre_imu_info,cur_imu_info,dt,dv);
        cur_imu_info.re=PositionUpdate(pre_imu_info,cur_imu_info.ve,cur_imu_info.t_tag.TimeDiff(pre_imu_info.t_tag.t_));

        Vector3d blh=Xyz2Blh(cur_imu_info.re);
        Matrix3d Cen=CalcCen(blh,COORD_NED);
        Matrix3d Cnb=cur_imu_info.Cbe.transpose()*Cen.transpose();
        cur_imu_info.rpy_n=RotationMatrix2Euler(Cnb);
        cur_imu_info.vn=Cen*cur_imu_info.ve;
        cur_imu_info.rn=Cen*cur_imu_info.re;

        TraceInsMechInfo(cur_imu_info,false,idx);
        cur_imu_info.dt=dt;
        cur_imu_info.pt=pre_imu_info.t_tag;

        return true;
    }

    Eigen::MatrixXd cInsMech:: StateTransferMat(tPPPLibConf C,PPPLib::tImuInfoUnit &pre_imu_info,
                                               PPPLib::tImuInfoUnit &cur_imu_info,int nx,double dt) {
        using Eigen::Matrix3d;
        using Eigen::MatrixXd;
        using Eigen::Vector3d;

        auto &vel=cur_imu_info.ve;
        auto &fb=cur_imu_info.cor_acce_rate;
        auto &wb=cur_imu_info.cor_gyro_rate;
        Vector3d wiee(0,0,OMGE_GPS);
        auto &Cbe=cur_imu_info.Cbe;

        MatrixXd F=MatrixXd::Zero(nx,nx);

        int ip=0;
        int iv=3;
        int ia=6;
        int iba=9;
        int ibg=12;
        int isa=0,isg=0,ira=0,irg=0,ilev=0;
        if(C.insC.est_sa) isa=ibg+3;
        else isa=ibg;
        if(C.insC.est_sg) isg=isa+3;
        else isg=isa;
        if(C.insC.est_ra) ira=isg+3;
        else ira=isg;
        if(C.insC.est_rg) irg=ira+6;
        else irg=ira;
        if(C.insC.est_level) ilev=irg+6;
        else ilev=irg;

        //position-velocity  Mpv
        F.block<3,3>(ip,iv)=Matrix3d::Identity();

        //velocity-velocity  Mvv
        F.block<3,3>(iv,iv)=(-2.0*VectorSkew(wiee));
        //velocity-attitude  Mva
        F.block<3,3>(iv,ia)=-VectorSkew(Cbe*fb);
        //velocity-ba
        F.block<3,3>(iv,iba)=-Cbe;

        //attitude-attitude  Maa
        F.block<3,3>(ia,ia)=-1.0*VectorSkew(wiee);
        //attitute-bg
        F.block<3,3>(ia,ibg)=Cbe;
//        cout<<F<<endl;

        //ba-ba
//        F.block<3,3>(iba,iba)=Matrix3d::Identity()*(-fabs(dt)/C.insC.correction_time_ba);
        F.block<3,3>(iba,iba)=Matrix3d::Zero();
        //bg-bg
//        F.block<3,3>(ibg,ibg)=Matrix3d::Identity()*(-fabs(dt)/C.insC.correction_time_bg);
        F.block<3,3>(ibg,ibg)=Matrix3d::Zero();

//        cout<<MatrixXd::Identity(nx,nx)+F*dt<<endl;

        return MatrixXd::Identity(nx,nx)+F*dt;
    }

    Eigen::Vector3d CalculateGravity(const Vector3d coord_blh,bool is_ecef){
        Vector3d blh=coord_blh;
        Vector3d re=Blh2Xyz(blh);
        Vector3d g(0,0,0);
        if(is_ecef){
            double mag_r=re.norm();
            if(mag_r==0,0) return g;
            else{

                /// 2.142
                double J2=1.082627E-3;
                double z_scale=5.0*SQR(re[2]/mag_r);
                Vector3d gamma(0,0,0);
                double fact=-MU_GPS/(mag_r*mag_r*mag_r);
                double fact1=1.5*J2*SQR(WGS84_EARTH_LONG_RADIUS/mag_r);
                gamma[0]=fact*(re[0]+fact1*(1.0-z_scale)*re[0]);
                gamma[1]=fact*(re[1]+fact1*(1.0-z_scale)*re[1]);
                gamma[2]=fact*(re[2]+fact1*(3.0-z_scale)*re[2]);

                /// add centripetal acceleration 2.133
                g[0]=gamma[0]+SQR(OMGE_GPS)*re[0];
                g[1]=gamma[1]+SQR(OMGE_GPS)*re[1];
                g[2]=gamma[2];

                return g;
            }

        }
        else{
            double gn = 9.7803267715 * (1 + 0.0052790414 * sin(coord_blh(0)) * sin(coord_blh(0)) + 0.0000232719 * pow(sin(coord_blh(0)),3) * sin(coord_blh(0)));
            gn += (-0.0000030876910891 + 0.0000000043977311 * sin(coord_blh(0)) * sin(coord_blh(0))) * coord_blh(2);
            gn += 0.0000000000007211 * coord_blh(2) * coord_blh(2);
            Vector3d gn_vec{0, 0, gn};
            return gn_vec;
        }
    }

    void cInsMech::TraceInsMechInfo(PPPLib::tImuInfoUnit &imu_info,bool prior,int idx) {
        LOG(DEBUG)<<"INS MECHANIZATION"<<(prior?"- ":"+ ")<< "("<<idx<<"): "<<imu_info.t_tag.GetTimeStr(4);
        LOG(DEBUG)<<"   "<<"GYRO VALUE: "<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.cor_gyro_rate.transpose()<<" rad/s";
        LOG(DEBUG)<<"   "<<"ACCE VALUE: "<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.cor_acce_rate.transpose()<<" m/s^2";
        LOG(DEBUG)<<"   "<<"ATTITUDE:(n)"<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.rpy_n.transpose()*R2D<<" deg";
        LOG(DEBUG)<<"   "<<"VELOCITY:(n)"<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.vn.transpose()<<" m/s";
        LOG(DEBUG)<<"   "<<"POSITION:(e)"<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.re.transpose()<<" m";
        LOG(DEBUG)<<"   "<<"GYRO BIAS:  "<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.bg.transpose()<<" rad/s";
        LOG(DEBUG)<<"   "<<"ACCE BIAS:  "<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.ba.transpose()<<" m/s^2";
    }

    cInsAlign::cInsAlign() {}

    cInsAlign::cInsAlign(cImuData imu_data,tPPPLibConf C) {
        imu_data_=&imu_data;
        C_=C;
    }

    cInsAlign::~cInsAlign() {}

    bool cInsAlign::CoarseAlign(tImuInfoUnit &imu_info) {

    }

    /// 程序在ENU-RFU框架下进行ins编排,
    /// imu数据应存储为RFU坐标系并表示为增量形式
    void AdjustImuData(tImuDataUnit& imu_data,IMU_COORD_TYPE coord_type,IMU_DATA_FORMAT data_format,GYRO_DATA_FORMAT gyro_val_format,double dt) {
        Vector3d gyro,acce;

#if 1
        if(coord_type==IMU_COORD_RFU){
            gyro=imu_data.gyro;
            acce=imu_data.acce;
            MatMul("NN",3,1,3,1.0,Crf,gyro.data(),0.0,imu_data.gyro.data());
            MatMul("NN",3,1,3,1.0,Crf,acce.data(),0.0,imu_data.acce.data());
        }

        if(data_format==IMU_FORMAT_INCR){
            for(int j=0;j<3;j++) imu_data.acce[j]/=dt;
            for(int j=0;j<3;j++) imu_data.gyro[j]/=dt;
        }
        if(gyro_val_format==GYRO_FORMAT_DEG){
            for(int j=0;j<3;j++) imu_data.gyro[j]*=D2R;
        }
#endif
    }

    cInsSim::cInsSim() {}

    cInsSim::cInsSim(double *ep,tPPPLibConf C) {
        sim_time_.Epoch2Time(ep);
        C_=C;
    }

    cInsSim::~cInsSim() {}

    void cInsSim::LoadSimTrj(string trj_file) {
        double data[11]={0};
        tSolInfoUnit sol_info={0};
        int i,week;
        double sow;
        ifstream inf;
        string line_str,part_str;
        cTime t_tag=sim_time_;
        Vector3d blh(0,0,0),vel(0,0,0),rpy(0,0,0);
        Matrix3d Cen;

        inf.open(trj_file);
        if(!inf.is_open()){
            LOG(ERROR)<<"FILE OPEN ERROR"<<trj_file;
        }

        int line_num=0;
        while(getline(inf,line_str)&&!inf.eof()){
            line_num++;
            istringstream read_str(line_str);
            t_tag=sim_time_;

            for(auto &j:data) j=0.0;

            for(i=0;i<11;i++){
                getline(read_str,part_str,',');
                if(line_num==1&&i==0) continue;
                data[i]=atof(part_str.c_str());
            }

#if 0
            sol_info.att[0]=data[0];  // roll
            sol_info.att[1]=data[1];  // pitch
            sol_info.att[2]=data[2];  // yaw
            sol_info.vel[0]=data[3];  // m/s
            sol_info.vel[1]=data[4];  // m/s
            sol_info.vel[2]=data[5];  // m/s
            blh[0]=data[6];           // lat    rad
            blh[1]=data[7];           // lon    rad
            blh[2]=data[8];           // height m
            sol_info.pos=Blh2Xyz(blh);
            sol_info.t_tag=t_tag+data[9];

#else
            sol_info.t_tag=t_tag+data[0];
            blh[0]=data[1]*D2R;
            blh[1]=data[2]*D2R;
            blh[2]=data[3];
            sol_info.pos=Blh2Xyz(blh);
            vel[0]=data[4];         //north velocity (m/s)
            vel[1]=data[5];
            vel[2]=data[6];
            Cen=CalcCen(blh,COORD_NED);
            sol_info.vel=Cen.transpose()*vel;
            sol_info.att[0]=data[7]*D2R;  //roll
            sol_info.att[1]=data[8]*D2R;  //pitch
            sol_info.att[2]=data[9]*D2R;  //yaw of body w.r.t NED
#endif
            sol_info.stat=SOL_SIM;
            sim_rej_.push_back(sol_info);
        }

        inf.close();
    }

    void cInsSim::LoadSimImu(string imu_file) {
        double data[8]={0};
        tImuDataUnit imu_data={0};
        ifstream inf;
        string line_str,part_str;
        cTime t_tag=sim_time_;
        int i;
        inf.open(imu_file);
        if(!inf.is_open()){
            LOG(ERROR)<<"FILE OPEN ERROR"<<imu_file;
        }

        int line_num=0;

        while(getline(inf,line_str)&&!inf.eof()){
            t_tag=sim_time_;
            line_num++;
            istringstream read_str(line_str);

            for(auto &j:data) j=0.0;

            for(i=0;i<8;i++){
                getline(read_str,part_str,',');
                data[i]=atof(part_str.c_str());
            }

            imu_data.gyro[0]=data[0];
            imu_data.gyro[1]=data[1];
            imu_data.gyro[2]=data[2];
            imu_data.acce[0]=data[3];
            imu_data.acce[1]=data[4];
            imu_data.acce[2]=data[5];
            imu_data.t_tag=t_tag+data[6];
            sim_imu_data_.data_.push_back(imu_data);
        }
        inf.close();
    }

    /// calculate specific force and angular rate
    void cInsSim::SimImuEcef(double dt,Matrix3d Cbe,Matrix3d pre_Cbe,Vector3d vel,Vector3d pre_vel,Vector3d re,tImuDataUnit& imu_data) {
        double alpha_ie=OMGE_GPS*dt;
        Vector3d blh=Xyz2Blh(re);

        Matrix3d C_earth=Matrix3d::Identity();
        C_earth(0,0)=cos(alpha_ie),C_earth(0,1)=sin(alpha_ie);
        C_earth(1,0)=-sin(alpha_ie),C_earth(1,1)=cos(alpha_ie);

        Matrix3d C_old_new; /// C_b(-)_b(+)
        C_old_new=Cbe.transpose()*C_earth*pre_Cbe;

        /// calculate the approximate angular rate w.r.t an inertial frame
        Vector3d alpha_ib_b;
        alpha_ib_b[0]=0.5*(C_old_new(1,2)-C_old_new(2,1));
        alpha_ib_b[1]=0.5*(C_old_new(2,0)-C_old_new(0,2));
        alpha_ib_b[2]=0.5*(C_old_new(0,1)-C_old_new(1,0));

        /// calculate and apply the scaling factor
        double temp=acos(0.5*(C_old_new(0,0)+C_old_new(1,1)+C_old_new(2,2)-1.0));
        if(temp>2E-05){
            alpha_ib_b*=(temp/sin(temp));
        }

        Vector3d omeg_ib_b=alpha_ib_b/dt;

        /// calculate the specific force resolved about ECEF-frame axes
        Vector3d w_ie(0,0,OMGE_GPS);
        Vector3d f_ib_e=((vel-pre_vel)/dt)-CalculateGravity(blh,true)+2.0*VectorSkew(w_ie)*pre_vel;

        /// calculate the average body-to-ecef coordinate transformation matrix over the update interval
        Matrix3d avg_Cbe;
        double mag_alpha=alpha_ib_b.norm();
        Matrix3d Alpha_ib_b=VectorSkew(alpha_ib_b);
        if(mag_alpha>1.0E-08){
            Matrix3d C_avgb_oldb;
            C_avgb_oldb=Matrix3d::Identity()+(1.0-cos(mag_alpha))/SQR(mag_alpha)*Alpha_ib_b+(1.0-sin(mag_alpha)/mag_alpha)/SQR(mag_alpha)*Alpha_ib_b*Alpha_ib_b;
            avg_Cbe=pre_Cbe*C_avgb_oldb-0.5*VectorSkew(w_ie)*dt*pre_Cbe;
        }
        else{
            avg_Cbe=pre_Cbe-0.5*VectorSkew(w_ie*dt)*pre_Cbe;
        }

        Vector3d f_ib_b=avg_Cbe.inverse()*f_ib_e;

        imu_data.acce=f_ib_b;
        imu_data.gyro=omeg_ib_b;
    }

    void cInsSim::SimImuObs() {
        int i;

        Vector3d blh,vel,pre_vel,re;
        Matrix3d pre_Cbe,Cbe,Cbn,Cen;

        /// first epoch info
        blh=Xyz2Blh(sim_rej_[0].pos);
        Cen=CalcCen(blh,COORD_NED);
        Rpy2Dcm(sim_rej_[0].att,Cbn);
        pre_Cbe=Cen.transpose()*Cbn;
        pre_vel=sim_rej_[0].vel;
        tImuDataUnit imu_data={0};

        /// epoch loop
        for(i=1;i<sim_rej_.size();i++){
            imu_data.t_tag=sim_rej_[i].t_tag;
            re=sim_rej_[i].pos;
            blh=Xyz2Blh(re);
            Cen=CalcCen(blh,COORD_NED);
            Rpy2Dcm(sim_rej_[i].att,Cbn);
            Cbe=Cen.transpose()*Cbn;
            vel=sim_rej_[i].vel;
            SimImuEcef(1.0/C_.insC.sample_rate,Cbe,pre_Cbe,vel,pre_vel,re,imu_data);
            SimImuErr(imu_data.acce,imu_data.gyro);

            pre_Cbe=Cbe;
            pre_vel=vel;
            sim_imu_data_.data_.push_back(imu_data);
        }
    }

    void cInsSim::ImuErrGrade(IMU_TYPE grade) {
        Vector3d err_gyro,err_acce;

        if(grade==IMU_SIM_AIAV){
            /// aviation-grade IMU error model
            imu_err_model_.acce_bias<<30*MICRO_G,-45*MICRO_G,26*MICRO_G;
            imu_err_model_.gyro_bias<<-0.0009*ARC_DEG_PER_HOUR,0.0013*ARC_DEG_PER_HOUR,-0.0008*ARC_DEG_PER_HOUR;
            imu_err_model_.Ma<<  100*PPM, -120*PPM,  80*PPM,
                                  -60*PPM, -120*PPM,-100*PPM,
                                 -100*PPM,   40*PPM,  90*PPM;
            imu_err_model_.Mg<<  8*PPM,-120*PPM, 100*PPM,
                                  0*PPM,  -6*PPM, -60*PPM,
                                  0*PPM,   0*PPM,  -7*PPM;
            double a=ARC_DEG_PER_HOUR/NORM_GRAVITY_VALUE;
            imu_err_model_.G_err<<0.0,0.0,0.0,
                                   0.0,0.0,0.0,
                                   0.0,0.0,0.0;
            imu_err_model_.vel_rw<<20*UG_PER_SQRT_HZ,20*UG_PER_SQRT_HZ,20*UG_PER_SQRT_HZ;
            imu_err_model_.ang_rw<<0.002*ARC_DEG_PER_SQRT_HOUR,0.002*ARC_DEG_PER_SQRT_HOUR,0.002*ARC_DEG_PER_SQRT_HOUR;
            imu_err_model_.acce_quant_level=5E-5;
            imu_err_model_.gyro_quant_level=1E-6;
        }
        else if(grade==IMU_SIM_TACTICAL){
            /// Tactical-grade IMU error model
            imu_err_model_.acce_bias<<900*MICRO_G,-1300*MICRO_G,800*MICRO_G;
            imu_err_model_.gyro_bias<<-9*ARC_DEG_PER_HOUR,13*ARC_DEG_PER_HOUR,-8*ARC_DEG_PER_HOUR;
            imu_err_model_.Ma<<  500*PPM, -300*PPM, 200*PPM,
                                -150*PPM, -600*PPM, 250*PPM,
                                 -250*PPM,  100*PPM, 450*PPM;
            imu_err_model_.Mg<<  400*PPM,-300*PPM, 250*PPM,
                                    0,-300*PPM,-150*PPM,
                                    0,   0,-350*PPM;
            double a=D2R/(3600.0*NORM_GRAVITY_VALUE);
            imu_err_model_.G_err<<0.9*a,-1.1*a,-0.6*a,
                                   -0.5*a,1.9*a,-1.6*a,
                                   0.3*a,1.1*a,-1.3*a;
            imu_err_model_.vel_rw<<100*UG_PER_SQRT_HZ,100*UG_PER_SQRT_HZ,100*UG_PER_SQRT_HZ;
            imu_err_model_.ang_rw<<0.1*ARC_DEG_PER_SQRT_HOUR,0.1*ARC_DEG_PER_SQRT_HOUR,0.1*ARC_DEG_PER_SQRT_HOUR;
            imu_err_model_.acce_quant_level=1E-2;
            imu_err_model_.gyro_quant_level=2E-4;
        }
        else if(grade==IMU_SIM_MEMS){
            /// consumer-grade IMU error model
            imu_err_model_.acce_bias<<9000*MICRO_G,-13000*MICRO_G,8000*MICRO_G;
            imu_err_model_.gyro_bias<<-180*ARC_DEG_PER_HOUR,260*ARC_DEG_PER_HOUR,-160*ARC_DEG_PER_HOUR;
            imu_err_model_.Ma<<   50000*PPM, -15000*PPM, 10000*PPM,
                                  -7500*PPM, -60000*PPM,-12500*PPM,
                                 -12500*PPM,   5000*PPM, 20000*PPM;
            imu_err_model_.Mg<<  40000*PPM,-14000*PPM, 12500*PPM,
                                      0,-30000*PPM, -7500*PPM,
                                      0,     0,-17500*PPM;
            double a=ARC_DEG_PER_HOUR/NORM_GRAVITY_VALUE;
            imu_err_model_.G_err<<90*a,-110*a,-60*a,
                                  -50*a,190*a,-160*a,
                                   30*a,110*a,-130*a;
            imu_err_model_.vel_rw<<1000*UG_PER_SQRT_HZ,1000*UG_PER_SQRT_HZ,1000*UG_PER_SQRT_HZ;
            imu_err_model_.ang_rw<<1.0*ARC_DEG_PER_SQRT_HOUR,1.0*ARC_DEG_PER_SQRT_HOUR,1.0*ARC_DEG_PER_SQRT_HOUR;
            imu_err_model_.acce_quant_level=1E-1;
            imu_err_model_.gyro_quant_level=2E-3;
        }

        imu_err_model_.gyro_quant_res<<0.0,0.0,0.0;
        imu_err_model_.acce_quant_res<<0.0,0.0,0.0;
    }

    void cInsSim::SimImuErr(Vector3d &acce, Vector3d &gyro) {
        double ts=1.0/sim_imu_data_.hz_;

        Vector3d acce_noise(0,0,0),gyro_noise(0,0,0);

        /// generate noise
        for(int i=0;i<3;i++){
            acce_noise[i]=imu_err_model_.vel_rw[i]/sqrt(ts)*RandNorm(0.5);
            gyro_noise[i]=imu_err_model_.ang_rw[i]/sqrt(ts)*RandNorm(0.5);
        }

        /// calculate accelerometer and gyro outputs
        Vector3d f_ib_b(0,0,0),omega_ib_b(0,0,0);
        f_ib_b=imu_err_model_.acce_bias+(Matrix3d::Identity()+imu_err_model_.Ma)*acce;
        omega_ib_b=imu_err_model_.gyro_bias+(Matrix3d::Identity()+imu_err_model_.Mg)*gyro+imu_err_model_.G_err*acce;

        /// quantize accelerometer outputs
        if(imu_err_model_.acce_quant_level>0.0){
            for(int i=0;i<3;i++){
                acce[i]=imu_err_model_.acce_quant_level*round((f_ib_b[i]+imu_err_model_.acce_quant_res[i])/imu_err_model_.acce_quant_level);
            }
            imu_err_model_.acce_quant_res=f_ib_b+imu_err_model_.acce_quant_res-acce;
        }
        else{
            acce=f_ib_b;
        }

        /// quantize gyro outputs
        if(imu_err_model_.gyro_quant_level>0.0){
            for(int i=0;i<3;i++){
                gyro[i]=imu_err_model_.gyro_quant_level*round((omega_ib_b[i]+imu_err_model_.gyro_quant_res[i])/imu_err_model_.gyro_quant_level);
            }
            imu_err_model_.gyro_quant_res=omega_ib_b+imu_err_model_.gyro_quant_res-gyro;
        }
        else{
            gyro=omega_ib_b;
        }
    }

    cIns::cIns() {}

    cIns::cIns(PPPLib::tPPPLibConf C) {
        C_=C;
        ts_=1.0/C_.insC.sample_rate;
        nts_=C_.insC.sample_number;
    }

    cIns::~cIns() {}

    void cIns::SolSync(PPPLib::tImuInfoUnit &imu_info,COORDINATE_TYPE coord) {
        Matrix3d Cen;
        if(coord==COORD_NED){
            /// NED to ECEF
            Cen=CalcCen(cur_imu_info_->rn,COORD_NED);
            cur_imu_info_->re=Blh2Xyz(cur_imu_info_->rn);
            cur_imu_info_->ve=Cen.transpose()*cur_imu_info_->vn;
            cur_imu_info_->Cbe=Cen.transpose()*cur_imu_info_->Cbn;
            Dcm2Rpy(cur_imu_info_->Cbe,cur_imu_info_->rpy_e);
        }
        else if(coord==COORD_XYZ){
            /// ECEF to NED
            cur_imu_info_->rn=Xyz2Blh(cur_imu_info_->re);
            Cen=CalcCen(cur_imu_info_->rn,COORD_NED);
            cur_imu_info_->vn=Cen*cur_imu_info_->ve;
            cur_imu_info_->Cbn=Cen*cur_imu_info_->Cbe;
            Dcm2Rpy(cur_imu_info_->Cbn,cur_imu_info_->rpy_n);
        }
    }

    void cIns::TraceInsMechInfo(tImuInfoUnit imu_info,bool prior,int idx) {
        LOG(DEBUG)<<"INS MECHANIZATION"<<(prior?"- ":"+ ")<< "("<<idx<<"): "<<imu_info.t_tag.GetTimeStr(4);
        LOG(DEBUG)<<"   "<<"GYRO VALUE(b):"<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.cor_gyro_rate.transpose()<<" rad/s";
        LOG(DEBUG)<<"   "<<"ACCE VALUE(b):"<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.cor_acce_rate.transpose()<<" m/s^2";
        LOG(DEBUG)<<"   "<<"ATTITUDE(ned):"<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.rpy_n.transpose()*R2D<<" deg";
        LOG(DEBUG)<<"   "<<"ATTITUDE(e)  :"<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.rpy_e.transpose()*R2D<<" deg";
        LOG(DEBUG)<<"   "<<"VELOCITY(ned):"<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.vn.transpose()<<" m/s";
        Vector3d rn=imu_info.rn;
        rn[0]*=R2D;rn[1]*=R2D;
        LOG(DEBUG)<<"   "<<"POSITION(ned):"<<setw(13)<<std::fixed<<setprecision(6)<<rn.transpose()<<" deg deg m";
        LOG(DEBUG)<<"   "<<"GYRO BIAS(b) :"<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.bg.transpose()/ARC_DEG_PER_HOUR<<" deg/h";
        LOG(DEBUG)<<"   "<<"ACCE BIAS(b) :"<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.ba.transpose()/MICRO_G<<" ug";
    }

    /* ins mech in ENU */
    void cIns::InsMech(PPPLib::tImuInfoUnit &cur_imu_info,const tImuInfoUnit pre_imu_info,int idx) {
        cur_imu_info_=&cur_imu_info;
        pre_imu_info_=pre_imu_info;

        TraceInsMechInfo(pre_imu_info_,true,idx);

        dt_=ts_*nts_;
        /// corrected
        if(C_.insC.sample_number==1){
            ConingScullingCompensation(2);
        }
        else{
            ConingScullingCompensation(0);
        }

        if(C_.insC.clp_bg&&C_.insC.clp_fact!=0.0){
            phim_-=pre_imu_info_.bg*dt_;
        }
        if(C_.insC.clp_ba&&C_.insC.clp_fact!=0.0){
            dvbm_-=pre_imu_info_.ba*dt_;
        }
        cur_imu_info_->cor_gyro_rate=phim_/dt_;
        cur_imu_info_->fb=cur_imu_info_->cor_acce_rate=dvbm_/dt_;

        InsUpdateEcef();

        TraceInsMechInfo(*cur_imu_info_,false,idx);

        gyros_.clear();acces_.clear();
    }

    /* 1: for poplinomial compensation method
     * 0: for optimal coning compensation method
     * 2: single sample+previous sample
     * */
    void cIns::ConingScullingCompensation(int type) {
        int i,j;
        phim_<<0,0,0;
        dvbm_<<0,0,0;
        if(nts_==0) nts_=2;
        MatrixXd cm(1,nts_-1),sm(1,nts_-1);
        // coning compensation
        Vector3d dphim(0,0,0);
        if(C_.insC.sample_number==1&&type==2){
            // single sample+previous sample
            phim_=cur_imu_info_->raw_gyro_incr;
            CrossVec3(pre_imu_info_.raw_gyro_incr.data(),cur_imu_info_->raw_gyro_incr.data(),dphim.data());
            dphim=1.0/12.0*dphim;
        }
        else if(type==1){
            // for optimal coning compensation method
        }
        else if(type==0){
            for(i=0;i<nts_;i++){
                phim_[0]+=gyros_[i][0];
                phim_[1]+=gyros_[i][1];
                phim_[2]+=gyros_[i][2];
            }

            MatrixXd a(1,nts_-1),b(nts_-1,nts_-1);
            cm=MatrixXd::Zero(1,nts_-1);
            a=MatrixXd::Zero(1,nts_-1);b=MatrixXd::Zero(nts_-1,3);

            for(i=0;i<nts_-1;i++){
                a(0,i)=kConingCoeff[nts_-2][i];
            }
            for(i=0;i<nts_-1;i++){
                for(j=0;j<3;j++){
                    b(i,j)=gyros_[i][j];
                }
            }
            cm=a*b;
            CrossVec3(cm.data(),gyros_[nts_-1].data(),dphim.data());
            phim_+=dphim;
            dphim_=dphim;
        }
        else{
            for(i=0;i<nts_;i++){
                phim_[0]+=gyros_[i][0];
                phim_[1]+=gyros_[i][1];
                phim_[2]+=gyros_[i][2];
            }
        }

        // sculling compensation
        Vector3d scull(0,0,0);
        if(C_.insC.sample_number==1&&type==2){
            // single sample+previous sample
            dvbm_=cur_imu_info_->raw_acce_incr;
            Vector3d s1(0,0,0),s2(0,0,0);
            CrossVec3(pre_imu_info_.raw_gyro_incr.data(),cur_imu_info_->raw_acce_incr.data(),s1.data());
            CrossVec3(pre_imu_info_.raw_acce_incr.data(),cur_imu_info_->raw_gyro_incr.data(),s2.data());
            scull=1/12.0*(s1+s2);
        }
        else if(type==1){

        }
        else if(type==0){
            // for optimal coning compensation method
            for(i=0;i<nts_;i++){
                dvbm_[0]+=acces_[i][0];
                dvbm_[1]+=acces_[i][1];
                dvbm_[2]+=acces_[i][2];
            }
            MatrixXd sm(1,nts_-1),a(1,nts_-1),b(nts_-1,nts_-1);
            sm=MatrixXd::Zero(1,nts_-1);
            a=MatrixXd::Zero(1,nts_-1);b=MatrixXd::Zero(nts_-1,3);

            for(i=0;i<nts_-1;i++){
                a(0,i)=kConingCoeff[nts_-2][i];
            }
            for(i=0;i<nts_-1;i++){
                for(j=0;j<3;j++){
                    b(i,j)=acces_[i][j];
                }
            }
            sm=a*b;
            CrossVec3(cm.data(),acces_[nts_-1].data(),scull.data());
            CrossVec3(sm.data(),gyros_[nts_-1].data(),a.data());
            scull+=a;
        }
        else{
            for(i=0;i<nts_;i++){
                dvbm_[0]+=acces_[i][0];
                dvbm_[1]+=acces_[i][1];
                dvbm_[2]+=acces_[i][2];
            }
        }

        Vector3d rotm(0,0,0);
        CrossVec3(phim_.data(),dvbm_.data(),rotm.data());
//        dvbm_+=1.0/2.0*rotm+scull;
//        dtheta_=1.0/2.0*rotm+scull;
    }

    void cIns::InsUpdateEcef() {
        /// attitude update
        double alpha_ie=OMGE_GPS*dt_;
        Matrix3d C_earth;
        C_earth=Matrix3d::Identity();
        C_earth(0,0)=cos(alpha_ie);C_earth(0,1)=sin(alpha_ie);
        C_earth(1,0)=-sin(alpha_ie);C_earth(1,1)=cos(alpha_ie);

        // calculate attitude increment
        Vector3d alpha_ib_b=cur_imu_info_->cor_gyro_rate*dt_;
        double mag_alpha=alpha_ib_b.norm();
        Matrix3d Alpha_ib_b=VectorSkew(alpha_ib_b);

        // obtain coordinate transformation matrix from the new attitude w.r.t an inertial frame to the old using Rodrigues' formula
        Matrix3d C_new_old;
        if(mag_alpha>1E-8){
            C_new_old=Matrix3d::Identity()+sin(mag_alpha)/mag_alpha*Alpha_ib_b+(1-cos(mag_alpha))/SQR(mag_alpha)*Alpha_ib_b*Alpha_ib_b;
        }
        else{
            C_new_old=Matrix3d::Identity()+Alpha_ib_b;
        }

        cur_imu_info_->Cbe=C_earth*pre_imu_info_.Cbe*C_new_old;
        Dcm2Rpy(cur_imu_info_->Cbe,cur_imu_info_->rpy_e);

        /// calculate the average body-to-ECEF-frame coordinate transformation matrix over the update interval
        Vector3d w_ie(0,0,OMGE_GPS);
        Matrix3d avg_Cbe;
        if(mag_alpha>1.0E-8){
            Matrix3d C_avgb_oldb;
            C_avgb_oldb=Matrix3d::Identity()+(1.0-cos(mag_alpha))/SQR(mag_alpha)*Alpha_ib_b+(1.0-sin(mag_alpha)/mag_alpha)/SQR(mag_alpha)*Alpha_ib_b*Alpha_ib_b;
            avg_Cbe=pre_imu_info_.Cbe*C_avgb_oldb-0.5*VectorSkew(w_ie)*dt_*pre_imu_info_.Cbe;
        }
        else{
            avg_Cbe=pre_imu_info_.Cbe-0.5*VectorSkew(w_ie)*dt_*pre_imu_info_.Cbe;
        }

        /// specific force frame transformation
        Vector3d f_ib_e=avg_Cbe*cur_imu_info_->cor_acce_rate;

        /// update velocity
        cur_imu_info_->ve=pre_imu_info_.ve+(f_ib_e+CalculateGravity(pre_imu_info_.rn,true)-2*VectorSkew(w_ie)*pre_imu_info_.ve)*dt_;

        /// update position
        cur_imu_info_->re=pre_imu_info_.re+(cur_imu_info_->ve+pre_imu_info_.ve)*0.5*dt_;

        SolSync(*cur_imu_info_,COORD_XYZ);
    }

    /* Error Transition Matrix in ECEF */
    Eigen::MatrixXd cIns::StateTransferMat_E(tPPPLibConf C, tImuInfoUnit &pre_imu_info, tImuInfoUnit &cur_imu_info, int nx, double dt) {
        Matrix3d O33=Matrix3d::Zero();
        Matrix3d I33=Matrix3d::Identity();
        Vector3d w_ie(0,0,OMGE_GPS);

        Matrix3d Mpv=I33;
        Matrix3d Mvv=-2.0*VectorSkew(w_ie);

        Vector3d fe=cur_imu_info.Cbe*cur_imu_info.cor_acce_rate;
        Matrix3d Mva=-VectorSkew(fe);

        Matrix3d Maa=-1.0*VectorSkew(w_ie);

        int ip=0;
        int iv=3;
        int ia=6;
        int iba=9;
        int ibg=12;
        int isa=0,isg=0,ira=0,irg=0,ilev=0;
        if(C.insC.est_sa) isa=ibg+3;
        else isa=ibg;
        if(C.insC.est_sg) isg=isa+3;
        else isg=isa;
        if(C.insC.est_ra) ira=isg+3;
        else ira=isg;
        if(C.insC.est_rg) irg=ira+6;
        else irg=ira;
        if(C.insC.est_level) ilev=irg+6;
        else ilev=irg;

        MatrixXd F;
        F=MatrixXd::Zero(nx,nx);
        /*
         * F = [ O33 I33 O33  O33  O33
         *       O33 Mvv Mva -Cbe  O33
         *       O33 O33 Maa  O33  Cbe
         *       O33 O33 O33  Ta   O33
         *       O33 O33 O33  O33  Tg]
         * */

        F.block<3,3>(ip,iv)=Mpv;
        F.block<3,3>(iv,iv)=Mvv;
        F.block<3,3>(iv,ia)=Mva;
        F.block<3,3>(iv,iba)=cur_imu_info.Cbe;
        F.block<3,3>(ia,ia)=Maa;
        F.block<3,3>(ia,ibg)=cur_imu_info.Cbe;

        MatrixXd Fk=F*dt;
        if(nts_*ts_>0.1){

        }
        else{
            Fk+=MatrixXd::Identity(nx,nx);
        }
        return Fk;
    }
}