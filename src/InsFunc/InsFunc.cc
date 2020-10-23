//
// Created by cc on 7/16/20.
//

#include "InsFunc.h"

extern const double Crf[9]={0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,-1.0};

namespace PPPLib {
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
        rpy[0]=-atan2(m(1,2),m(2,2));
        rpy[1]=asin(m(0,2));
        rpy[2]=-atan2(m(0,1),m(0,0));
        return rpy;
    }

    Eigen::Vector3d Quaternion2Euler(const Quaterniond& q){
#if 0
        return RotationMatrix2Euler(q.toRotationMatrix());
#endif
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
        cur_imu_info.rpy=RotationMatrix2Euler(Cnb);
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
        if(is_ecef){
            const double constant_J2=0.00108263;
            const double constant_J4=-2.37091222e-6;
            const double constant_J6=6.08347e-9;
            double p = sqrt(coord_blh(0) * coord_blh(0) + coord_blh(1) * coord_blh(1) + coord_blh(2) * coord_blh(2));
            double t = coord_blh(2) / p;
            double a_p = WGS84_EARTH_LONG_RADIUS/ p;
            double a1 = -WGS84_GM / p / p;
            double a2 = 1 + 1.5 * constant_J2 * a_p * a_p - (15.0 / 8) * constant_J4 * pow(a_p,3) * a_p + (35.0 / 16) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double a3 = -4.5 * constant_J2 * a_p * a_p + (75.0 / 4) * constant_J4 * pow(a_p,3) * a_p - (735.0 / 16) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double a4 = -(175.0 / 8) * constant_J4 * pow(a_p,3) * a_p + (2205.0 / 16) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double a5 = -(1617.0 / 16) * constant_J6 * pow(a_p,3) * pow(a_p,3);

            double b1 = 3 * constant_J2 * a_p * a_p - (15.0 / 2) * constant_J4 * pow(a_p,3) * a_p + (105.0 / 8) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double b2 = (35.0 / 2) * constant_J4 * pow(a_p,3) * a_p - (945.0 / 12) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double b3 = (693.0 / 8) * constant_J6 * pow(a_p,3) * pow(a_p,3);

            double c1 = a2;
            double c2 = a3 - b1;
            double c3 = a4 - b2;
            double c4 = a5 - b3;
            double d1 = a2 + b1;
            double d2 = c2 + b2;
            double d3 = c3 + b3;
            double d4 = c4;
            Vector3d ge_vec;
            ge_vec(0) = (c1 + c2 * t * t + c3 * pow(t,3) * t + c4 * pow(t,3) * pow(t,3)) * coord_blh(0) * a1 / p + OMGE_GPS * OMGE_GPS * coord_blh(0);
            ge_vec(1) = (c1 + c2 * t * t + c3 * pow(t,3) * t + c4 * pow(t,3) * pow(t,3)) * coord_blh(1) * a1 / p + OMGE_GPS * OMGE_GPS * coord_blh(1);
            ge_vec(2) = (d1 + d2 * t * t + d3 * pow(t,3) * t + d4 * pow(t,3) * pow(t,3)) * coord_blh(2) * a1 / p;
            return ge_vec;
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
        LOG(DEBUG)<<"   "<<"ATTITUDE:(n)"<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.rpy.transpose()*R2D<<" deg";
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

    void AdjustImuData(tImuDataUnit& imu_data,IMU_COORD_TYPE coord_type,IMU_DATA_FORMAT data_format,GYRO_DATA_FORMAT gyro_val_format,double dt) {
        Vector3d gyro,acce;

        if(coord_type==IMU_COORD_RFU){
            gyro=imu_data.gyro;
            acce=imu_data.acce;
            MatMul("NN",3,1,3,1.0,Crf,gyro.data(),0.0,imu_data.gyro.data());
            MatMul("NN",3,1,3,1.0,Crf,acce.data(),0.0,imu_data.acce.data());
        }

        //
        if(data_format==IMU_FORMAT_RATE){
            for(int j=0;j<3;j++) imu_data.acce[j]*=dt;
            for(int j=0;j<3;j++) imu_data.gyro[j]*=dt;
        }
        if(gyro_val_format==GYRO_FORMAT_DEG){
            for(int j=0;j<3;j++) imu_data.gyro[j]*=D2R;
        }
    }

    cInsSim::cInsSim() {}

    cInsSim::cInsSim(double *ep) {
        sim_time_.Epoch2Time(ep);
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
        Vector3d blh(0,0,0);

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
//                if(line_num==1&&i==0) continue;
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
//        ImuErrSim(IMU_GRADE_INERTIAL);
    }

    /* IMU ERROR SETTING, refer to PSINS
     * gyro constant bias (deg/h) -->rad/s
     * acce constant bias (ug)
     * angular random walk (deg/sqrt(h))
     * velocity random walk (ug/sqrt(Hz))
     * */
    void cInsSim::ImuErrSim(PPPLib::IMU_GRADE grade) {
        Vector3d err_gyro,err_acce;

        if(grade==IMU_GRADE_INERTIAL){
            gyro_bias_<<0.03*ARC_DEG_PER_HOUR,0.03*ARC_DEG_PER_HOUR,0.03*ARC_DEG_PER_HOUR;      // ==> rad/s
            acce_bias_<<100*MILLI_G,100*MILLI_G,100*MILLI_G;                  // ==> m/s^2
            ang_rw_<<0.001*ARC_DEG_PER_SQRT_HOUR,0.001*ARC_DEG_PER_SQRT_HOUR,0.001*ARC_DEG_PER_SQRT_HOUR; // ==> rad/s^{1/2}
            vel_rw_<<5.0*UG_PER_SQRT_HZ,5.0*UG_PER_SQRT_HZ,5.0*UG_PER_SQRT_HZ;              // ==> m/s^{2.5}
        }
        else if(grade==IMU_GRADE_TACTICAL){

        }
        else if(grade==IMU_GRADE_MEMS){

        }

        double ts=1.0/sim_imu_data_.hz_;
        double sts=sqrt(1.0/sim_imu_data_.hz_);
        err_gyro[0]=ts*gyro_bias_[0]+sts*ang_rw_[0]*RandNorm(1.0);
        err_gyro[1]=ts*gyro_bias_[1]+sts*ang_rw_[1]*RandNorm(1.0);
        err_gyro[2]=ts*gyro_bias_[2]+sts*ang_rw_[2]*RandNorm(1.0);
        err_acce[0]=ts*acce_bias_[0]+sts*vel_rw_[0]*RandNorm(1.0);
        err_acce[1]=ts*acce_bias_[1]+sts*vel_rw_[1]*RandNorm(1.0);
        err_acce[2]=ts*acce_bias_[2]+sts*vel_rw_[2]*RandNorm(1.0);

        int i=0;
        for(i=0;i<sim_imu_data_.data_.size();i++){
            sim_imu_data_.data_.at(i).gyro+=err_gyro;
            sim_imu_data_.data_.at(i).acce+=err_acce;

            cout<<"GYRO SIM VALUE: "<<sim_imu_data_.data_.at(i).gyro.transpose()<<endl;
            cout<<"ACCE SIM VALUE: "<<sim_imu_data_.data_.at(i).acce.transpose()<<endl;
        }

        init_att_err_<<30*60,-30*60,20;
        init_vel_err_<<0.1,0.1,0.1;
        init_pos_err_<<1.0/WGS84_EARTH_LONG_RADIUS,1.0/WGS84_EARTH_LONG_RADIUS,3;
    }

    cIns::cIns() {}

    cIns::cIns(PPPLib::tPPPLibConf C) {
        C_=C;
        ts_=1.0/C_.insC.sample_rate;
        nts_=C_.insC.sample_number;
    }

    cIns::~cIns() {}

    void cIns::InitInsStat() {
        rpy0_<<-7.0513E-10,7.0513E-10,-1.6923E-06;
        vn0_<<0.1,0.1,0.1;
        rn0_<<0.5977,1.9008,383;

        UpdateEarthPar(rn0_,vn0_);
        Mpv_<<0,         1/eth_.RMh, 0,
              1/eth_.cos_lat_RNh,         0, 0,
              0,                  0, 1;
    }

    void cIns::UpdateVel_N() {
        cur_imu_info_->fn=Cbn_*cur_imu_info_->fb;
        Vector3d aaa=RotateVec(-eth_.w_n_in*nts_*ts_/2,cur_imu_info_->fn);
        cur_imu_info_->an=RotateVec(-eth_.w_n_in*nts_*ts_/2,cur_imu_info_->fn)+eth_.gcc_n;
        cur_imu_info_->vn=pre_imu_info_.vn+cur_imu_info_->an*nts_*ts_;
    }

    void cIns::UpdatePos_N() {
        Mpv_<<0,                   1/eth_.RMh, 0,
              1/eth_.cos_lat_RNh,  0,          0,
                0,                 0,          1;

        Vector3d dv=Mpv_*(cur_imu_info_->vn+pre_imu_info_.vn)/2.0;
        cur_imu_info_->rn=pre_imu_info_.rn+dv*nts_*ts_;
    }

    void cIns::UpdateAtt_N() {
        Qbn_=AttUpdateRotVec(Qbn_,phim_,eth_.w_n_in*ts_*nts_);
        cur_imu_info_->Qbn=Qbn_;
        cur_imu_info_->Cbn=Cbn_=Quaternion2RotationMatrix(Qbn_);
        cur_imu_info_->rpy=Quaternion2Euler(Qbn_);
    }

    Eigen::Quaterniond cIns::AttUpdateRotVec(Eigen::Quaterniond qnb, Vector3d rv_ib, Vector3d rv_in) {

        double n2=SQR(rv_ib.norm()),n,n_2;
        double rv_ib0=0.0,s=0.0;
        if(n2<1.0E-08){
            rv_ib0=1-n2*(1/8.0-n2/384.0);s=1.0/2-n2*(1.0/48-n2/3840.0);
        }
        else{
            n=sqrt(n2);n_2=n/2;
            rv_ib0=cos(n_2);s=sin(n_2)/n;
        }
        rv_ib*=s;
        double qb1,qb2,qb3,qb4;
        qb1=qnb.w()*rv_ib0-qnb.x()*rv_ib[0]-qnb.y()*rv_ib[1]-qnb.z()*rv_ib[2];
        qb2=qnb.w()*rv_ib[0]+qnb.x()*rv_ib0+qnb.y()*rv_ib[2]-qnb.z()*rv_ib[1];
        qb3=qnb.w()*rv_ib[1]+qnb.y()*rv_ib0+qnb.z()*rv_ib[0]-qnb.x()*rv_ib[2];
        qb4=qnb.w()*rv_ib[2]+qnb.z()*rv_ib0+qnb.x()*rv_ib[1]-qnb.y()*rv_ib[0];

        double rv_in0=0.0;
        n2=SQR(rv_in.norm());
        Eigen::Quaterniond qnb1;
        if(n2<1.0E-08){
            rv_in0=1-n2*(1.0/8-n2/384);s=-1.0/2+n2*(1.0/48-n2/3840);
        }
        else{
            n=sqrt(n2);n_2=n/2;
            rv_in0=cos(n_2);s=-sin(n_2)/n;
        }
        rv_in*=s;
        qnb1.w()=rv_in0*qb1-rv_in[0]*qb2-rv_in[1]*qb3-rv_in[2]*qb4;
        qnb1.x()=rv_in0*qb2+rv_in[0]*qb1+rv_in[1]*qb4-rv_in[2]*qb3;
        qnb1.y()=rv_in0*qb3+rv_in[1]*qb1+rv_in[2]*qb2-rv_in[0]*qb4;
        qnb1.z()=rv_in0*qb4+rv_in[2]*qb1+rv_in[0]*qb3-rv_in[1]*qb2;

        n2=SQR(qnb1.norm());
        if(n2>1.000001||n2<0.999999){
            double nq=1.0/sqrt(n2);
            qnb1.w()*=nq;
            qnb1.x()*=nq;
            qnb1.y()*=nq;
            qnb1.z()*=nq;
        }
        return qnb1;
    }

    void cIns::TraceInsMechInfo(tImuInfoUnit imu_info,bool prior,int idx) {
        LOG(DEBUG)<<"INS MECHANIZATION"<<(prior?"- ":"+ ")<< "("<<idx<<"): "<<imu_info.t_tag.GetTimeStr(4);
        LOG(DEBUG)<<"   "<<"GYRO VALUE: "<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.cor_gyro_rate.transpose()<<" rad/s";
        LOG(DEBUG)<<"   "<<"ACCE VALUE: "<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.cor_acce_rate.transpose()<<" m/s^2";
        LOG(DEBUG)<<"   "<<"ATTITUDE:(n)"<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.rpy.transpose()*R2D<<" deg";
        LOG(DEBUG)<<"   "<<"VELOCITY:(n)"<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.vn.transpose()<<" m/s";
        LOG(DEBUG)<<"   "<<"POSITION:(n)"<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.rn.transpose()<<" m";
        LOG(DEBUG)<<"   "<<"GYRO BIAS:  "<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.bg.transpose()/ARC_DEG_PER_HOUR<<" deg/h";
        LOG(DEBUG)<<"   "<<"ACCE BIAS:  "<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.ba.transpose()/MICRO_G<<" ug";
    }

    /* ins mech in ENU */
    void cIns::InsMech(PPPLib::tImuInfoUnit &cur_imu_info,const tImuInfoUnit pre_imu_info,int idx) {
        cur_imu_info_=&cur_imu_info;
        pre_imu_info_=pre_imu_info;

        TraceInsMechInfo(pre_imu_info_,true,idx);

        ConingScullingCompensation(0);
        phim_-=pre_imu_info_.bg*ts_*nts_;
        dvbm_-=pre_imu_info_.ba*ts_*nts_;
        cur_imu_info_->cor_gyro_rate=phim_/(ts_*nts_);
        cur_imu_info_->fb=cur_imu_info_->cor_acce_rate=dvbm_/(ts_*nts_);

        if(C_.insC.mech_coord==MECH_ENU){
            Vector3d vn_1=pre_imu_info_.vn+pre_imu_info_.an*(ts_*nts_/2.0);
            Vector3d rn_1=pre_imu_info_.rn+Mpv_*vn_1*(ts_*nts_/2.0);
            UpdateEarthPar(rn_1,vn_1);

            Qbn_=pre_imu_info_.Qbn;
            Cbn_=pre_imu_info_.Cbn;
            UpdateVel_N();
            UpdatePos_N();
            UpdateAtt_N();
        }
        else if(C_.insC.mech_coord==MECH_ECEF){
            Vector3d ve_1=pre_imu_info_.ve+pre_imu_info_.ae*(ts_*nts_/2.0);
            Vector3d re_1=pre_imu_info_.re+ve_1*(ts_*nts_/2.0);

            Qbe_=pre_imu_info_.Qbe;
            Cbe_=pre_imu_info_.Cbe;

            UpdateVel_E();
            UpdatePos_E();
            UpdateAtt_E();
        }

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

        Vector3d rotm(0,0,0);
        CrossVec3(phim_.data(),dvbm_.data(),rotm.data());
        dvbm_+=1.0/2.0*rotm+scull;
    }

    void cIns::UpdateEarthPar(Vector3d pos,Vector3d vel) {
        eth_.rn=pos;eth_.vn=vel;
        eth_.sin_lat=sin(pos[0]);eth_.cos_lat=cos(pos[0]);eth_.tan_lat=eth_.sin_lat/eth_.cos_lat;
        eth_.sin_lat2=SQR(eth_.sin_lat);
        eth_.sin_lat4=SQR(eth_.sin_lat2);
        eth_.sq=1.0-WGS84_FIRST_E2*eth_.sin_lat2;
        eth_.RN=WGS84_EARTH_LONG_RADIUS/sqrt(eth_.sq);
        eth_.RNh=eth_.RN+pos[2]; eth_.cos_lat_RNh=eth_.cos_lat*eth_.RNh;
        eth_.RMh=eth_.RN*(1.0-WGS84_FIRST_E2)/eth_.sq+pos[2];
        eth_.w_n_ie<<0.0,OMGE_GPS*eth_.cos_lat,OMGE_GPS*eth_.sin_lat;
        eth_.w_e_ie<<0.0,0.0,OMGE_GPS;
        eth_.w_n_en<<-vel[1]/eth_.RMh,vel[0]/eth_.RNh,vel[0]/eth_.RNh*eth_.tan_lat;
        eth_.w_n_in=eth_.w_n_ie+eth_.w_n_en;
        eth_.w_n_ie_n=eth_.w_n_ie+eth_.w_n_in;
        eth_.g=NORM_GRAVITY_VALUE*(1+5.27094E-03*eth_.sin_lat2+2.32718E-05*eth_.sin_lat4)-3.086E-06*pos[2];
        eth_.gn<<0,0,-eth_.g;
        eth_.ge=CalculateGravity(pos,true);
        eth_.gcc_n[0]=eth_.w_n_ie_n[2]*vel[1]-eth_.w_n_ie_n[1]*vel[2];
        eth_.gcc_n[1]=eth_.w_n_ie_n[0]*vel[2]-eth_.w_n_ie_n[2]*vel[0];
        eth_.gcc_n[2]=eth_.w_n_ie_n[1]*vel[0]-eth_.w_n_ie_n[0]*vel[1]+eth_.gn[2];
        eth_.gcc_e=-2.0*VectorSkew(eth_.w_e_ie)*vel+eth_.ge;
    }

    /* Error Transition Matrix in ENU
     * Ft = [Mpp Mpv O33 O33 O33
     *       Mvp Mvv Mva -Cbn O33
     *       Map Mav Maa O33 Cbn
     *       O33 O33 O33 Ta  O33
     *       O33 O33 O33 O33 Tg]
     * */
    Eigen::MatrixXd cIns::StateTransferMat_N(tPPPLibConf C,tImuInfoUnit& pre_imu_info,tImuInfoUnit& cur_imu_info,int nx,double dt) {
        double tan_lat=eth_.tan_lat,sec_lat=1.0/eth_.cos_lat;
        double f_RMh=1.0/eth_.RMh, f_RNh=1.0/eth_.RNh,f_clRNh=1.0/eth_.cos_lat_RNh;
        double f_RMh2=f_RMh*f_RMh, f_RNh2=f_RNh*f_RNh;
        Vector3d vn=cur_imu_info.vn;
        double ve_cl_RNh=vn[0]*f_clRNh,ve_RNh2=vn[0]*f_RNh2,vn_RMh2=vn[1]*f_RMh2;
        Matrix3d O33=Matrix3d::Zero();

        Matrix3d Mp1,Mp2;
        Mp1=Matrix3d::Zero();Mp2=Matrix3d::Zero();
        Mp1(1,0)=-eth_.w_n_ie[2];Mp1(2,0)=eth_.w_n_ie[1];
        Mp2(2,0)=ve_cl_RNh*sec_lat;
        Mp2(0,2)=vn_RMh2;
        Mp2(1,2)=-ve_RNh2;
        Mp2(2,2)=-ve_RNh2*tan_lat;

        Matrix3d Avn=VectorSkew(vn);
        Matrix3d Awn=VectorSkew(eth_.w_n_ie_n);
        Matrix3d Maa;
        Maa=Matrix3d::Zero();
        Maa(1,0)=-eth_.w_n_in[2];
        Maa(2,0)=eth_.w_n_in[1];
        Maa(0,1)=eth_.w_n_in[2];
        Maa(2,1)=-eth_.w_n_in[0];
        Maa(0,2)=-eth_.w_n_in[1];
        Maa(1,2)=eth_.w_n_in[0];

        Matrix3d Mav,Map;
        Mav=Matrix3d::Zero();
        Mav(1,0)=f_RNh;
        Mav(2,0)=f_RNh*tan_lat;
        Mav(0,1)=-f_RNh;

        Map=Mp1+Mp2;

        Matrix3d Mva,Mvv,Mvp;
        Mva=Matrix3d::Zero();
        Vector3d fn=Cbn_*cur_imu_info.cor_acce_rate;
        Mva=VectorSkew(fn);
        Mvv=Avn*Mav-Awn;
        Mvp=Avn*(Mp1+Map);

        double g0=9.7803267714;
        double sc_lat=eth_.sin_lat*eth_.cos_lat;
        Mvp(2,0)-=g0*(5.27094e-3*2+2.32718e-5*4*eth_.sin_lat2)*sc_lat;
        Mvp(2,2)+=3.086e-6;

        Matrix3d Mpv=Mpv_,Mpp;
        Mpp=Matrix3d::Zero();
        Mpp(1,0)=ve_cl_RNh*tan_lat;
        Mpp(0,2)=-vn_RMh2;
        Mpp(1,2)=-ve_RNh2*sec_lat;

        MatrixXd F;
        F=MatrixXd::Zero(nx,nx);

        /*
         * Ft=[ Mpp Mpv O33  O33  O33
         *      Mvp Mvv Mva -Cbn  O33
         *      Map Mav Maa  O33  Cbn
         *      O33 O33 O33  Ta   O33
         *      O33 O33 O33  O33  Tg ]
         * */
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

        F.block<3,3>(ip,ip)=Mpp;
        F.block<3,3>(ip,iv)=Mpv;

        F.block<3,3>(iv,ip)=Mvp;
        F.block<3,3>(iv,iv)=Mvv;
        F.block<3,3>(iv,ia)=Mva;
        F.block<3,3>(iv,iba)=-cur_imu_info.Cbn;

        F.block<3,3>(ia,ip)=Map;
        F.block<3,3>(ia,iv)=Mav;
        F.block<3,3>(ia,ia)=Maa;
        F.block<3,3>(ia,ibg)=cur_imu_info.Cbn;

        MatrixXd Fk=F*(nts_*ts_);
        if(nts_*ts_>0.1){

        }
        else{
            Fk+=MatrixXd::Identity(nx,nx);
        }

        return Fk;
    }

    void cIns::UpdateVel_E() {
        Vector3d fe=Cbe_*cur_imu_info_->cor_acce_rate;
        cur_imu_info_->ae=RotateVec(-eth_.w_e_ie*nts_*ts_/2,fe)+eth_.gcc_e;
        cur_imu_info_->ve=pre_imu_info_.ve+cur_imu_info_->ae*nts_*ts_;
    }

    void cIns::UpdatePos_E() {
        cur_imu_info_->re=pre_imu_info_.re+(pre_imu_info_.ve+cur_imu_info_->ve)*nts_*ts_/2.0;
    }

    void cIns::UpdateAtt_E() {
        Qbe_=AttUpdateRotVec(Qbe_,phim_,eth_.w_e_ie*ts_*nts_);
        cur_imu_info_->Qbe=Qbe_;
        cur_imu_info_->Cbe=Cbe_=Quaternion2RotationMatrix(Qbe_);
        cur_imu_info_->rpy=Quaternion2Euler(Qbe_);
    }

    /* Error Transition Matrix in ECEF */
    Eigen::MatrixXd cIns::StateTransferMat_E(tPPPLibConf C, tImuInfoUnit &pre_imu_info, tImuInfoUnit &cur_imu_info, int nx, double dt) {
        Matrix3d O33=Matrix3d::Zero();
        Matrix3d I33=Matrix3d::Identity();

        Matrix3d Mpv=I33;
        Matrix3d Mvv=-2.0*VectorSkew(eth_.w_e_ie);

        Vector3d fe=Cbe_*cur_imu_info.cor_acce_rate;
        Matrix3d Mva=VectorSkew(fe);

        Matrix3d Maa=-VectorSkew(eth_.w_e_ie);

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
        F.block<3,3>(iv,iba)=-cur_imu_info.Cbe;
        F.block<3,3>(ia,ia)=Maa;

        MatrixXd Fk=F*(nts_*ts_);
        if(nts_*ts_>0.1){

        }
        else{
            Fk+=MatrixXd::Identity(nx,nx);
        }

        return Fk;
    }


}