#ifndef ROV_MODEL_H
#define ROV_MODEL_H

#include <QObject>

#define ANPA_MOD_CNT 23

extern QVector<double> K;
extern double X[2000][2];

struct InitData {
    double k_gamma;
    double m;
    double delta_m;
    double cv1[3];
    double cv2[3];
    double cw1[3];
    double cw2[3];
    double lambda[6][6];
    double J[3];
    double kd;
    double h;
    double Td;
    double l1,l2;
    double depth_limit;
    double max_depth;
}; //struct InitData


class ROV_Model : public QObject {
    Q_OBJECT
public:
    explicit ROV_Model(QObject *parent = 0);
    virtual ~ROV_Model();
signals:

public slots:
private:
    void start();

public:
    void model(const float Umvl,const float Umnl,const float Umvp,const float Umnp);
    void runge(const float Umvl,const float Umnl,
               const float Umvp,const float Umnp,const float Ttimer,const float dt=0.01);

    double a[ANPA_MOD_CNT];
    double da[ANPA_MOD_CNT];
    double delta_f;
    //константы
    double k_gamma;
    double m;
    double g;
    double G;
    double delta_m;
    double cv1[3];
    double cv2[3];
    double cw1[3];
    double cw2[3];
    double lambda[7][7];
    double J[3];
    double kd;
    double h;
    double Td;
    double l1,l2;
    double depth_limit;
    double max_depth;
    //переменные
    double sumX, sumZ;
    double cur_depth, Wx, Wy, Wz;
    double Psi_g, Gamma_g, Tetta_g;

    double Psi_gi, W_Psi_g, W_Gamma_g, W_Tetta_g;
    int N;
    double deltaSx, deltaSz;

    double Plz,Plp,Pmpr,Pml; //упоры движителей
    double Ppnp, Ppnl, Pznp, Pznl, Ppvp, Ppvl, Pzvl, Pzvp;
    double Ppnp_x, Ppnl_x, Pznp_x, Pznl_x, Ppvp_x, Ppvl_x, Pzvl_x, Pzvp_x;
    double Ppnp_y, Ppnl_y, Pznp_y, Pznl_y, Ppvp_y, Ppvl_y, Pzvl_y, Pzvp_y;
    double Ppnp_z, Ppnl_z, Pznp_z, Pznl_z, Ppvp_z, Ppvl_z, Pzvl_z, Pzvp_z;
    double Upnp, Upnl, Uznp, Uznl, Upvp, Upvl, Uzvl, Uzvp; //напряжения движителей

    double FloatageX, FloatageY, FloatageZ, Fdx, Fdy, Fdz, Fgx, Fgy, Fgz;
    double Mdx, Mdy, Mdz, Mgx, Mgy, Mgz;
    double Max,May,Maz; // моменты от силы Архимеда

    double x_global, y_global, z_global;
    double vx_local,  vy_local, vz_local;  //lineinye skorosti SPA v svyazannyh osyah
    double vx_global, vy_global, vz_global;

public:
    void resetModel();
    void tick(const float Umvl,const float Umnl,const float Umvp,const float Umnp,const float Ttimer);
    float Fx,Fy,Fz; //total forces for XYZ-axis
    float Mx,My,Mz; //total moments for XYZ-axis
};

#endif // ROV_MODEL_H
