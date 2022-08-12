#include "su_rov.h"

SU_ROV::SU_ROV(QObject *parent) : QObject(parent)
{
    X_protocol = new x_protocol("kx_pult.conf", "x",X);
    K_protocol = new Qkx_coeffs("kx_pult.conf","k");

    X[1][0]=32;

    connect(&timer, &QTimer::timeout,[this](){
        X[2][0]=K[32];
    });
    timer.start(1000);
}
