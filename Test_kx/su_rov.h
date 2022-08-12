#ifndef SU_ROV_H
#define SU_ROV_H

#include <QObject>
#include <QVector>
#include "qkx_coeffs.h"
#include "kx_protocol.h"
#include <QTimer>
extern double X[2000][2];
extern QVector<double> K;

class SU_ROV : public QObject
{
    Q_OBJECT
public:
    explicit SU_ROV(QObject *parent = 0);

signals:

public slots:
protected:
    Qkx_coeffs * K_protocol=nullptr;
    x_protocol * X_protocol=nullptr;
    QTimer timer;
};

#endif // SU_ROV_H
