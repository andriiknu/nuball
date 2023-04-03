#include "TMatrixD.h"
#include "TVector3.h"

TMatrixD GetTransMatrix (double theta, double phi) {

    double aphi[] = {cos(phi), sin(phi), 0, 
                    -sin(phi), cos(phi), 0, 
                    0, 0, 1};
    TMatrixD Aphi(3, 3, aphi);
    double atheta[] = {cos(theta), 0, -sin(theta),
                        0, 1, 0,
                        sin(theta), 0, cos(theta)};
    TMatrixD Atheta(3, 3, atheta);
    auto A = Atheta*Aphi;
    return A;
}

TMatrixD GetKetVector (double theta, double phi) {
    double c[] = {sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
    TMatrixD A(3, 1, c);
    return A;
}
