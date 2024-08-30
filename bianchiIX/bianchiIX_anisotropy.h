//
// Created by rafkat on 8/26/24.
//

#ifndef ODECALCULATIONS_BIANCHIIX_ANISOTROPY_H
#define ODECALCULATIONS_BIANCHIIX_ANISOTROPY_H

#include "../constants.h"
#include "../nr/nr3.h"
#include "../nr/odeint.h"
#include <cmath>


struct bianchiIX_anisotropic_rhs {
    Doub s_eta;
    Doub s_omega_2;

    bianchiIX_anisotropic_rhs(Doub eta, Doub omega_2) : s_eta(eta), s_omega_2(omega_2) {}

    void operator()(const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        Doub numeratorDyDx1 = (
                (-432.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                 - 72.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                 + 9.0 * s_eta) * std::pow(y[0], 4.0) * std::pow(y[3], 4.0)
                + (-96.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   + 36.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2])
                   + 3.0) * std::pow(y[0], 4.0) * std::pow(y[3], 2.0)
                + (16.0 * s_eta * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   + 4.0 * Pi * std::pow(EulerConst, 2.0 * y[2])) * std::pow(y[0], 4.0)
                + (-144.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   - 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                   + 3.0 * s_eta) * std::pow(y[0], 3.0) * std::pow(y[1], 2.0) * std::pow(y[3], 2.0)
                + (16.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   + 8.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2])
                   + 1.0) * std::pow(y[0], 3.0) * std::pow(y[1], 2.0)
                + (432.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   + 72.0 * std::pow(s_eta, 2) * Pi * std::pow(EulerConst, 2.0 * y[2])
                   - 9.0 * s_eta) * std::pow(y[0], 2.0) * std::pow(y[1], 2.0) * std::pow(y[3], 2.0)
                + (-240.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   - 60.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2])) * std::pow(y[0], 2.0) * std::pow(y[1], 2.0)
                + (-120.0 * Pi * s_omega_2 * std::pow(s_eta, 2.0) * std::pow(EulerConst, 2.0 * y[2])
                   - 6.0 * s_eta * s_omega_2) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                + (-64.0 * std::pow(s_eta, 2.0) * s_omega_2 * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   - 12.0 * s_omega_2 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2])
                   + s_omega_2) * std::pow(y[0], 2.0)
                + (-48.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   - 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                   - 3.0 * s_eta) * y[0] * std::pow(y[1], 4.0)
                + (-48.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   - 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                   - 3.0 * s_eta) * s_omega_2 * y[0] * std::pow(y[1], 2.0)
                + (192.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   + 48.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])) * std::pow(y[1], 4.0)
                + (240.0 * std::pow(Pi, 2.0) * std::pow(s_eta, 3.0) * s_omega_2 * std::pow(EulerConst, 4.0 * y[2])
                   + 48.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                   - 3.0 * s_omega_2 * s_eta) * std::pow(y[1], 2.0)
                + 48.0 * std::pow(EulerConst, 4.0 * y[2]) * std::pow(Pi, 2.0) * std::pow(s_omega_2, 2.0) *
                  std::pow(s_eta, 3.0)
                - 3.0 * std::pow(s_omega_2, 2.0) * s_eta
        );

        Doub denominatorDyDx1 = (
                96.0 * (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0 - 1.0 / 12.0) * std::pow(y[0], 2)
                        + (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                           + s_eta / 4.0) * std::pow(y[1], 2.0)
                        + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                        + s_omega_2 * s_eta / 4.0
                ) * y[0] * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
        );

        Doub numeratorDyDx2 = (
                -2.0 * y[1] * (
                        (-s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]) - 1.0 / 8.0) * std::pow(y[0], 2.0)
                        + (-std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0 - s_eta / 8.0) * y[0]
                          * std::pow(y[1], 2.0)
                        + (std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0 + s_eta / 8.0)
                          * std::pow(y[1], 2.0)
                        + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                )
        );

        Doub denominatorDyDx2 = (
                y[0] * (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0 - 1.0 / 12.0) * std::pow(y[0], 2.0)
                        + (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) + s_eta / 4.0)
                          * std::pow(y[1], 2.0)
                        + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                        + s_omega_2 * s_eta / 4.0
                )
        );


        Doub numeratorDyDx3 = (
                y[1] * y[3] * (
                        (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                         - 9.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                         + 3.0 * s_eta / 16.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (-3.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                           + 1.0 / 16.0) * std::pow(y[0], 2.0)
                        + (-std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                           - 2.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]))
                          * std::pow(y[1], 2.0) * y[0]
                        + (2.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                           + 11.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                           - 3.0 * s_eta / 16.0) * std::pow(y[1], 2.0)
                        - 3.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0) / 2.0
                        + std::pow(EulerConst, 4.0 * y[2]) * std::pow(Pi, 2.0) * s_omega_2 * std::pow(s_eta, 3.0)
                        - 3.0 * s_omega_2 * s_eta / 16.0
                )
        );


        Doub denominatorDyDx3 = (
                (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0 - 1.0 / 12.0) * std::pow(y[0], 2)
                        + (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                           + s_eta / 4.0) * std::pow(y[1], 2.0)
                        + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                        + s_omega_2 * s_eta / 4.0
                ) * y[0] * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
        );

        dydx[0] = y[1];
        dydx[1] = numeratorDyDx1 / denominatorDyDx1;
        dydx[2] = numeratorDyDx2 / denominatorDyDx2;
        dydx[3] = numeratorDyDx3 / denominatorDyDx3;
    }

    void jacobian(const Doub x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy) {
        Doub numeratorDyDx1 = (
                (-432.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                 - 72.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                 + 9.0 * s_eta) * std::pow(y[0], 4.0) * std::pow(y[3], 4.0)
                + (-96.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   + 36.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2])
                   + 3.0) * std::pow(y[0], 4.0) * std::pow(y[3], 2.0)
                + (16.0 * s_eta * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   + 4.0 * Pi * std::pow(EulerConst, 2.0 * y[2])) * std::pow(y[0], 4.0)
                + (-144.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   - 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                   + 3.0 * s_eta) * std::pow(y[0], 3.0) * std::pow(y[1], 2.0) * std::pow(y[3], 2.0)
                + (16.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   + 8.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2])
                   + 1.0) * std::pow(y[0], 3.0) * std::pow(y[1], 2.0)
                + (432.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   + 72.0 * std::pow(s_eta, 2) * Pi * std::pow(EulerConst, 2.0 * y[2])
                   - 9.0 * s_eta) * std::pow(y[0], 2.0) * std::pow(y[1], 2.0) * std::pow(y[3], 2.0)
                + (-240.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   - 60.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2])) * std::pow(y[0], 2.0) * std::pow(y[1], 2.0)
                + (-120.0 * Pi * s_omega_2 * std::pow(s_eta, 2.0) * std::pow(EulerConst, 2.0 * y[2])
                   - 6.0 * s_eta * s_omega_2) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                + (-64.0 * std::pow(s_eta, 2.0) * s_omega_2 * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   - 12.0 * s_omega_2 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2])
                   + s_omega_2) * std::pow(y[0], 2.0)
                + (-48.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   - 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                   - 3.0 * s_eta) * y[0] * std::pow(y[1], 4.0)
                + (-48.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   - 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                   - 3.0 * s_eta) * s_omega_2 * y[0] * std::pow(y[1], 2.0)
                + (192.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   + 48.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])) * std::pow(y[1], 4.0)
                + (240.0 * std::pow(Pi, 2.0) * std::pow(s_eta, 3.0) * s_omega_2 * std::pow(EulerConst, 4.0 * y[2])
                   + 48.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                   - 3.0 * s_omega_2 * s_eta) * std::pow(y[1], 2.0)
                + 48.0 * std::pow(EulerConst, 4.0 * y[2]) * std::pow(Pi, 2.0) * std::pow(s_omega_2, 2.0) *
                  std::pow(s_eta, 3.0)
                - 3.0 * std::pow(s_omega_2, 2.0) * s_eta
        );

        Doub denominatorDyDx1 = (
                96.0 * (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0 - 1.0 / 12.0) * std::pow(y[0], 2)
                        + (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                           + s_eta / 4.0) * std::pow(y[1], 2.0)
                        + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                        + s_omega_2 * s_eta / 4.0
                ) * y[0] * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
        );

        Doub numeratorDyDx2 = (
                -2.0 * y[1] * (
                        (-s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]) - 1.0 / 8.0) * std::pow(y[0], 2.0)
                        + (-std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0 - s_eta / 8.0) * y[0]
                          * std::pow(y[1], 2.0)
                        + (std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0 + s_eta / 8.0)
                          * std::pow(y[1], 2.0)
                        + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                )
        );

        Doub denominatorDyDx2 = (
                y[0] * (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0 - 1.0 / 12.0) * std::pow(y[0], 2.0)
                        + (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) + s_eta / 4.0)
                          * std::pow(y[1], 2.0)
                        + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                        + s_omega_2 * s_eta / 4.0
                )
        );


        Doub numeratorDyDx3 = (
                y[1] * y[3] * (
                        (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                         - 9.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                         + 3.0 * s_eta / 16.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (-3.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                           + 1.0 / 16.0) * std::pow(y[0], 2.0)
                        + (-std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                           - 2.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]))
                          * std::pow(y[1], 2.0) * y[0]
                        + (2.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                           + 11.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]))
                          * std::pow(y[1], 2.0)
                        - 3.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0) / 2.0
                        + std::pow(EulerConst, 4.0 * y[2]) * std::pow(Pi, 2.0) * s_omega_2 * std::pow(s_eta, 3.0)
                        - 3.0 * s_omega_2 * s_eta / 16.0
                )
        );


        Doub denominatorDyDx3 = (
                (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0 - 1.0 / 12.0) * std::pow(y[0], 2)
                        + (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                           + s_eta / 4.0) * std::pow(y[1], 2.0)
                        + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                        + s_omega_2 * s_eta / 4.0
                ) * y[0] * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
        );

        Int n = y.size();
        for (Int i = 0; i < n; i++) dfdx[i] = 0.0;
        dfdy[0][0] = 0.0;
        dfdy[0][1] = 1.0;
        dfdy[0][2] = 0.0;
        dfdy[0][3] = 0.0;

        dfdy[1][0] = (
                             (-432.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                              - 72.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                              + 9.0 * s_eta) * 4.0 * std::pow(y[0], 3.0) * std::pow(y[3], 4.0)
                             + (-96.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                + 36.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2])
                                + 3.0) * 4.0 * std::pow(y[0], 3.0) * std::pow(y[3], 2.0)
                             + (16.0 * s_eta * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                + 4.0 * Pi * std::pow(EulerConst, 2.0 * y[2]))
                               * 4.0 * std::pow(y[0], 3.0)
                             + (-144.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                + 3.0 * s_eta)
                               * 3.0 * std::pow(y[0], 2.0) * std::pow(y[1], 2.0) * std::pow(y[3], 2.0)
                             + (16.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                + 8.0 * s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi
                                + 1.0) * 3.0 * std::pow(y[0], 2.0) * std::pow(y[1], 2.0)
                             + (432.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                + 72.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                - 9.0 * s_eta) * 2.0 * y[0] * std::pow(y[1], 2.0) * std::pow(y[3], 2.0)
                             + (-240.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 60.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]))
                               * 2.0 * y[0] * std::pow(y[1], 2.0)
                             + (-120.0 * Pi * s_omega_2 * std::pow(s_eta, 2.0) * std::pow(EulerConst, 2.0 * y[2])
                                - 6.0 * s_omega_2 * s_eta) * 2.0 * y[0] * std::pow(y[3], 2.0)
                             + (-64.0 * std::pow(s_eta, 2.0) * s_omega_2 * std::pow(Pi, 2.0) *
                                std::pow(EulerConst, 4.0 * y[2])
                                - 12.0 * s_omega_2 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2])
                                + s_omega_2) * 2.0 * y[0]
                             + (-48.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                - 3.0 * s_eta) * std::pow(y[1], 4.0)
                             + (-48.0 * std::pow(EulerConst, 4.0 * y[2]) * std::pow(Pi, 2.0) * s_omega_2 *
                                std::pow(s_eta, 3.0)
                                - 24.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                                - 3.0 * s_omega_2 * s_eta) * std::pow(y[1], 2.0)
                     ) / denominatorDyDx1
                     + (-1.0) * numeratorDyDx1 / std::pow(denominatorDyDx1, 2.0) * (
                96.0 * (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * 2.0 * y[0] * std::pow(y[3], 2.0)
                        + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0 - 1.0 / 12.0) * 2.0 * y[0]
                ) * y[0] * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
                + 96.0 * (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0 - 1.0 / 12.0) * std::pow(y[0], 2.0)
                        + (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) + s_eta / 4.0)
                          * std::pow(y[1], 2.0)
                        + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                        + s_omega_2 * s_eta / 4.0
                ) * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
        );

        dfdy[1][1] = (
                             (-144.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                              - 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                              + 3.0 * s_eta) * 2.0 * std::pow(y[0], 3.0) * y[1] * std::pow(y[3], 2.0)
                             + (16.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                + 8.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2])
                                + 1.0) * 2.0 * std::pow(y[0], 3.0) * y[1]
                             + (432.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                + 72.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                - 9.0 * s_eta) * 2.0 * std::pow(y[0], 2.0) * y[1] * std::pow(y[3], 2.0)
                             + (-240.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 60.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]))
                               * 2.0 * std::pow(y[0], 2.0) * y[1]
                             + (-48.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                - 3.0 * s_eta) * 4.0 * y[0] * std::pow(y[1], 3.0)
                             + (-48.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                - 3.0 * s_eta) * s_omega_2 * 2.0 * y[0] * y[1]
                             + (192.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                + 48.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]))
                               * 4.0 * std::pow(y[1], 3.0)
                             + (240.0 * std::pow(EulerConst, 4.0 * y[2]) * std::pow(Pi, 2.0) * s_omega_2 *
                                std::pow(s_eta, 3.0)
                                + 48.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                                - 3.0 * s_omega_2 * s_eta) * 2.0 * y[1]
                     ) / denominatorDyDx1
                     + (-1.0) * numeratorDyDx1 / std::pow(denominatorDyDx1, 2.0) * (
                96.0 * (
                        (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) + s_eta / 4.0) * 2.0 * y[1]
                ) * y[0] * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
        );

        dfdy[1][2] = (
                             (-4.0 * 432.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                              - 2.0 * 72.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]))
                             * std::pow(y[0], 4.0) * std::pow(y[3], 4.0)
                             + (-4.0 * 96.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0)
                                * std::pow(EulerConst, 4.0 * y[2])
                                + 2.0 * 36.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]))
                               * std::pow(y[0], 4.0) * std::pow(y[3], 2.0)
                             + (4.0 * 16.0 * s_eta * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                + 2.0 * 4.0 * Pi * std::pow(EulerConst, 2.0 * y[2])) * std::pow(y[0], 4.0)
                             + (-4.0 * 144.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0)
                                * std::pow(EulerConst, 4.0 * y[2])
                                - 2.0 * 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]))
                               * std::pow(y[0], 3.0) * std::pow(y[1], 2.0) * std::pow(y[3], 2.0)
                             + (4.0 * 16.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                + 2.0 * 8.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]))
                               * std::pow(y[0], 3.0) * std::pow(y[1], 2.0)
                             + (4.0 * 432.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0)
                                * std::pow(EulerConst, 4.0 * y[2])
                                + 2.0 * 72.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]))
                               * std::pow(y[0], 2.0) * std::pow(y[1], 2.0) * std::pow(y[3], 2.0)
                             + (-4.0 * 240.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0)
                                * std::pow(EulerConst, 4.0 * y[2])
                                - 2.0 * 60.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]))
                               * std::pow(y[0], 2.0) * std::pow(y[1], 2.0)
                             + (-120.0 * 2.0 * Pi * s_omega_2 * std::pow(s_eta, 2.0) * std::pow(EulerConst, 2.0 * y[2]))
                               * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                             + (-4.0 * 64.0 * std::pow(s_eta, 2.0) * s_omega_2 * std::pow(Pi, 2.0)
                                * std::pow(EulerConst, 4.0 * y[2])
                                - 2.0 * 12.0 * s_omega_2 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]))
                               * std::pow(y[0], 2.0)
                             + (-4.0 * 48.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0)
                                * std::pow(EulerConst, 4.0 * y[2])
                                - 2.0 * 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]))
                               * y[0] * std::pow(y[1], 4.0)
                             + (-4.0 * 48.0 * std::pow(EulerConst, 4.0 * y[2]) * std::pow(Pi, 2.0)
                                * s_omega_2 * std::pow(s_eta, 3.0)
                                - 2.0 * 24.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0))
                               * y[0] * std::pow(y[1], 2.0)
                             + (4.0 * 192.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0)
                                * std::pow(EulerConst, 4.0 * y[2])
                                + 2.0 * 48.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]))
                               * std::pow(y[1], 4.0)
                             + (4.0 * 240.0 * std::pow(EulerConst, 4.0 * y[2]) * std::pow(Pi, 2.0)
                                * s_omega_2 * std::pow(s_eta, 3.0)
                                + 2.0 * 48.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0))
                               * std::pow(y[1], 2.0)
                             + 4.0 * 48.0 * std::pow(EulerConst, 4.0 * y[2]) * std::pow(Pi, 2.0)
                               * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0)
                     ) / denominatorDyDx1
                     + (-1.0) * numeratorDyDx1 / std::pow(denominatorDyDx1, 2.0) * (
                96.0 * (
                        3.0 * 2.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[0], 2.0)
                        * std::pow(y[3], 2.0)
                        - 2.0 * s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi * std::pow(y[0], 2.0) / 3.0
                        - 3.0 * 2.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 2.0)
                        + 2.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                ) * y[0] * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
                + 96.0 * (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0 - 1.0 / 12.0) * std::pow(y[0], 2.0)
                        + (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                           + s_eta / 4.0) * std::pow(y[1], 2.0)
                        + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                        + s_omega_2 * s_eta / 4.0) * y[0] * 2.0 * s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi
        );

        dfdy[1][3] = (
                             (-432.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                              - 72.0 * Pi * std::pow(s_eta, 2.0) * std::pow(EulerConst, 2.0 * y[2])
                              + 9.0 * s_eta) * 4.0 * std::pow(y[0], 4.0) * std::pow(y[3], 3.0)
                             + (-96.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                + 36.0 * s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi
                                + 3.0) * 2.0 * std::pow(y[0], 4.0) * y[3]
                             + (-144.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 24.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                + 3.0 * s_eta) * 2.0 * std::pow(y[0], 3.0) * std::pow(y[1], 2.0) * y[3]
                             + (-120.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                                - 6.0 * s_omega_2 * s_eta) * 2.0 * std::pow(y[0], 2.0) * y[3]
                             + (432.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                + 72.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                - 9.0 * s_eta) * 2.0 * std::pow(y[0], 2.0) * std::pow(y[1], 2.0) * y[3]
                     ) / denominatorDyDx1
                     + (-1.0) * numeratorDyDx1 / std::pow(denominatorDyDx1, 2.0) * (
                96.0 * (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * 2.0 * std::pow(y[0], 2.0) * y[3]
                ) * y[0] * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
        );


        dfdy[2][0] = -2.0 * y[1] * (
                (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi - 1.0 / 8.0) * 2.0 * y[0]
                + (-std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0 - s_eta / 8.0)
                  * std::pow(y[1], 2.0)
        ) / denominatorDyDx2
                     + (-1.0) * numeratorDyDx2 / std::pow(denominatorDyDx2, 2.0) * (
                (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0 - 1.0 / 12.0) * std::pow(y[0], 2.0)
                        + (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                           + s_eta / 4.0) * std::pow(y[1], 2.0)
                        + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                        + s_omega_2 * s_eta / 4.0
                ) + y[0] * (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * 2.0 * y[0] * std::pow(y[3], 2.0)
                        + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0
                           - 1.0 / 12.0) * 2.0 * y[0]
                )
        );

        dfdy[2][1] = (
                             -2.0 * (
                                     (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi - 1.0 / 8.0) * std::pow(y[0], 2.0)
                                     + (-std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                                        - s_eta / 8.0) * std::pow(y[1], 2.0) * y[0]
                                     + (std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                                        + s_eta / 8.0) * std::pow(y[1], 2.0)
                                     + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                             ) - 2.0 * y[1] * (
                                     (-std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                                      - s_eta / 8.0) * 2.0 * y[1] * y[0]
                                     + (std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                                        + s_eta / 8.0) * 2.0 * y[1]
                             )
                     ) / denominatorDyDx2
                     + (-1.0) * numeratorDyDx2 / std::pow(denominatorDyDx2, 2) * (
                y[0] * (
                        -3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) + s_eta / 4.0
                ) * 2.0 * y[1]
        );

        dfdy[2][2] = -2.0 * y[1] * std::pow(EulerConst, 2.0 * y[2]) * (
                -2.0 * s_eta * Pi * std::pow(y[0], 2.0)
                - std::pow(s_eta, 2.0) * Pi * std::pow(y[1], 2.0) * y[0]
                + std::pow(s_eta, 2.0) * Pi * std::pow(y[1], 2.0)
                + 2.0 * Pi * s_omega_2 * std::pow(s_eta, 2.0)
        ) / denominatorDyDx2
                     + (-1.0) * numeratorDyDx2 / std::pow(denominatorDyDx2, 2.0) * (
                y[0] * std::pow(EulerConst, 2.0 * y[2]) * (
                        2.0 * 3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        - 2.0 * s_eta * Pi * std::pow(y[0], 2.0) / 3.0
                        - 2.0 * 3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(y[1], 2.0)
                        + 2.0 * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                )
        );

        dfdy[2][3] = (-1.0) * numeratorDyDx2 / std::pow(denominatorDyDx2, 2.0) * (
                y[0] * (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * 2.0 * std::pow(y[0], 2.0) * y[3]
                )
        );

        dfdy[3][0] = (
                             y[1] * y[3] * (
                                     (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                                      -
                                      9.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                      + 3.0 * s_eta / 16.0) * 2.0 * y[0] * std::pow(y[3], 2.0)
                                     +
                                     (-3.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                      + 1.0 / 16.0) * 2.0 * y[0]
                                     + (-std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                                        - 2.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) *
                                          std::pow(EulerConst, 4.0 * y[2]))
                                       * std::pow(y[1], 2.0)
                             )
                     ) / denominatorDyDx3
                     + (-1.0) * numeratorDyDx3 / std::pow(denominatorDyDx3, 2.0) * (
                (
                        (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                         - s_eta / 4.0) * 2.0 * y[0] * std::pow(y[3], 2.0)
                        + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0 - 1.0 / 12.0) * 2.0 * y[0]
                ) * y[0] * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
                + (
                          (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                           - s_eta / 4.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                          + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0 - 1.0 / 12.0) * std::pow(y[0], 2.0)
                          + (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                             + s_eta / 4.0) * std::pow(y[1], 2.0)
                          + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                          + s_omega_2 * s_eta / 4.0
                  ) * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
        );

        dfdy[3][1] = (
                             y[3] * (
                                     (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                                      -
                                      9.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                      + 3.0 * s_eta / 16.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                                     +
                                     (-3.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                      + 1.0 / 16.0) * std::pow(y[0], 2.0)
                                     + (-std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                                        - 2.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) *
                                          std::pow(EulerConst, 4.0 * y[2]))
                                       * std::pow(y[1], 2.0) * y[0]
                                     + (2.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                        + 11.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) *
                                          std::pow(EulerConst, 4.0 * y[2])
                                        - 3.0 * s_eta / 16.0) * std::pow(y[1], 2.0)
                                     - 3.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0) /
                                       2.0
                                     + std::pow(EulerConst, 4.0 * y[2]) * std::pow(Pi, 2.0) * s_omega_2 *
                                       std::pow(s_eta, 3.0)
                                     - 3.0 * s_omega_2 * s_eta / 16.0
                             ) + y[1] * y[3] * (
                                     (-std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                                      -
                                      2.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]))
                                     * 2.0 * y[1] * y[0]
                                     + (2.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                        + 11.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) *
                                          std::pow(EulerConst, 4.0 * y[2])
                                        - 3.0 * s_eta / 16.0) * 2.0 * y[1]
                             )
                     ) / denominatorDyDx3
                     + (-1.0) * numeratorDyDx3 / std::pow(denominatorDyDx3, 2.0) * (
                2.0 * y[1] * (
                        -3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                        + s_eta / 4.0
                ) * y[0] * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
        );

        dfdy[3][2] = (
                             y[1] * y[3] * (
                                     (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                      - 4.0 * 9.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) *
                                        std::pow(EulerConst, 4.0 * y[2]))
                                     * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                                     - 4.0 * 3.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) *
                                       std::pow(EulerConst, 4.0 * y[2])
                                       * std::pow(y[0], 2.0)
                                     + (-std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                        - 4.0 * 2.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) *
                                          std::pow(EulerConst, 4.0 * y[2]))
                                       * std::pow(y[1], 2.0) * y[0]
                                     + (2.0 * 2.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                        + 11.0 * 4.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) *
                                          std::pow(EulerConst, 4.0 * y[2]))
                                       * std::pow(y[1], 2.0)
                                     - 3.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                                     + 4.0 * std::pow(EulerConst, 4.0 * y[2]) * std::pow(Pi, 2.0) * s_omega_2 *
                                       std::pow(s_eta, 3.0)
                             )
                     ) / denominatorDyDx3
                     + (-1.0) * numeratorDyDx3 / std::pow(denominatorDyDx3, 2.0) * (
                (
                        2.0 * 3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                        * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        - 2.0 * std::pow(EulerConst, 2.0 * y[2]) * s_eta * Pi * std::pow(y[0], 2.0) / 3.0
                        - 2.0 * 3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 2.0)
                        + 2.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                ) * y[0] * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
                + (
                          (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                           - s_eta / 4.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                          + (-s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi / 3.0
                             - 1.0 / 12.0) * std::pow(y[0], 2.0)
                          + (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                             + s_eta / 4.0) * std::pow(y[1], 2.0)
                          + std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                          + s_omega_2 * s_eta / 4.0
                  ) * y[0] * 2.0 * s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi
        );

        dfdy[3][3] = (
                             y[1] * (
                                     (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                                      -
                                      9.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                      + 3.0 * s_eta / 16.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                                     +
                                     (-3.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                      + 1.0 / 16.0) * std::pow(y[0], 2.0)
                                     + (-std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                                        - 2.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) *
                                          std::pow(EulerConst, 4.0 * y[2]))
                                       * std::pow(y[1], 2.0) * y[0]
                                     + (2.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                        + 11.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) *
                                          std::pow(EulerConst, 4.0 * y[2])
                                        - 3.0 * s_eta / 16.0) * std::pow(y[1], 2.0)
                                     - 3.0 / 2.0 * std::pow(EulerConst, 2.0 * y[2]) * Pi * s_omega_2 *
                                       std::pow(s_eta, 2.0)
                                     + std::pow(EulerConst, 4.0 * y[2]) * std::pow(Pi, 2.0) * s_omega_2 *
                                       std::pow(s_eta, 3.0)
                                     - 3.0 * s_omega_2 * s_eta / 16.0
                             ) + y[1] * y[3] * (
                                     (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 2.0
                                      -
                                      9.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                      + 3.0 * s_eta / 16.0) * 2.0 * std::pow(y[0], 2.0) * y[3]
                             )
                     ) / denominatorDyDx3
                     + (-1.0) * numeratorDyDx3 / std::pow(denominatorDyDx3, 2.0) * (
                (3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) - s_eta / 4.0)
                * 2.0 * std::pow(y[0], 3.0) * y[3] * (s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi + 1.0 / 4.0)
        );


    }

};

void get_bianchiIX_anisotropy_points();

#endif //ODECALCULATIONS_BIANCHIIX_ANISOTROPY_H
