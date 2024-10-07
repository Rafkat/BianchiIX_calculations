//
// Created by rafkat on 8/26/24.
//

#ifndef ODECALCULATIONS_BIANCHIV_ANISOTROPY_H
#define ODECALCULATIONS_BIANCHIV_ANISOTROPY_H

#include "../constants.h"
#include "../nr/nr3.h"
#include "../nr/odeint.h"
#include <cmath>


struct bianchiV_anisotropic_rhs {
    Doub s_eta;
    Doub s_omega_2;

    bianchiV_anisotropic_rhs(Doub eta, Doub omega_2) : s_eta(eta), s_omega_2(omega_2) {}

    void operator()(const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        Doub common_denominator1 = (
                (
                        (Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0)
                         - s_eta / 12.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta / 9.0
                           - 1.0 / 36.0) * std::pow(y[0], 2.0)
                        + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0)
                           + s_eta / 12.0) * std::pow(y[1], 2.0)
                        - Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0) / 3.0
                        - s_omega_2 * s_eta / 12.0
                )
        );

        Doub common_denominator2 = (
                Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta + 0.25
        );

        dydx[0] = y[1];
        dydx[1] = (
                          (-432.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                           - 72.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0)
                           + 9.0 * s_eta) * std::pow(y[0], 4.0) * std::pow(y[3], 4.0)
                          + (-96.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                             + 36.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                             + 3.0) * std::pow(y[3], 2.0) * std::pow(y[0], 4.0)
                          + (16.0 * s_eta * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                             + 4.0 * Pi * std::pow(EulerConst, 2.0 * y[2])) * std::pow(y[0], 4.0)
                          + (288.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                             + 48.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0)
                             - 6.0 * s_eta) * std::pow(y[0], 2.0) * std::pow(y[1], 2.0) * std::pow(y[3], 2.0)
                          + (-224.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                             - 52.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                             + 1.0) * std::pow(y[0], 2.0) * std::pow(y[1], 2.0)
                          + (120.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                             + 6.0 * s_omega_2 * s_eta) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                          + (64.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * s_omega_2
                             * std::pow(EulerConst, 4.0 * y[2])
                             + 12.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * s_eta
                             - s_omega_2) * std::pow(y[0], 2.0)
                          + (144.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                             + 24.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0)
                             - 3.0 * s_eta) * std::pow(y[1], 4.0)
                          + (-192.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * s_omega_2
                             * std::pow(EulerConst, 4.0 * y[2])
                             - 24.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                             + 6.0 * s_omega_2 * s_eta) * std::pow(y[1], 2.0)
                          + 48.0 * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                            * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0)
                          - 3.0 * std::pow(s_omega_2, 2.0) * s_eta
                  ) / 288.0 / y[0] / (common_denominator1 * common_denominator2);

        dydx[2] = (
                          2.0 * y[1] * (
                                  (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta + 1.0 / 8.0) * std::pow(y[0], 2.0)
                                  + Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                          )
                  ) / 3.0 / y[0] / common_denominator1;

        dydx[3] = -(
                3 * y[1] * y[3] * (
                        (Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0) / 6.0
                         + std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                         - s_eta / 48.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 3.0
                           - 1.0 / 144.0) * std::pow(y[0], 2.0)
                        + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0) / 6.0
                           - std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                           + s_eta / 48.0) * std::pow(y[1], 2.0)
                        - Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0) / 6.0
                        + std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * s_omega_2 * std::pow(EulerConst, 4.0 * y[2]) / 9.0
                        - s_omega_2 * s_eta / 48.0
                )
        ) / common_denominator1 / common_denominator2 / y[0];
    }

    void jacobian(const Doub x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy) {
        Doub common_denominator1 = (
                (
                        (Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0)
                         - s_eta / 12.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                        + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta / 9.0
                           - 1.0 / 36.0) * std::pow(y[0], 2.0)
                        + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0)
                           + s_eta / 12.0) * std::pow(y[1], 2.0)
                        - Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0) / 3.0
                        - s_omega_2 * s_eta / 12.0
                )
        );

        Doub common_denominator2 = (
                Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta + 0.25
        );


        Int n = y.size();
        for (Int i = 0; i < n; i++) dfdx[i] = 0.0;
        dfdy[0][0] = 0.0;
        dfdy[0][1] = 1.0;
        dfdy[0][2] = 0.0;
        dfdy[0][3] = 0.0;

        dfdy[1][0] = (
                             (-15552.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                              - 1296.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                              + 540.0 * std::pow(s_eta, 3.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                              - 27.0 * std::pow(s_eta, 2.0)) * std::pow(y[0], 6.0) * std::pow(y[3], 6.0)
                             + (-1728.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                                + 2304.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                + 36.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0)
                                - 18.0 * s_eta) * std::pow(y[0], 6.0) * std::pow(y[3], 4.0)
                             + (960.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                                + 48.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 60.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                                - 3.0) * std::pow(y[0], 6.0) * std::pow(y[3], 2.0)
                             + (-64.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                                - 32.0 * std::pow(Pi, 2.0) * s_eta * std::pow(EulerConst, 4.0 * y[2])
                                - 4.0 * Pi * std::pow(EulerConst, 2.0 * y[2])) * std::pow(y[0], 6.0)
                             + (36288.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                                + 3024.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 1260.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2])
                                + 63.0 * std::pow(s_eta, 2.0)) * std::pow(y[0], 4.0) * std::pow(y[1], 2.0)
                               * std::pow(y[3], 4.0)
                             + (19584.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                                - 3072.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 168.0 * Pi * std::pow(s_eta, 2.0) * std::pow(EulerConst, 2.0 * y[2])
                                + 24.0 * s_eta) * std::pow(y[0], 4.0) * std::pow(y[1], 2.0) * std::pow(y[3], 2.0)
                             + (-2624.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                                - 720.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 12.0 * Pi * s_eta * std::pow(EulerConst, 2.0 * y[2])
                                + 1.0) * std::pow(y[0], 4.0) * std::pow(y[1], 2.0)
                             + (15552.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * s_omega_2
                                * std::pow(EulerConst, 6.0 * y[2])
                                + 2160.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                                  * std::pow(EulerConst, 4.0 * y[2])
                                + 468.0 * Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2])
                                - 63.0 * s_omega_2 * std::pow(s_eta, 2.0)) * std::pow(y[0], 4.0) * std::pow(y[3], 4.0)
                             + (1152.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * s_omega_2
                                * std::pow(EulerConst, 6.0 * y[2])
                                - 192.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * s_omega_2
                                  * std::pow(EulerConst, 4.0 * y[2])
                                - 216.0 * Pi * s_omega_2 * std::pow(s_eta, 2.0) * std::pow(EulerConst, 2.0 * y[2])
                                - 24.0 * s_omega_2 * s_eta) * std::pow(y[3], 2.0) * std::pow(y[0], 4.0)
                             + (-320.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 3.0) * s_omega_2
                                * std::pow(EulerConst, 6.0 * y[2])
                                - 176.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * s_omega_2
                                  * std::pow(EulerConst, 4.0 * y[2])
                                - 28.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * s_eta
                                - s_omega_2) * std::pow(y[0], 4.0)
                             + (-25920.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                                - 2160.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                + 900.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2])
                                - 45.0 * std::pow(s_eta, 2.0)) * std::pow(y[0], 2.0) * std::pow(y[1], 4.0)
                               * std::pow(y[3], 2.0)
                             + (9792.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                                + 1920.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 156.0 * Pi * std::pow(s_eta, 2.0) * std::pow(EulerConst, 2.0 * y[2])
                                - 6.0 * s_eta) * std::pow(y[0], 2.0) * std::pow(y[1], 4.0)
                             + (17280.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * s_omega_2
                                * std::pow(EulerConst, 6.0 * y[2])
                                - 4896.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                                  * std::pow(EulerConst, 4.0 * y[2])
                                - 792.0 * Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2])
                                + 90.0 * s_omega_2 * std::pow(s_eta, 2.0)) * std::pow(y[0], 2.0)
                               * std::pow(y[1], 2.0) * std::pow(y[3], 2.0)
                             + (-1920.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * s_omega_2
                                * std::pow(EulerConst, 6.0 * y[2])
                                + 192.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * s_omega_2
                                  * std::pow(EulerConst, 4.0 * y[2])
                                + 216.0 * Pi * s_omega_2 * std::pow(s_eta, 2.0) * std::pow(EulerConst, 2.0 * y[2])
                                + 12.0 * s_omega_2 * s_eta) * std::pow(y[0], 2.0) * std::pow(y[1], 2.0)
                             + (-5184.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(s_omega_2, 2.0)
                                * std::pow(EulerConst, 6.0 * y[2])
                                - 1008.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                                  * std::pow(EulerConst, 4.0 * y[2])
                                - 108.0 * Pi * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0)
                                  * std::pow(EulerConst, 2.0 * y[2])
                                - 45.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0))
                               * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                             + (-192.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(s_omega_2, 2.0)
                                * std::pow(EulerConst, 6.0 * y[2])
                                - 192.0 * std::pow(Pi, 2.0) * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0)
                                  * std::pow(EulerConst, 4.0 * y[2])
                                - 60.0 * Pi * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0)
                                  * std::pow(EulerConst, 2.0 * y[2])
                                - 6.0 * std::pow(s_omega_2, 2.0) * s_eta) * std::pow(y[0], 2.0)
                             + (5184.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                                + 432.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 180.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2])
                                + 9.0 * std::pow(s_eta, 2.0)) * std::pow(y[1], 6.0)
                             + (-5184.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * s_omega_2
                                * std::pow(EulerConst, 6.0 * y[2])
                                + 432.0 * std::pow(s_eta, 4.0) * s_omega_2 * std::pow(Pi, 2.0)
                                  * std::pow(EulerConst, 4.0 * y[2])
                                + 324.0 * Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2])
                                - 27.0 * s_omega_2 * std::pow(s_eta, 2.0)) * std::pow(y[1], 4.0)
                             + (-576.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(s_omega_2, 2.0)
                                * std::pow(EulerConst, 6.0 * y[2])
                                - 1008.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                                  * std::pow(EulerConst, 4.0 * y[2])
                                - 108.0 * Pi * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0)
                                  * std::pow(EulerConst, 2.0 * y[2])
                                + 27.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0)) * std::pow(y[1], 2.0)
                             + 576.0 * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                               * std::pow(s_omega_2, 3.0) * std::pow(s_eta, 5.0)
                             + 144.0 * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                               * std::pow(s_omega_2, 3.0) * std::pow(s_eta, 4.0)
                             - 36.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_omega_2, 3.0)
                               * std::pow(s_eta, 3.0)
                             - 9.0 * std::pow(s_omega_2, 3.0) * std::pow(s_eta, 2.0)
                     ) / 10368.0 / std::pow(y[0], 2.0)
                     / std::pow(common_denominator1, 2.0) / common_denominator2;

        dfdy[1][1] = -y[1] * (
                (-5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 144.0
                 + std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 12.0
                 + std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                 + std::pow(s_eta, 2.0) / 576.0) * std::pow(y[0], 4.0) * std::pow(y[3], 4.0)
                + (-7.0 * Pi * std::pow(s_eta, 2.0) * std::pow(EulerConst, 2.0 * y[2]) / 216.0
                   + std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 54.0
                   + 22.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                   + s_eta / 864.0) * std::pow(y[3], 2.0) * std::pow(y[0], 4.0)
                + (-Pi * s_eta * std::pow(EulerConst, 2.0 * y[2]) / 144.0
                   - 11.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 108.0
                   - 23.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 81.0
                   + 1.0 / 5184.0) * std::pow(y[0], 4.0)
                + (5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   - std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 6.0
                   - 2.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                   - std::pow(s_eta, 2.0) / 288.0) * std::pow(y[0], 2.0) * std::pow(y[1], 2.0) * std::pow(y[3], 2.0)
                + (Pi * std::pow(s_eta, 2.0) * std::pow(EulerConst, 2.0 * y[2]) / 216.0
                   + 5.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 54.0
                   + 2.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                   - s_eta / 864.0) * std::pow(y[0], 2.0) * std::pow(y[1], 2.0)
                + (-Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   - s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 2.0
                   + 2.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) * s_omega_2
                   + s_omega_2 * std::pow(s_eta, 2.0) / 288.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0) / 72.0
                   - 19.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * s_omega_2
                     * std::pow(EulerConst, 4.0 * y[2]) / 54.0
                   - 10.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * s_omega_2
                     * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                   + s_omega_2 * s_eta / 864.0) * std::pow(y[0], 2.0)
                + (-5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 144.0
                   + std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 12.0
                   + std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                   + std::pow(s_eta, 2.0) / 576.0) * std::pow(y[1], 4.0)
                + (Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   + 5.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 18.0
                   + 2.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * s_omega_2
                     * std::pow(EulerConst, 6.0 * y[2]) / 3.0
                   - s_omega_2 * std::pow(s_eta, 2.0) / 288.0) * std::pow(y[1], 2.0)
                + Pi * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 48.0
                - 5.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                  * std::pow(EulerConst, 4.0 * y[2]) / 36.0
                - 7.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(s_omega_2, 2.0)
                  * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                + std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0) / 576.0
        ) / std::pow(common_denominator1, 2.0) / common_denominator2 / y[0];

        dfdy[1][2] = -Pi * std::pow(EulerConst, 2.0 * y[2]) * (
                (-Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 6.0
                 + std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                 + std::pow(s_eta, 2.0) / 144.0) * std::pow(y[0], 6.0) * std::pow(y[3], 4.0)
                + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0) / 27.0
                   - 2.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 9.0
                   + s_eta / 216.0) * std::pow(y[0], 6.0) * std::pow(y[3], 2.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta / 162.0
                   + std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 81.0
                   + 1.0 / 1296.0) * std::pow(y[0], 6.0)
                + (-Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 9.0
                   - 14.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 9.0
                   - std::pow(s_eta, 2.0) / 24.0) * std::pow(y[0], 4.0) * std::pow(y[1], 2.0) * std::pow(y[3], 2.0)
                + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0) / 9.0
                   - 2.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 9.0
                   - s_eta / 72.0) * std::pow(y[0], 4.0) * std::pow(y[1], 2.0)
                + (-std::pow(s_eta, 4.0) * s_omega_2 * Pi * std::pow(EulerConst, 2.0 * y[2]) / 3.0
                   + 2.0 * std::pow(s_eta, 5.0) * s_omega_2 * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   + std::pow(s_eta, 3.0) * s_omega_2 / 72.0) * std::pow(y[0], 4.0) * std::pow(y[3], 4.0)
                + (-5.0 * Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 27.0
                   - 10.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 9.0
                   + 5.0 * s_omega_2 * std::pow(s_eta, 2.0) / 216.0) * std::pow(y[3], 2.0) * std::pow(y[0], 4.0)
                + (4.0 * Pi * s_omega_2 * std::pow(s_eta, 2.0) * std::pow(EulerConst, 2.0 * y[2]) / 81.0
                   + 8.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                     * s_omega_2 / 81.0
                   + s_omega_2 * s_eta / 162.0) * std::pow(y[0], 4.0)
                + (5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 18.0
                   + 5.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 9.0
                   + 5.0 * std::pow(s_eta, 2.0) / 144.0) * std::pow(y[0], 2.0) * std::pow(y[1], 4.0)
                + (2.0 * std::pow(s_eta, 4.0) * s_omega_2 * Pi * std::pow(EulerConst, 2.0 * y[2]) / 9.0
                   - 20.0 * std::pow(s_eta, 5.0) * s_omega_2 * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 9.0
                   - std::pow(s_eta, 3.0) * s_omega_2 / 36.0) * std::pow(y[0], 2.0) * std::pow(y[1], 2.0)
                  * std::pow(y[3], 2.0)
                + 11.0 / 27.0 * (-Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2])
                                 - 2.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                                   * std::pow(EulerConst, 4.0 * y[2])
                                 - s_omega_2 * std::pow(s_eta, 2.0) / 8.0) * std::pow(y[0], 2.0) * std::pow(y[1], 2.0)
                + (-2.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 4.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 9.0
                   - 4.0 * std::pow(s_eta, 5.0) * std::pow(s_omega_2, 2.0) * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 3.0
                   + std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0) / 36.0) * std::pow(y[0], 2.0) * std::pow(y[3], 2.0)
                + (7.0 * Pi * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 54.0
                   + 7.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 27.0
                   + 7.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0) / 432.0) * std::pow(y[0], 2.0)
                + (std::pow(s_eta, 4.0) * s_omega_2 * Pi * std::pow(EulerConst, 2.0 * y[2]) / 9.0
                   + 2.0 * std::pow(s_eta, 5.0) * s_omega_2 * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 9.0
                   + std::pow(s_eta, 3.0) * s_omega_2 / 72.0) * std::pow(y[1], 4.0)
                + (-2.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 4.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 9.0
                   - 4.0 * std::pow(s_eta, 5.0) * std::pow(s_omega_2, 2.0) * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 9.0
                   - std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0) / 36.0) * std::pow(y[1], 2.0)
                + Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_omega_2, 3.0) * std::pow(s_eta, 4.0) / 9.0
                + 2.0 * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) * std::pow(s_omega_2, 3.0)
                  * std::pow(s_eta, 5.0) / 9.0
                + std::pow(s_omega_2, 3.0) * std::pow(s_eta, 3.0) / 72.0
        ) / std::pow(common_denominator1 * common_denominator2, 2.0) / 4.0 / y[0];

        dfdy[1][3] = -3.0 * y[0] * y[3] * (
                (-5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 144.0
                 + std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 12.0
                 + std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                 + std::pow(s_eta, 2.0) / 576.0) * std::pow(y[0], 4.0) * std::pow(y[3], 4.0)
                + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0) / 216.0
                   - 5.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 54.0
                   - 2.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                   + s_eta / 864.0) * std::pow(y[3], 2.0) * std::pow(y[0], 4.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta / 432.0
                   + std::pow(s_eta * Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 108.0
                   + std::pow(s_eta * Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 81.0
                   + 1.0 / 5184.0) * std::pow(y[0], 4.0)
                + (5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   - std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 6.0
                   - 2.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                   - std::pow(s_eta, 2.0) / 288.0) * std::pow(y[0] * y[1] * y[3], 2.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0) / 72.0
                   + std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 18.0
                   - 2.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 3.0
                   - s_eta / 864.0) * std::pow(y[0] * y[1], 2.0)
                + (-Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   - 5.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 18.0
                   - 2.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * s_omega_2
                     * std::pow(EulerConst, 6.0 * y[2]) / 3.0
                   + s_omega_2 * std::pow(s_eta, 2.0) / 288.0) * std::pow(y[0] * y[3], 2.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0) / 72.0
                   + std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * s_omega_2 * std::pow(EulerConst, 4.0 * y[2]) / 18.0
                   + 2.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * s_omega_2
                     * std::pow(EulerConst, 6.0 * y[2]) / 27.0
                   + s_omega_2 * s_eta / 864.0) * std::pow(y[0], 2.0)
                + (-5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 144.0
                   + std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 12.0
                   + std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                   + std::pow(s_eta, 2.0) / 576.0) * std::pow(y[1], 4.0)
                + (Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   + 19.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 54.0
                   - 2.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * s_omega_2
                     * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                   - s_omega_2 * std::pow(s_eta, 2.0) / 288.0) * std::pow(y[1], 2.0)
                + Pi * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 48.0
                + std::pow(s_omega_2, 2.0) * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                  * std::pow(EulerConst, 4.0 * y[2]) / 12.0
                + std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(s_omega_2, 2.0)
                  * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                + std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0) / 576.0
        ) / std::pow(common_denominator1, 2.0) / common_denominator2;

        dfdy[2][0] = -2.0 * y[1] * (
                (Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0) / 24.0
                 + std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                 - s_eta / 96.0) * std::pow(y[3], 2.0) * std::pow(y[0], 4.0)
                + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta / 24.0
                   - std::pow(s_eta * Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 9.0
                   - 1.0 / 288.0) * std::pow(y[0], 4.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0) / 24.0
                   + std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   - s_eta / 96.0) * std::pow(y[1] * y[0], 2.0)
                + (-Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 4.0
                   + 3.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]))
                  * std::pow(y[0] * y[3], 2.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0) / 24.0
                   + s_omega_2 * s_eta * 96.0) * std::pow(y[0], 2.0)
                + (Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 12.0
                   - s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]))
                  * std::pow(y[1], 2.0)
                - Pi * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 12.0
                - std::pow(s_omega_2, 2.0) * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                  * std::pow(EulerConst, 4.0 * y[2]) / 3.0
        ) / std::pow(common_denominator1, 2.0) / 3.0 / std::pow(y[0], 2.0);

        dfdy[2][1] = (
                             2.0 * (
                                     (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta + 1.0 / 8.0) * std::pow(y[0], 2.0)
                                     + Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                             ) * (
                                     (Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0)
                                      - s_eta / 12.0) * std::pow(y[0] * y[3], 2.0)
                                     + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta / 9.0
                                        - 1.0 / 36.0) * std::pow(y[0], 2.0)
                                     + (Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0)
                                        - s_eta / 12.0) * std::pow(y[1], 2.0)
                                     - Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0) / 3.0
                                     - s_omega_2 * s_eta / 12.0
                             )
                     ) / 3.0 / y[0]
                     / std::pow(common_denominator1, 2.0);

        dfdy[2][2] = -5.0 * y[1] * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]) * (
                2.0 / 5.0 * s_omega_2 * std::pow(s_eta * y[0] * y[3], 2.0)
                - 2.0 / 5.0 * s_omega_2 * std::pow(s_eta * y[1], 2.0)
                + 2.0 / 5.0 * std::pow(s_omega_2 * s_eta, 2.0)
                + s_eta * std::pow(y[0], 4.0) * std::pow(y[3], 2.0)
                - s_eta * std::pow(y[0] * y[1], 2.0)
                + 1.0 / 3.0 * s_omega_2 * s_eta * std::pow(y[0], 2.0)
                + 1.0 / 15.0 * std::pow(y[0], 4.0)
        ) / 18.0 / y[0] / std::pow(common_denominator1, 2.0);

        dfdy[2][3] = -4.0 * y[1] * s_eta * y[0] * y[3] * (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                                                          - 1.0 / 12.0) * (
                             (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta + 1.0 / 8.0) * std::pow(y[0], 2.0)
                             + Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                     ) / 3.0 / std::pow(common_denominator1, 2.0);

        dfdy[3][0] = 3.0 * y[1] * y[3] * (
                (-5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 144.0
                 + std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 12.0
                 + std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                 + std::pow(s_eta, 2.0) / 576.0) * std::pow(y[0] * y[3], 4.0)
                + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0) / 108.0
                   - 2.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 27.0
                   + 2.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                   + s_eta / 864.0) * std::pow(y[3], 2.0) * std::pow(y[0], 4.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta / 1296.0
                   - std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 108.0
                   - std::pow(s_eta, 3.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 27.0
                   + 1.0 / 5184.0) * std::pow(y[0], 4.0)
                + (5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   - std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 6.0
                   - 2.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                   - std::pow(s_eta, 2.0) / 288.0) * std::pow(y[1] * y[0] * y[3], 2.0)
                + (std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 9.0
                   + 2.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 3.0
                   - s_eta / 864.0) * std::pow(y[0] * y[1], 2.0)
                + (-Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   - 7.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 18.0
                   + 2.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * s_omega_2
                     * std::pow(EulerConst, 6.0 * y[2]) / 3.0
                   + s_omega_2 * std::pow(s_eta, 2.0) / 288.0) * std::pow(y[0] * y[3], 2.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0) / 54.0
                   + 2.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * s_omega_2
                     * std::pow(EulerConst, 4.0 * y[2]) / 27.0
                   + 2.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0)
                     * s_omega_2 * std::pow(EulerConst, 6.0 * y[2]) / 27.0
                   + s_omega_2 * s_eta / 864.0) * std::pow(y[0], 2.0)
                + (-5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 144.0
                   + std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 12.0
                   + std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                   + std::pow(s_eta, 2.0) / 576.0) * std::pow(y[1], 4.0)
                + (Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   + 17.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 54.0
                   + 2.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * s_omega_2
                     * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                   - s_omega_2 * std::pow(s_eta, 2.0) / 288.0) * std::pow(y[1], 2.0)
                + Pi * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 48.0
                + 5.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                  * std::pow(EulerConst, 4.0 * y[2]) / 108.0
                - std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(s_omega_2, 2.0)
                  * std::pow(EulerConst, 6.0 * y[2]) / 27.0
                + std::pow(s_omega_2 * s_eta, 2.0) / 576.0
        ) / y[0] / y[0] / std::pow(common_denominator1, 2.0) / common_denominator2;

        dfdy[3][1] = -3.0 * y[3] * (
                (-5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 144.0
                 + std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 12.0
                 + std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                 + std::pow(s_eta, 2.0) / 576.0) * std::pow(y[0] * y[3], 4.0)
                + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0) / 108.0
                   - 2.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 27.0
                   + 2.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                   + s_eta / 864.0) * std::pow(y[3], 2.0) * std::pow(y[0], 4.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta / 1296.0
                   - std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 108.0
                   - std::pow(s_eta, 3.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 27.0
                   + 1.0 / 5184.0) * std::pow(y[0], 4.0)
                + (5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   - std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 6.0
                   - 2.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                   - std::pow(s_eta, 2.0) / 288.0) * std::pow(y[1] * y[0] * y[3], 2.0)
                + (std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 9.0
                   + 2.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 3.0
                   - s_eta / 864.0) * std::pow(y[0] * y[1], 2.0)
                + (-Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   - 17.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 54.0
                   - 2.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * s_omega_2
                     * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                   + s_omega_2 * std::pow(s_eta, 2.0) / 288.0) * std::pow(y[0] * y[3], 2.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0) / 108.0
                   - std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * s_omega_2
                     * std::pow(EulerConst, 4.0 * y[2]) / 81.0
                   - 10.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0)
                     * s_omega_2 * std::pow(EulerConst, 6.0 * y[2]) / 81.0
                   + s_omega_2 * s_eta / 864.0) * std::pow(y[0], 2.0)
                + (-5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 144.0
                   + std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 12.0
                   + std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                   + std::pow(s_eta, 2.0) / 576.0) * std::pow(y[1], 4.0)
                + (Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   + 13.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 54.0
                   + 10.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * s_omega_2
                     * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                   - s_omega_2 * std::pow(s_eta, 2.0) / 288.0) * std::pow(y[1], 2.0)
                + Pi * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 48.0
                + 5.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                  * std::pow(EulerConst, 4.0 * y[2]) / 108.0
                - std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(s_omega_2, 2.0)
                  * std::pow(EulerConst, 6.0 * y[2]) / 27.0
                + std::pow(s_omega_2 * s_eta, 2.0) / 576.0
        ) / y[0] / std::pow(common_denominator1, 2.0) / common_denominator2;

        dfdy[3][2] = -y[1] * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]) * y[3] * (
                (-Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0)
                 + std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                 - s_eta / 16.0) * std::pow(y[3], 2.0) * std::pow(y[0], 4.0)
                + (-Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta / 3.0
                   - std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   - 1.0 / 48.0) * std::pow(y[0], 4.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0)
                   - std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   + s_eta / 16.0) * std::pow(y[0] * y[1], 2.0)
                + (4.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   - Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]))
                  * std::pow(y[0] * y[3], 2.0)
                + (-4.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0) / 3.0
                   - 13.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * s_omega_2
                     * std::pow(EulerConst, 4.0 * y[2]) / 3.0
                   - s_omega_2 * s_eta / 16.0) * std::pow(y[0], 2.0)
                + (-4.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                   + Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2])) * std::pow(y[1], 2.0)
                - 4.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                  * std::pow(EulerConst, 4.0 * y[2])
                - Pi * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2])
        ) / 9.0 / y[0] / std::pow(common_denominator1 * common_denominator2, 2.0);

        dfdy[3][3] = -3.0 * y[1] * (
                (-5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 144.0
                 + std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 12.0
                 + std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                 + std::pow(s_eta, 2.0) / 576.0) * std::pow(y[0] * y[3], 4.0)
                + (-std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 9.0
                   - 2.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 3.0
                   + s_eta / 864.0) * std::pow(y[3], 2.0) * std::pow(y[0], 4.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_eta / 1296.0
                   - std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 108.0
                   - std::pow(s_eta, 3.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 27.0
                   + 1.0 / 5184.0) * std::pow(y[0], 4.0)
                + (5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   - std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 6.0
                   - 2.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                   - std::pow(s_eta, 2.0) / 288.0) * std::pow(y[1] * y[0] * y[3], 2.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_eta, 2.0) / 108.0
                   + 2.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 27.0
                   - 2.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                   - s_eta / 864.0) * std::pow(y[0] * y[1], 2.0)
                + (-Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   - 13.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 54.0
                   - 10.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * s_omega_2
                     * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                   + s_omega_2 * std::pow(s_eta, 2.0) / 288.0) * std::pow(y[0] * y[3], 2.0)
                + (Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0) / 108.0
                   - std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * s_omega_2
                     * std::pow(EulerConst, 4.0 * y[2]) / 81.0
                   - 10.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 3.0)
                     * s_omega_2 * std::pow(EulerConst, 6.0 * y[2]) / 81.0
                   + s_omega_2 * s_eta / 864.0) * std::pow(y[0], 2.0)
                + (-5.0 * Pi * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 144.0
                   + std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 12.0
                   + std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(EulerConst, 6.0 * y[2])
                   + std::pow(s_eta, 2.0) / 576.0) * std::pow(y[1], 4.0)
                + (Pi * s_omega_2 * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 72.0
                   + 17.0 * s_omega_2 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                     * std::pow(EulerConst, 4.0 * y[2]) / 54.0
                   + 2.0 * std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * s_omega_2
                     * std::pow(EulerConst, 6.0 * y[2]) / 9.0
                   - s_omega_2 * std::pow(s_eta, 2.0) / 288.0) * std::pow(y[1], 2.0)
                + Pi * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0) * std::pow(EulerConst, 2.0 * y[2]) / 48.0
                + 5.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                  * std::pow(EulerConst, 4.0 * y[2]) / 108.0
                - std::pow(s_eta, 5.0) * std::pow(Pi, 3.0) * std::pow(s_omega_2, 2.0)
                  * std::pow(EulerConst, 6.0 * y[2]) / 27.0
                + std::pow(s_omega_2 * s_eta, 2.0) / 576.0
        ) / y[0] / std::pow(common_denominator1, 2.0) / common_denominator2;


    }

};

void get_bianchiV_anisotropy_points(int);

#endif //ODECALCULATIONS_BIANCHIV_ANISOTROPY_H
