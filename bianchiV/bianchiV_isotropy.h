//
// Created by rafkat on 8/19/24.
//

#ifndef ODECALCULATIONS_BIANCHIV_ISOTROPY_H
#define ODECALCULATIONS_BIANCHIV_ISOTROPY_H

#include "../constants.h"
#include "../nr/nr3.h"
#include "../nr/odeint.h"
#include <cmath>


struct bianchiV_rhs {
    Doub s_eta;
    Doub s_omega_2;

    bianchiV_rhs(Doub eta, Doub omega_2) : s_eta(eta), s_omega_2(omega_2) {}

    void operator()(const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        Doub common_denominator = (
                y[0] * (
                        (4.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]) + 1.0) * std::pow(y[0], 2.0)
                        + (36.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) - 3.0 * s_eta)
                          * std::pow(y[1], 2.0)
                        + 12.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                        + 3.0 * s_omega_2 * s_eta
                )
        );

        dydx[0] = y[1];
        dydx[1] = -(
                4.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[0], 4.0)
                + (-56.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]) + 1.0)
                  * std::pow(y[1], 2.0) * std::pow(y[0], 2.0)
                + (16.0 * s_omega_2 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]) - s_omega_2) * std::pow(y[0], 2.0)
                + (36.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) - 3.0 * s_eta)
                  * std::pow(y[1], 4.0)
                + (-48.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                   + 6.0 * s_omega_2 * s_eta) * std::pow(y[1], 2.0)
                + 12.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0)
                - 3.0 * std::pow(s_omega_2, 2.0) * s_eta
        ) / 2.0 / common_denominator;

        dydx[2] = -(
                3.0 * y[1] * (
                        (8.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]) + 1.0) * std::pow(y[0], 2.0)
                        + 8.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                )
        ) / common_denominator;
    }

    void jacobian(const Doub x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy) {
        Doub common_denominator = (
                y[0] * (
                        (4.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]) + 1.0) * std::pow(y[0], 2.0)
                        + (36.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) - 3.0 * s_eta)
                          * std::pow(y[1], 2.0)
                        + 12.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                        + 3.0 * s_omega_2 * s_eta
                )
        );

        Int n = y.size();
        for (Int i = 0; i < n; i++) dfdx[i] = 0.0;
        dfdy[0][0] = 0.0;
        dfdy[0][1] = 1.0;
        dfdy[0][2] = 0.0;

        dfdy[1][0] = (
                             (-16.0 * s_eta * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                              - 4.0 * Pi * std::pow(EulerConst, 2.0 * y[2])) * std::pow(y[0], 6.0)
                             + (656.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 16.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2])
                                + 1.0) * std::pow(y[0], 4.0) * std::pow(y[1], 2.0)
                             + (80.0 * std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * s_omega_2
                                * std::pow(EulerConst, 4.0 * y[2])
                                - 24.0 * s_omega_2 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2])
                                - s_omega_2) * std::pow(y[0], 4.0)
                             + (2448.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 132.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                - 6.0 * s_eta) * std::pow(y[0], 2.0) * std::pow(y[1], 4.0)
                             + (-480.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * s_omega_2
                                * std::pow(EulerConst, 4.0 * y[2])
                                + 168.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                                + 12.0 * s_omega_2 * s_eta) * std::pow(y[1], 2.0) * std::pow(y[0], 2.0)
                             + (-48.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * std::pow(s_omega_2, 2.0)
                                * std::pow(EulerConst, 4.0 * y[2])
                                - 36.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_omega_2, 2.0)
                                  * std::pow(s_eta, 2.0)
                                - 6.0 * std::pow(s_omega_2, 2.0) * s_eta) * std::pow(y[0], 2.0)
                             + (1296.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                - 216.0 * std::pow(s_eta, 3.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                                + 9.0 * std::pow(s_eta, 2.0)) * std::pow(y[1], 6.0)
                             + (-1296.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * s_omega_2
                                * std::pow(EulerConst, 4.0 * y[2])
                                + 432.0 * std::pow(s_eta, 3.0) * s_omega_2 * Pi * std::pow(EulerConst, 2.0 * y[2])
                                - 27.0 * s_omega_2 * std::pow(s_eta, 2.0)) * std::pow(y[1], 4.0)
                             + (-144.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(s_omega_2, 2.0)
                                * std::pow(EulerConst, 4.0 * y[2])
                                - 216.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0) * Pi
                                  * std::pow(EulerConst, 2.0 * y[2])
                                + 27.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0)) * std::pow(y[1], 2.0)
                             + 144.0 * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2]) * std::pow(s_omega_2, 3.0)
                               * std::pow(s_eta, 4.0)
                             - 9.0 * std::pow(s_omega_2, 3.0) * std::pow(s_eta, 2.0)
                     ) / 16.0 / 32.0 / std::pow(common_denominator, 2.0);

        dfdy[1][1] = (
                             23.0 * y[1] * (
                                     (5.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]) / 46.0
                                      + std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                      - 1.0 / 368.0) * std::pow(y[0], 4.0)
                                     + (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 23.0
                                        - 18.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0)
                                          * std::pow(EulerConst, 4.0 * y[2]) / 23.0
                                        + 3.0 * s_eta / 184.0) * std::pow(y[1], 2.0) * std::pow(y[0], 2.0)
                                     + (6.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2
                                        * std::pow(s_eta, 2.0) / 23.0
                                        + 90.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0) * s_omega_2
                                          * std::pow(EulerConst, 4.0 * y[2]) / 23.0
                                        - 3.0 * s_omega_2 * s_eta / 184.0) * std::pow(y[0], 2.0)
                                     + (27.0 * std::pow(s_eta, 3.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 46.0
                                        - 81.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                                          * std::pow(EulerConst, 4.0 * y[2]) / 23.0
                                        - 9.0 * std::pow(s_eta, 2.0) / 368.0) * std::pow(y[1], 4.0)
                                     + (-9.0 * s_omega_2 * std::pow(s_eta, 3.0) * Pi
                                        * std::pow(EulerConst, 2.0 * y[2]) / 23.0
                                        - 54.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * s_omega_2
                                          * std::pow(EulerConst, 4.0 * y[2]) / 23.0
                                        + 9.0 * s_omega_2 * std::pow(s_eta, 2.0) / 184.0) * std::pow(y[1], 2.0)
                                     - 9.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0) * Pi
                                       * std::pow(EulerConst, 2.0 * y[2]) / 46.0
                                     + 63.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0)
                                       * std::pow(s_omega_2, 2.0) * std::pow(EulerConst, 4.0 * y[2]) / 23.0
                                     - 9.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0) / 368.0
                             )
                     ) / 16.0 / (std::pow(common_denominator, 2.0) / y[0]);

        dfdy[1][2] = -(
                (-3.0 * std::pow(y[1], 2.0) * s_eta + 3.0 * s_omega_2 * s_eta + std::pow(y[0], 2.0))
                * std::pow(EulerConst, 2.0 * y[2]) * Pi
                * (-6.0 * s_omega_2 * std::pow(s_eta, 2.0) * std::pow(y[1], 2.0)
                   - 15.0 * s_eta * std::pow(y[0], 2.0) * std::pow(y[1], 2.0)
                   + 6.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0)
                   + 5.0 * s_omega_2 * s_eta * std::pow(y[0], 2.0)
                   + std::pow(y[0], 4.0))
        ) / 4.0 / 16.0 / (std::pow(common_denominator, 2.0) / y[0]);

        dfdy[2][0] = (
                             6.0 * y[1] * (
                                     (3.0 * s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]) / 8.0
                                      + std::pow(s_eta, 2.0) * std::pow(Pi, 2.0) * std::pow(EulerConst, 4.0 * y[2])
                                      + 1.0 / 32.0) * std::pow(y[0], 4.0)
                                     + (-3.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2]) / 8.0
                                        - 9.0 * std::pow(s_eta, 3.0) * std::pow(Pi, 2.0)
                                          * std::pow(EulerConst, 4.0 * y[2])
                                        + 3.0 * s_eta / 32.0) * std::pow(y[1], 2.0) * std::pow(y[0], 2.0)
                                     + (-3.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2
                                        * std::pow(s_eta, 2.0) / 8.0
                                        - 3.0 * s_omega_2 * s_eta / 32.0) * std::pow(y[0], 2.0)
                                     + (-3.0 * s_omega_2 * std::pow(s_eta, 3.0) * Pi
                                        * std::pow(EulerConst, 2.0 * y[2]) / 4.0
                                        + 9.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * s_omega_2
                                          * std::pow(EulerConst, 4.0 * y[2])) * std::pow(y[1], 2.0)
                                     + 3.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 3.0) * Pi
                                       * std::pow(EulerConst, 2.0 * y[2]) / 4.0
                                     + 3.0 * std::pow(s_eta, 4.0) * std::pow(Pi, 2.0) * std::pow(s_omega_2, 2.0)
                                       * std::pow(EulerConst, 4.0 * y[2])
                             )
                     ) / 16.0 / std::pow(common_denominator, 2.0);

        dfdy[2][1] = -6.0 * (
                ((s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]) + 1.0 / 8.0) * std::pow(y[0], 2.0)
                 + Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0))
                * (
                        (s_eta * Pi * std::pow(EulerConst, 2.0 * y[2]) + 1.0 / 4.0) * std::pow(y[0], 2.0)
                        + (-9.0 * std::pow(s_eta, 2.0) * Pi * std::pow(EulerConst, 2.0 * y[2])
                           + 3.0 * s_eta / 4.0) * std::pow(y[1], 2.0)
                        + 3.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                        + 3.0 * s_omega_2 * s_eta / 4.0
                )
        ) / 16.0 / (std::pow(common_denominator, 2.0) / y[0]);

        dfdy[2][2] = -3.0 * y[1] * s_eta * std::pow(EulerConst, 2.0 * y[2]) * Pi * (
                -6.0 * s_omega_2 * std::pow(s_eta, 2.0) * std::pow(y[1], 2.0)
                - 15.0 * s_eta * std::pow(y[0], 2.0) * std::pow(y[1], 2.0)
                + 6.0 * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0)
                + 5.0 * s_omega_2 * s_eta * std::pow(y[0], 2.0)
                + std::pow(y[0], 4.0)
        ) / 16.0 / 2.0 / (std::pow(common_denominator, 2.0) / y[0]);
    }

};

void get_bianchiV_isotropy_points(int);

#endif //ODECALCULATIONS_BIANCHIV_ISOTROPY_H
