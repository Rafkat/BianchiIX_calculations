//
// Created by rafkat on 8/19/24.
//

#ifndef ODECALCULATIONS_BIANCHIIX_ISOTROPY_H
#define ODECALCULATIONS_BIANCHIIX_ISOTROPY_H

#include "../constants.h"
#include "../nr/nr3.h"
#include "../nr/odeint.h"
#include <cmath>


struct bianchiIX_rhs {
    Doub s_eta;
    Doub s_omega_2;

    bianchiIX_rhs(Doub eta, Doub omega_2) : s_eta(eta), s_omega_2(omega_2) {}

    void operator()(const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        Doub numeratorDyDx1 = -(
                36.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 4.0) * std::pow(s_eta, 2.0)
                - 56.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 2.0) * s_eta
                + 48.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 2.0) * s_omega_2 * std::pow(s_eta, 2.0)
                + 4.0 * Pi * std::pow(y[0], 4.0) * std::pow(EulerConst, 2.0 * y[2])
                - 16.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * s_eta
                + 12.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0)
                - 3.0 * std::pow(y[1], 4.0) * s_eta
                + std::pow(y[0], 2.0) * std::pow(y[1], 2.0)
                - 6.0 * std::pow(y[1], 2.0) * s_omega_2 * s_eta
                + std::pow(y[0], 2.0) * s_omega_2
                - 3.0 * std::pow(s_omega_2, 2.0) * s_eta
        );
        Doub denominatorDyDx1 = 2.0 * y[0] * (
                36.0 * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 2.0) * Pi * std::pow(s_eta, 2.0)
                + 4.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                - 12.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                - 3.0 * std::pow(y[1], 2.0) * s_eta
                + std::pow(y[0], 2.0)
                - 3.0 * s_omega_2 * s_eta
        );

        Doub numeratorDyDx2 = -3.0 * y[1] * (
                8.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                - 8.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                + std::pow(y[0], 2.0)
        );

        Doub denominatorDyDx2 = y[0] * (
                36.0 * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 2.0) * Pi * std::pow(s_eta, 2.0)
                + 4.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                - 12.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                - 3.0 * y[1] * s_eta
                + std::pow(y[0], 2.0)
                - 3.0 * s_omega_2 * s_eta
        );


        dydx[0] = y[1];
        dydx[1] = numeratorDyDx1 / (denominatorDyDx1);
        dydx[2] = numeratorDyDx2 / (denominatorDyDx2);
    }

    void jacobian(const Doub x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy) {
        Doub numeratorDyDx1 = -(
                36.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 4.0) * std::pow(s_eta, 2.0)
                - 56.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 2.0) * s_eta
                + 48.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 2.0) * s_omega_2 * std::pow(s_eta, 2.0)
                + 4.0 * Pi * std::pow(y[0], 4.0) * std::pow(EulerConst, 2.0 * y[2])
                - 16.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * s_eta
                + 12.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0)
                - 3.0 * std::pow(y[1], 4.0) * s_eta
                + std::pow(y[0], 2.0) * std::pow(y[1], 2.0)
                - 6.0 * std::pow(y[1], 2.0) * s_omega_2 * s_eta
                + std::pow(y[0], 2.0) * s_omega_2
                - 3.0 * std::pow(s_omega_2, 2.0) * s_eta
        );
        Doub denominatorDyDx1 = 2.0 * y[0] * (
                36.0 * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 2.0) * Pi * std::pow(s_eta, 2.0)
                + 4.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                - 12.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                - 3.0 * std::pow(y[1], 2.0) * s_eta
                + std::pow(y[0], 2.0)
                - 3.0 * s_omega_2 * s_eta
        );

        Doub numeratorDyDx2 = -3.0 * y[1] * (
                8.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                - 8.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                + std::pow(y[0], 2.0)
        );

        Doub denominatorDyDx2 = y[0] * (
                36.0 * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 2.0) * Pi * std::pow(s_eta, 2.0)
                + 4.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                - 12.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                - 3.0 * y[1] * s_eta
                + std::pow(y[0], 2.0)
                - 3.0 * s_omega_2 * s_eta
        );

        Int n = y.size();
        for (Int i = 0; i < n; i++) dfdx[i] = 0.0;
        dfdy[0][0] = 0.0;
        dfdy[0][1] = 1.0;
        dfdy[0][2] = 0.0;

        dfdy[1][0] = -(
                -112.0 * Pi * y[0] * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 2) * s_eta
                + 16.0 * Pi * std::pow(y[0], 3.0) * std::pow(EulerConst, 2.0 * y[2])
                - 32.0 * Pi * y[0] * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * s_eta
                + 2.0 * y[0] * std::pow(y[1], 2.0)
                + 2.0 * y[0] * s_omega_2
        ) / (denominatorDyDx1)
                     + (-1.0) * numeratorDyDx1 / std::pow(denominatorDyDx1, 2.0) * (
                2.0 * (
                        36.0 * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 2.0) * Pi * std::pow(s_eta, 2.0)
                        + 4.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                        - 12.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                        - 3.0 * std::pow(y[1], 2.0) * s_eta
                        + y[0] * y[0]
                        - 3.0 * s_omega_2 * s_eta
                ) + 2.0 * y[0] * (
                        8.0 * Pi * y[0] * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                        + 2.0 * y[0]
                )
        );

        dfdy[1][1] = -(
                144.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 3.0) * std::pow(s_eta, 2.0)
                - 112.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * y[1] * s_eta
                + 96.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * y[1] * s_omega_2 * std::pow(s_eta, 2.0)
                - 12.0 * std::pow(y[1], 3.0) * s_eta
                + 2.0 * std::pow(y[0], 2.0) * y[1]
                - 12.0 * y[1] * s_omega_2 * s_eta
        ) / (denominatorDyDx1)
                     + (-1.0) * numeratorDyDx1 / std::pow(denominatorDyDx1, 2.0) * (
                2.0 * y[0] * (
                        72.0 * std::pow(EulerConst, 2.0 * y[2]) * y[1] * Pi * std::pow(s_eta, 2.0)
                        - 6.0 * y[1] * s_eta
                )
        );

        dfdy[1][2] = -2.0 * std::pow(EulerConst, 2.0 * y[2]) * (
                36.0 * Pi * std::pow(y[1], 4.0) * std::pow(s_eta, 2.0)
                - 56.0 * Pi * std::pow(y[0], 2.0) * std::pow(y[1], 2.0) * s_eta
                + 48.0 * Pi * std::pow(y[1], 2.0) * s_omega_2 * std::pow(s_eta, 2.0)
                + 4.0 * Pi * std::pow(y[0], 4.0)
                - 16.0 * Pi * std::pow(y[0], 2.0) * s_omega_2 * s_eta
                + 12.0 * Pi * std::pow(s_omega_2, 2.0) * std::pow(s_eta, 2.0)
        ) / (denominatorDyDx1)
                     + (-1.0) * numeratorDyDx1 / std::pow(denominatorDyDx1, 2.0) * (
                4.0 * y[0] * std::pow(EulerConst, 2.0 * y[2]) * (
                        36.0 * std::pow(y[1], 2.0) * Pi * std::pow(s_eta, 2.0)
                        + 4.0 * Pi * std::pow(y[0], 2.0) * s_eta
                        - 12.0 * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                )
        );

        dfdy[2][0] = -3.0 * y[1] * (
                16.0 * Pi * y[0] * std::pow(EulerConst, 2.0 * y[2]) * s_eta + 2.0 * y[0]) / (denominatorDyDx2)
                     + numeratorDyDx2 / std::pow(denominatorDyDx2, 2.0) * (
                2.0 * y[0] * y[0]
                + 8.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                + (
                        36.0 * std::pow(EulerConst, 2.0 * y[2]) * std::pow(y[1], 2.0) * Pi * std::pow(s_eta, 2.0)
                        + 4.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                        - 12.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                        - 3.0 * y[1] * s_eta
                        + y[0] * y[0]
                        - 3.0 * s_omega_2 * s_eta
                )
        );

        dfdy[2][1] = -3.0 * (8.0 * Pi * std::pow(y[0], 2.0) * std::pow(EulerConst, 2.0 * y[2]) * s_eta
                             - 8.0 * Pi * std::pow(EulerConst, 2.0 * y[2]) * s_omega_2 * std::pow(s_eta, 2.0)
                             + std::pow(y[0], 2.0)) / (denominatorDyDx2)
                     + numeratorDyDx2 / std::pow(denominatorDyDx2, 2) * (
                y[0] * (
                        72.0 * std::pow(EulerConst, 2.0 * y[2]) * y[1] * Pi * std::pow(s_eta, 2.0)
                        - 3.0 * s_eta
                )
        );

        dfdy[2][2] = -6.0 * y[1] * std::pow(EulerConst, 2.0 * y[2]) * (
                8.0 * Pi * std::pow(y[0], 2.0) * s_eta
                - 8.0 * Pi * s_omega_2 * std::pow(s_eta, 2.0)
        ) / (denominatorDyDx2)
                     + numeratorDyDx2 / std::pow(denominatorDyDx2, 2.0) * (
                2.0 * y[0] * std::pow(EulerConst, 2.0 * y[2]) * (
                        36.0 * std::pow(y[1], 2.0) * Pi * std::pow(s_eta, 2.0)
                        + 4.0 * Pi * std::pow(y[0], 2.0) * s_eta
                        - 12.0 * Pi * s_omega_2 * std::pow(s_eta, 2.0)
                )
        );
    }

};

void get_bianchiIX_isotropy_points();

#endif //ODECALCULATIONS_BIANCHIIX_ISOTROPY_H
