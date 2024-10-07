//
// Created by rafkat on 8/27/24.
//

#ifndef ODECALCULATIONS_BIX_ANISOTROPIC_INITIAL_CONDITIONS_H
#define ODECALCULATIONS_BIX_ANISOTROPIC_INITIAL_CONDITIONS_H

#include "../nr/nr3.h"
#include "../constants.h"

namespace BIX_anisotropic_initial_conditions {
    Doub a_iv = 1.0;
    Doub da_iv = 1.0;
    Doub omega_2 = 1.0e-5;
    Doub eta = 18.0 * 1.0e-3;
    Doub dbp_iv = 1.0e-3;
    Doub dbm_iv = 1.0e-3;
    Doub bp_iv = 1.0e-4;
    Doub bm_iv = 1.0e-4;
    Doub KIX = -0.333333333333333e0 * exp(0.40e1 * bp_iv + 0.692820323027552e1 * bm_iv) -
               0.333333333333333e0 * exp(0.40e1 * bp_iv - 0.692820323027552e1 * bm_iv) -
               0.333333333333333e0 * exp(-0.80e1 * bp_iv) +
               0.666666666666666e0 * exp(-0.20e1 * bp_iv + 0.346410161513776e1 * bm_iv) +
               0.666666666666666e0 * exp(-0.20e1 * bp_iv - 0.346410161513776e1 * bm_iv) +
               0.666666666666666e0 * exp(0.40e1 * bp_iv);

    Doub phi_iv =
            0.5 * std::log(3.0 / 4.0 / Pi *
                           (da_iv * da_iv / a_iv / a_iv - dbp_iv * dbp_iv - dbm_iv * dbm_iv + omega_2 * KIX / a_iv / a_iv) /
                           (1.0 - 3.0 * eta *
                                  (3.0 * da_iv * da_iv / a_iv / a_iv - 3.0 * dbp_iv * dbp_iv - 3.0 * dbm_iv * dbm_iv +
                                   omega_2 * KIX/ a_iv / a_iv)));
}


namespace BIX_ode_parameters {
    const Int anisotropic_nvar{7};
    const Doub atol{1.0e-11};
    const Doub rtol{1.0e-11};
    const Doub h1{1.0e-5};
    const Doub hmin{1.0e-25};
    const Doub x1{0.0};
}

namespace bIX_anisotropic_up {
    const Doub x2_bw{-3.421218299};
    const Doub x2_fw{19.7477};
    VecDoub ystart(BIX_ode_parameters::anisotropic_nvar);
}

namespace bIX_anisotropic_down {
    const Doub x2_fw{-bIX_anisotropic_up::x2_bw};
    const Doub x2_bw{-bIX_anisotropic_up::x2_fw};
    VecDoub ystart(BIX_ode_parameters::anisotropic_nvar);
}


#endif //ODECALCULATIONS_BIX_ANISOTROPIC_INITIAL_CONDITIONS_H
