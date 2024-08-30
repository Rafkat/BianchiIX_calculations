//
// Created by rafkat on 8/19/24.
//

#ifndef ODECALCULATIONS_ISOTROPIC_INITIAL_CONDITIONS_H
#define ODECALCULATIONS_ISOTROPIC_INITIAL_CONDITIONS_H

#include "nr/nr3.h"
#include "constants.h"

namespace isotropic_initial_conditions {
    Doub a_iv = 1.0;
    Doub da_iv = 1.0;
    Doub omega_2 = 1.0e-3;
    Doub eta = 18.0 * 1.0e-3;
    Doub s_iv = 1.0e-3;

    Doub phi_iv = 0.5 * std::log(3.0 / 4.0 / Pi * (da_iv * da_iv / a_iv / a_iv + omega_2 / a_iv / a_iv) /
                                 (1 - 3 * eta * (3 * da_iv * da_iv / a_iv / a_iv + omega_2 / a_iv / a_iv)));
}


namespace ode_parameters {
    const Int isotropic_nvar{3};
    const Int anisotropic_nvar{4};
    const Doub atol{1.0e-11};
    const Doub rtol{1.0e-11};
    const Doub h1{1.0e-8};
    const Doub hmin{1.0e-19};
    const Doub x1{0.0};
}

namespace bIX_isotropic_up {
    const Doub x2_bw{-2.362292};
    const Doub x2_fw{5.0};
    VecDoub ystart(ode_parameters::isotropic_nvar);
}

namespace bIX_isotropic_down {
    const Doub x2_fw{-bIX_isotropic_up::x2_bw};
    const Doub x2_bw{-bIX_isotropic_up::x2_fw};
    VecDoub ystart(ode_parameters::isotropic_nvar);
}


#endif //ODECALCULATIONS_ISOTROPIC_INITIAL_CONDITIONS_H
