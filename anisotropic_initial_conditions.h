//
// Created by rafkat on 8/27/24.
//

#ifndef ODECALCULATIONS_ANISOTROPIC_INITIAL_CONDITIONS_H
#define ODECALCULATIONS_ANISOTROPIC_INITIAL_CONDITIONS_H

#include "nr/nr3.h"
#include "constants.h"

namespace anisotropic_initial_conditions {
    Doub a_iv = 1.0;
    Doub da_iv = 1.0;
    Doub omega_2 = 1.0e-3;
    Doub eta = 18.0 * 1.0e-3;
    Doub s_iv = 1.0e-3;

    Doub phi_iv =
            0.5 * std::log(3.0 / 4.0 / Pi * (da_iv * da_iv / a_iv / a_iv - s_iv * s_iv + omega_2 / a_iv / a_iv) /
                           (1.0 - 3.0 * eta * (3.0 * da_iv * da_iv / a_iv / a_iv - 3.0 * s_iv * s_iv + omega_2 / a_iv / a_iv)));
}


namespace ode_parameters {
    const Int anisotropic_nvar{4};
    const Doub atol{1.0e-11};
    const Doub rtol{1.0e-11};
    const Doub h1{1.0e-5};
    const Doub hmin{1.0e-25};
    const Doub x1{0.0};
}

namespace bIX_anisotropic_up {
    const Doub x2_bw{-2.44};
    const Doub x2_fw{18.4};
    VecDoub ystart(ode_parameters::anisotropic_nvar);
}

namespace bIX_anisotropic_down {
    const Doub x2_fw{-bIX_anisotropic_up::x2_bw};
    const Doub x2_bw{-bIX_anisotropic_up::x2_fw};
    VecDoub ystart(ode_parameters::anisotropic_nvar);
}


#endif //ODECALCULATIONS_ANISOTROPIC_INITIAL_CONDITIONS_H
