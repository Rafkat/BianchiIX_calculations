#include "../isotropic_initial_conditions.h"
#include "bianchiIX_isotropy.h"
#include "../nr/stepperross.h"
#include "../output/output_writer.h"

void get_bianchiIX_isotropy_points() {
    bIX_isotropic_up::ystart[0] = isotropic_initial_conditions::a_iv;
    bIX_isotropic_up::ystart[1] = isotropic_initial_conditions::da_iv;
    bIX_isotropic_up::ystart[2] = isotropic_initial_conditions::phi_iv;

    bIX_isotropic_down::ystart[0] = -isotropic_initial_conditions::a_iv;
    bIX_isotropic_down::ystart[1] = isotropic_initial_conditions::da_iv;
    bIX_isotropic_down::ystart[2] = isotropic_initial_conditions::phi_iv;

    int n_points{100};

    bianchiIX_rhs d(isotropic_initial_conditions::eta, isotropic_initial_conditions::omega_2);

    std::vector<std::vector<Doub>> ip{{ode_parameters::x1, bIX_isotropic_up::x2_bw},
                                      {ode_parameters::x1, bIX_isotropic_up::x2_fw},
                                      {ode_parameters::x1, bIX_isotropic_down::x2_bw},
                                      {ode_parameters::x1, bIX_isotropic_down::x2_fw}};

    std::vector<VecDoub> ystarts{bIX_isotropic_up::ystart, bIX_isotropic_up::ystart,
                                 bIX_isotropic_down::ystart, bIX_isotropic_down::ystart};

    std::vector<string_view> file_names{"../bianchiIX_" + std::to_string(n_points) + "p_fw.txt",
                                        "../bianchiIX_" + std::to_string(n_points) + "p_fw.txt",
                                        "../bianchiIX_" + std::to_string(n_points) + "p_bw.txt",
                                        "../bianchiIX_" + std::to_string(n_points) + "p_bw.txt"};

    ofstream bianchiIX_upside;
    ofstream bianchiIX_downside;

    bianchiIX_upside.open("../output/bianchiIX_isotropic_" + std::to_string(n_points) + "p_upside.txt");
    bianchiIX_downside.open("../output/bianchiIX_isotropic" + std::to_string(n_points) + "p_downside.txt");

    for (std::size_t i{0}; i < ip.size(); i++) {
        Output out(n_points);
        Odeint<StepperRoss<bianchiIX_rhs> > ode(ystarts[i],
                                                ip[i][0],
                                                ip[i][1],
                                                ode_parameters::atol,
                                                ode_parameters::rtol,
                                                ode_parameters::h1,
                                                ode_parameters::hmin,
                                                out,
                                                d);
        ode.integrate();

        switch (i) {
            case 0:
                write_to_file(out, bianchiIX_upside, true);
                break;
            case 1:
                write_to_file(out, bianchiIX_upside, false);
                break;
            case 2:
                write_to_file(out, bianchiIX_downside, true);
                break;
            case 3:
                write_to_file(out, bianchiIX_downside, false);
                break;
            default:
                break;
        }

    }

    bianchiIX_downside.close();
    bianchiIX_upside.close();
}