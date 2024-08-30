#include "../anisotropic_initial_conditions.h"
#include "bianchiIX_anisotropy.h"
//#include "../nr/stepperross.h"
//#include "../nr/stepperdopr853.h"
#include "../nr/stepperbs.h"
#include "../output/output_writer.h"


void get_bianchiIX_anisotropy_points() {
    bIX_anisotropic_up::ystart[0] = anisotropic_initial_conditions::a_iv;
    bIX_anisotropic_up::ystart[1] = anisotropic_initial_conditions::da_iv;
    bIX_anisotropic_up::ystart[2] = anisotropic_initial_conditions::phi_iv;
    bIX_anisotropic_up::ystart[3] = anisotropic_initial_conditions::s_iv;

    bIX_anisotropic_down::ystart[0] = anisotropic_initial_conditions::a_iv;
    bIX_anisotropic_down::ystart[1] = -anisotropic_initial_conditions::da_iv;
    bIX_anisotropic_down::ystart[2] = anisotropic_initial_conditions::phi_iv;
    bIX_anisotropic_down::ystart[3] = -anisotropic_initial_conditions::s_iv;

    int n_points{300};

    bianchiIX_anisotropic_rhs d(anisotropic_initial_conditions::eta, anisotropic_initial_conditions::omega_2);

    std::vector<std::vector<Doub>> ip{{ode_parameters::x1, bIX_anisotropic_up::x2_bw},
                                      {ode_parameters::x1, bIX_anisotropic_up::x2_fw},
                                      {ode_parameters::x1, bIX_anisotropic_down::x2_bw},
                                      {ode_parameters::x1, bIX_anisotropic_down::x2_fw}};

    std::vector<VecDoub> ystarts{bIX_anisotropic_up::ystart, bIX_anisotropic_up::ystart,
                                 bIX_anisotropic_down::ystart, bIX_anisotropic_down::ystart};

    std::vector<string_view> file_names{"../bianchiIX_" + std::to_string(n_points) + "p_fw.txt",
                                        "../bianchiIX_" + std::to_string(n_points) + "p_fw.txt",
                                        "../bianchiIX_" + std::to_string(n_points) + "p_bw.txt",
                                        "../bianchiIX_" + std::to_string(n_points) + "p_bw.txt"};

    ofstream bianchiIX_upside;
    ofstream bianchiIX_downside;

    bianchiIX_upside.open("../output/bianchiIX_anisotropic_" + std::to_string(n_points) + "p_upside.txt");
    bianchiIX_downside.open("../output/bianchiIX_anisotropic_" + std::to_string(n_points) + "p_downside.txt");

    for (std::size_t i{0}; i < ip.size(); i++) {
        Output out(n_points);
        Odeint<StepperBS<bianchiIX_anisotropic_rhs> > ode(ystarts[i],
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
                write_to_file(out, bianchiIX_upside, true, true);
                break;
            case 1:
                write_to_file(out, bianchiIX_upside, false, true);
                break;
            case 2:
                write_to_file(out, bianchiIX_downside, true, true);
                break;
            case 3:
                write_to_file(out, bianchiIX_downside, false, true);
                break;
            default:
                break;
        }

    }

    bianchiIX_downside.close();
    bianchiIX_upside.close();
}