#include "../initial_conditions/BV_anisotropic_initial_conditions.h"
#include "bianchiV_anisotropy.h"
#include "../nr/stepperbs.h"
#include "../output/output_writer.h"


void get_bianchiV_anisotropy_points(int n_points) {
    BV_anisotropic_up::ystart[0] = BV_anisotropic_initial_conditions::a_iv;
    BV_anisotropic_up::ystart[1] = BV_anisotropic_initial_conditions::da_iv;
    BV_anisotropic_up::ystart[2] = BV_anisotropic_initial_conditions::phi_iv;
    BV_anisotropic_up::ystart[3] = BV_anisotropic_initial_conditions::s_iv;

    BV_anisotropic_down::ystart[0] = BV_anisotropic_initial_conditions::a_iv;
    BV_anisotropic_down::ystart[1] = -BV_anisotropic_initial_conditions::da_iv;
    BV_anisotropic_down::ystart[2] = BV_anisotropic_initial_conditions::phi_iv;
    BV_anisotropic_down::ystart[3] = BV_anisotropic_initial_conditions::s_iv;

    bianchiV_anisotropic_rhs d(BV_anisotropic_initial_conditions::eta, BV_anisotropic_initial_conditions::omega_2);

    std::vector<std::vector<Doub>> ip{{BV_ode_parameters::x1, BV_anisotropic_up::x2_bw},
                                      {BV_ode_parameters::x1, BV_anisotropic_up::x2_fw},
                                      {BV_ode_parameters::x1, BV_anisotropic_down::x2_bw},
                                      {BV_ode_parameters::x1, BV_anisotropic_down::x2_fw}};

    std::vector<VecDoub> ystarts{BV_anisotropic_up::ystart, BV_anisotropic_up::ystart,
                                 BV_anisotropic_down::ystart, BV_anisotropic_down::ystart};

    std::vector<string_view> file_names{"../bianchiV_" + std::to_string(n_points) + "p_fw.txt",
                                        "../bianchiV_" + std::to_string(n_points) + "p_fw.txt",
                                        "../bianchiV_" + std::to_string(n_points) + "p_bw.txt",
                                        "../bianchiV_" + std::to_string(n_points) + "p_bw.txt"};

    ofstream bianchiV_upside;
    ofstream bianchiV_downside;

    bianchiV_upside.open("../output/bianchiV_anisotropic_eta_" + std::to_string(BV_anisotropic_initial_conditions::eta)
                          + "_omega_2_" + std::to_string(BV_anisotropic_initial_conditions::omega_2)
                          + "_s_" + std::to_string(BV_anisotropic_initial_conditions::s_iv)
                          + "_" + std::to_string(n_points) + "p_upside.txt");
    bianchiV_downside.open("../output/bianchiV_anisotropic_eta_" + std::to_string(BV_anisotropic_initial_conditions::eta)
                            + "_omega_2_" + std::to_string(BV_anisotropic_initial_conditions::omega_2)
                            + "_s_" + std::to_string(BV_anisotropic_initial_conditions::s_iv)
                            + "_" + std::to_string(n_points) + "p_downside.txt");

    for (std::size_t i{0}; i < ip.size(); i++) {
        Output out(n_points);
        Odeint<StepperBS<bianchiV_anisotropic_rhs> > ode(ystarts[i],
                                                          ip[i][0],
                                                          ip[i][1],
                                                          BV_ode_parameters::atol,
                                                          BV_ode_parameters::rtol,
                                                          BV_ode_parameters::h1,
                                                          BV_ode_parameters::hmin,
                                                          out,
                                                          d);
        ode.integrate();

        switch (i) {
            case 0:
                write_to_file(out, bianchiV_upside, true, true);
                break;
            case 1:
                write_to_file(out, bianchiV_upside, false, true);
                break;
            case 2:
                write_to_file(out, bianchiV_downside, true, true);
                break;
            case 3:
                write_to_file(out, bianchiV_downside, false, true);
                break;
            default:
                break;
        }

    }

    bianchiV_downside.close();
    bianchiV_upside.close();
}