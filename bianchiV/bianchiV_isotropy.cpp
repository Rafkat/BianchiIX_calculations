#include "../initial_conditions/BV_isotropic_initial_conditions.h"
#include "bianchiV_isotropy.h"
#include "../nr/stepperbs.h"
#include "../output/output_writer.h"

void get_bianchiV_isotropy_points(int n_points) {
    BV_isotropic_up::ystart[0] = BV_isotropic_initial_conditions::a_iv;
    BV_isotropic_up::ystart[1] = BV_isotropic_initial_conditions::da_iv;
    BV_isotropic_up::ystart[2] = BV_isotropic_initial_conditions::phi_iv;

    BV_isotropic_down::ystart[0] = BV_isotropic_initial_conditions::a_iv;
    BV_isotropic_down::ystart[1] = -BV_isotropic_initial_conditions::da_iv;
    BV_isotropic_down::ystart[2] = BV_isotropic_initial_conditions::phi_iv;

    bianchiV_rhs d(BV_isotropic_initial_conditions::eta, BV_isotropic_initial_conditions::omega_2);

    std::vector<std::vector<Doub>> ip{{BV_ode_parameters::x1, BV_isotropic_up::x2_bw},
                                      {BV_ode_parameters::x1, BV_isotropic_up::x2_fw},
                                      {BV_ode_parameters::x1, BV_isotropic_down::x2_bw},
                                      {BV_ode_parameters::x1, BV_isotropic_down::x2_fw}};

    std::vector<VecDoub> ystarts{BV_isotropic_up::ystart, BV_isotropic_up::ystart,
                                 BV_isotropic_down::ystart, BV_isotropic_down::ystart};

    std::vector<string_view> file_names{"../bianchiV_" + std::to_string(n_points) + "p_fw.txt",
                                        "../bianchiV_" + std::to_string(n_points) + "p_fw.txt",
                                        "../bianchiV_" + std::to_string(n_points) + "p_bw.txt",
                                        "../bianchiV_" + std::to_string(n_points) + "p_bw.txt"};

    ofstream bianchiV_upside;
    ofstream bianchiV_downside;

    bianchiV_upside.open("../output/bianchiV_isotropic_eta_" + std::to_string(BV_isotropic_initial_conditions::eta)
                          + "_omega_2_" + std::to_string(BV_isotropic_initial_conditions::omega_2)
                          + "_" + std::to_string(n_points) + "p_upside.txt");
    bianchiV_downside.open("../output/bianchiV_isotropic_eta_" + std::to_string(BV_isotropic_initial_conditions::eta)
                            + "_omega_2_" + std::to_string(BV_isotropic_initial_conditions::omega_2)
                            + "_" + std::to_string(n_points) + "p_downside.txt");

    for (std::size_t i{0}; i < ip.size(); i++) {
        Output out(n_points);
        Odeint<StepperBS<bianchiV_rhs> > ode(ystarts[i],
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
                write_to_file(out, bianchiV_upside, true);
                break;
            case 1:
                write_to_file(out, bianchiV_upside, false);
                break;
            case 2:
                write_to_file(out, bianchiV_downside, true);
                break;
            case 3:
                write_to_file(out, bianchiV_downside, false);
                break;
            default:
                break;
        }

    }

    bianchiV_downside.close();
    bianchiV_upside.close();
}