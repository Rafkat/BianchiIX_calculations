#include "../initial_conditions/BIX_isotropic_initial_conditions.h"
#include "bianchiIX_isotropy.h"
#include "../nr/stepperbs.h"
#include "../output/output_writer.h"

void get_bianchiIX_isotropy_points(int n_points) {
    BIX_isotropic_up::ystart[0] = BIX_isotropic_initial_conditions::a_iv;
    BIX_isotropic_up::ystart[1] = BIX_isotropic_initial_conditions::da_iv;
    BIX_isotropic_up::ystart[2] = BIX_isotropic_initial_conditions::phi_iv;

    BIX_isotropic_down::ystart[0] = -BIX_isotropic_initial_conditions::a_iv;
    BIX_isotropic_down::ystart[1] = BIX_isotropic_initial_conditions::da_iv;
    BIX_isotropic_down::ystart[2] = BIX_isotropic_initial_conditions::phi_iv;

    bianchiIX_rhs d(BIX_isotropic_initial_conditions::eta, BIX_isotropic_initial_conditions::omega_2);

    std::vector<std::vector<Doub>> ip{{BIX_ode_parameters::x1, BIX_isotropic_up::x2_bw},
                                      {BIX_ode_parameters::x1, BIX_isotropic_up::x2_fw},
                                      {BIX_ode_parameters::x1, BIX_isotropic_down::x2_bw},
                                      {BIX_ode_parameters::x1, BIX_isotropic_down::x2_fw}};

    std::vector<VecDoub> ystarts{BIX_isotropic_up::ystart, BIX_isotropic_up::ystart,
                                 BIX_isotropic_down::ystart, BIX_isotropic_down::ystart};

    std::vector<string_view> file_names{"../bianchiIX_" + std::to_string(n_points) + "p_fw.txt",
                                        "../bianchiIX_" + std::to_string(n_points) + "p_fw.txt",
                                        "../bianchiIX_" + std::to_string(n_points) + "p_bw.txt",
                                        "../bianchiIX_" + std::to_string(n_points) + "p_bw.txt"};

    ofstream bianchiIX_upside;
    ofstream bianchiIX_downside;

    bianchiIX_upside.open("../output/bianchiIX_isotropic_eta_" + std::to_string(BIX_isotropic_initial_conditions::eta)
                          + "_omega_2_" + std::to_string(BIX_isotropic_initial_conditions::omega_2)
                          + "_" + std::to_string(n_points) + "p_upside.txt");
    bianchiIX_downside.open("../output/bianchiIX_isotropic_eta_" + std::to_string(BIX_isotropic_initial_conditions::eta)
                            + "_omega_2_" + std::to_string(BIX_isotropic_initial_conditions::omega_2)
                            + "_" + std::to_string(n_points) + "p_downside.txt");

    for (std::size_t i{0}; i < ip.size(); i++) {
        Output out(n_points);
        Odeint<StepperBS<bianchiIX_rhs> > ode(ystarts[i],
                                              ip[i][0],
                                              ip[i][1],
                                              BIX_ode_parameters::atol,
                                              BIX_ode_parameters::rtol,
                                              BIX_ode_parameters::h1,
                                              BIX_ode_parameters::hmin,
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