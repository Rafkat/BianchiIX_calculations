#include "../initial_conditions/BIX_anisotropic_initial_conditions.h"
#include "bianchiIX_anisotropy.h"
//#include "../nr/stepperross.h"
//#include "../nr/stepperdopr853.h"
#include "../nr/stepperbs.h"
#include "../output/output_writer.h"


void get_bianchiIX_anisotropy_points(int n_points) {
    bIX_anisotropic_up::ystart[0] = BIX_anisotropic_initial_conditions::a_iv;
    bIX_anisotropic_up::ystart[1] = BIX_anisotropic_initial_conditions::da_iv;
    bIX_anisotropic_up::ystart[2] = BIX_anisotropic_initial_conditions::phi_iv;
    bIX_anisotropic_up::ystart[3] = BIX_anisotropic_initial_conditions::dbp_iv;
    bIX_anisotropic_up::ystart[4] = BIX_anisotropic_initial_conditions::dbm_iv;
    bIX_anisotropic_up::ystart[5] = BIX_anisotropic_initial_conditions::bp_iv;
    bIX_anisotropic_up::ystart[6] = BIX_anisotropic_initial_conditions::bm_iv;

    bIX_anisotropic_down::ystart[0] = BIX_anisotropic_initial_conditions::a_iv;
    bIX_anisotropic_down::ystart[1] = -BIX_anisotropic_initial_conditions::da_iv;
    bIX_anisotropic_down::ystart[2] = BIX_anisotropic_initial_conditions::phi_iv;
    bIX_anisotropic_down::ystart[3] = -BIX_anisotropic_initial_conditions::dbp_iv;
    bIX_anisotropic_down::ystart[4] = -BIX_anisotropic_initial_conditions::dbm_iv;
    bIX_anisotropic_down::ystart[5] = -BIX_anisotropic_initial_conditions::bp_iv;
    bIX_anisotropic_down::ystart[6] = -BIX_anisotropic_initial_conditions::bm_iv;

    bianchiIX_anisotropic_rhs d(BIX_anisotropic_initial_conditions::eta, BIX_anisotropic_initial_conditions::omega_2);

    std::vector<std::vector<Doub>> ip{{BIX_ode_parameters::x1, bIX_anisotropic_up::x2_bw},
                                      {BIX_ode_parameters::x1, bIX_anisotropic_up::x2_fw},
                                      {BIX_ode_parameters::x1, bIX_anisotropic_down::x2_bw},
                                      {BIX_ode_parameters::x1, bIX_anisotropic_down::x2_fw}};

    std::vector<VecDoub> ystarts{bIX_anisotropic_up::ystart, bIX_anisotropic_up::ystart,
                                 bIX_anisotropic_down::ystart, bIX_anisotropic_down::ystart};

    std::vector<string_view> file_names{"../bianchiIX_" + std::to_string(n_points) + "p_fw.txt",
                                        "../bianchiIX_" + std::to_string(n_points) + "p_fw.txt",
                                        "../bianchiIX_" + std::to_string(n_points) + "p_bw.txt",
                                        "../bianchiIX_" + std::to_string(n_points) + "p_bw.txt"};

    ofstream bianchiIX_upside;
    ofstream bianchiIX_downside;

    bianchiIX_upside.open(
            "../output/bianchiIX_anisotropic_eta_" + std::to_string(BIX_anisotropic_initial_conditions::eta)
            + "_omega_2_" + std::to_string(BIX_anisotropic_initial_conditions::omega_2)
            + "_dbp_" + std::to_string(BIX_anisotropic_initial_conditions::dbp_iv)
            + "_dbm_" + std::to_string(BIX_anisotropic_initial_conditions::dbm_iv)
            + "_" + std::to_string(n_points) + "p_upside.txt");
    bianchiIX_downside.open(
            "../output/bianchiIX_anisotropic_eta_" + std::to_string(BIX_anisotropic_initial_conditions::eta)
            + "_omega_2_" + std::to_string(BIX_anisotropic_initial_conditions::omega_2)
            + "_dbp_" + std::to_string(BIX_anisotropic_initial_conditions::dbp_iv)
            + "_dbm_" + std::to_string(BIX_anisotropic_initial_conditions::dbm_iv)
            + "_" + std::to_string(n_points) + "p_downside.txt");

    for (std::size_t i{0}; i < ip.size(); i++) {
        Output out(n_points);
        Odeint<StepperBS<bianchiIX_anisotropic_rhs> > ode(ystarts[i],
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