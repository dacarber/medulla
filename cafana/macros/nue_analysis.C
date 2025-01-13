/**
 * @file analysis.C
 * @brief The main analysis macro for the muon2024 analysis.
 * @details This macro drives the analysis by configuring the variables, cuts,
 * and samples to be used in the analysis. This is accomplished through the use
 * of the Analysis class, which containerizes the configuration of the analysis
 * and reduces the amount of boilerplate code needed to run the analysis.
 * @author mueller@fnal.gov
*/
#include "include/variables.h"
#include "include/muon2024/variables_muon2024.h"
#include "include/cuts.h"
#include "include/muon2024/cuts_muon2024.h"
#include "include/preprocessor.h"
#include "include/analysis.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

#include <algorithm>

void analysis()
{
    ana::Analysis analysis("muon2024_1muNp_data");

    ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/mc_v09_84_00_01/flat/*.root");
    analysis.AddLoader("mc", &mc, true);

    ana::SpectrumLoader onbeam("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/data_v09_84_00_01/onbeam/flat/*.root");
    analysis.AddLoader("onbeam", &onbeam, false);

    ana::SpectrumLoader offbeam("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/data_v09_84_00_01/offbeam/flat/*.root");
    analysis.AddLoader("offbeam", &offbeam, false);

    ana::SpectrumLoader var00("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var00_nominal.flat.root");
    analysis.AddLoader("nominal", &var00, true);
    
    ana::SpectrumLoader var01("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var01_untunedtpcsigshape.flat.root");
    analysis.AddLoader("var01_untunedtpcsigshape", &var01, true);

    ana::SpectrumLoader var03low("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var03_tpcind1decreasegain.flat.root");
    analysis.AddLoader("var03_tpcind1decreasegain", &var03low, true);

    ana::SpectrumLoader var03high("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var03_tpcind1increasegain.flat.root");
    analysis.AddLoader("var03_tpcind1increasegain", &var03high, true);

    ana::SpectrumLoader var04("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var04_pmtdecreaseqe.flat.root");
    analysis.AddLoader("var04_pmtdecreaseqe", &var04, true);

    ana::SpectrumLoader var05("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var05_ellipsoidalrecomb.flat.root");
    analysis.AddLoader("var05_ellipsoidalrecomb", &var05, true);

    ana::SpectrumLoader var06low("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var06_tpccohnoisem1sigma.flat.root");
    analysis.AddLoader("var06_tpccohnoisem1sigma", &var06low, true);

    ana::SpectrumLoader var06high("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var06_tpccohnoisep1sigma.flat.root");
    analysis.AddLoader("var06_tpccohnoisep1sigma", &var06high, true);

    ana::SpectrumLoader var07low("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var07_tpcintnoisem1sigma.flat.root");
    analysis.AddLoader("var07_tpcintnoisem1sigma", &var07low, true);

    ana::SpectrumLoader var07high("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var07_tpcintnoisep1sigma.flat.root");
    analysis.AddLoader("var07_tpcintnoisep1sigma", &var07high, true);

    ana::SpectrumLoader var08low("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var08_tpclowlifetime.flat.root");
    analysis.AddLoader("var08_tpclowlifetime", &var08low, true);

    ana::SpectrumLoader var08high("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var08_tpchighlifetime.flat.root");
    analysis.AddLoader("var08_tpchighlifetime", &var08high, true);

    ana::SpectrumLoader var09("/pnfs/icarus/persistent/users/mueller/spinereco2024/allplanes/detsys_v09_89_01_01/var09_null.flat.root");
    analysis.AddLoader("var09_null", &var09, true);

    /**
     * @brief Add a set of variables for selected interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define CUT cuts::muon2024::all_1muNp_cut
    std::map<std::string, ana::SpillMultiVar> vars_selected;
    vars_selected.insert({"nu_id", ana::SpillMultiVar(SPINEVAR_RT(vars::neutrino_id, CUT, cuts::no_cut))});
    vars_selected.insert({"baseline", ana::SpillMultiVar(SPINEVAR_RT(vars::true_neutrino_baseline, CUT, cuts::no_cut))});
    vars_selected.insert({"pdg", ana::SpillMultiVar(SPINEVAR_RT(vars::true_neutrino_pdg, CUT, cuts::no_cut))});
    vars_selected.insert({"cc", ana::SpillMultiVar(SPINEVAR_RT(vars::true_neutrino_cc, CUT, cuts::no_cut))});
    vars_selected.insert({"category", ana::SpillMultiVar(SPINEVAR_RT(vars::muon2024::category, CUT, cuts::no_cut))});
    vars_selected.insert({"interaction_mode", ana::SpillMultiVar(SPINEVAR_RT(vars::neutrino_interaction_mode, CUT, cuts::no_cut))});
    vars_selected.insert({"true_edep", ana::SpillMultiVar(SPINEVAR_RT(vars::true_neutrino_energy, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_edep", ana::SpillMultiVar(SPINEVAR_RR(vars::visible_energy, CUT, cuts::no_cut))});
    vars_selected.insert({"true_muon_x", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_end_x, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_muon_x", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_end_x, CUT, cuts::no_cut))});
    vars_selected.insert({"true_muon_y", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_end_y, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_muon_y", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_end_y, CUT, cuts::no_cut))});
    vars_selected.insert({"true_muon_z", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_end_z, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_muon_z", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_end_z, CUT, cuts::no_cut))});
    vars_selected.insert({"true_proton_x", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_end_x, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_proton_x", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_end_x, CUT, cuts::no_cut))});
    vars_selected.insert({"true_proton_y", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_end_y, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_proton_y", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_end_y, CUT, cuts::no_cut))});
    vars_selected.insert({"true_proton_z", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_end_z, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_proton_z", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_end_z, CUT, cuts::no_cut))});
    vars_selected.insert({"true_tmuon", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_ke, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_tmuon", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_ke, CUT, cuts::no_cut))});
    vars_selected.insert({"true_tproton", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_ke, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_tproton", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_ke, CUT, cuts::no_cut))});
    vars_selected.insert({"true_ptmuon", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_muon_pt, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_ptmuon", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_pt, CUT, cuts::no_cut))});
    vars_selected.insert({"true_ptproton", ana::SpillMultiVar(SPINEVAR_RT(vars::leading_proton_pt, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_ptproton", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_pt, CUT, cuts::no_cut))});
    vars_selected.insert({"true_theta_mu", ana::SpillMultiVar(SPINEVAR_RT(vars::muon_polar_angle, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_theta_mu", ana::SpillMultiVar(SPINEVAR_RR(vars::muon_polar_angle, CUT, cuts::no_cut))});
    vars_selected.insert({"true_phi_mu", ana::SpillMultiVar(SPINEVAR_RT(vars::muon_azimuthal_angle, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_phi_mu", ana::SpillMultiVar(SPINEVAR_RR(vars::muon_azimuthal_angle, CUT, cuts::no_cut))});
    vars_selected.insert({"true_opening_angle", ana::SpillMultiVar(SPINEVAR_RT(vars::muon2024::opening_angle, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_opening_angle", ana::SpillMultiVar(SPINEVAR_RR(vars::muon2024::opening_angle, CUT, cuts::no_cut))});
    vars_selected.insert({"true_dpT", ana::SpillMultiVar(SPINEVAR_RT(vars::interaction_pt, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_dpT", ana::SpillMultiVar(SPINEVAR_RR(vars::interaction_pt, CUT, cuts::no_cut))});
    vars_selected.insert({"true_dphiT", ana::SpillMultiVar(SPINEVAR_RT(vars::phiT, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_dphiT", ana::SpillMultiVar(SPINEVAR_RR(vars::phiT, CUT, cuts::no_cut))});
    vars_selected.insert({"true_edalphaT", ana::SpillMultiVar(SPINEVAR_RT(vars::alphaT, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_edalphaT", ana::SpillMultiVar(SPINEVAR_RR(vars::alphaT, CUT, cuts::no_cut))});
    vars_selected.insert({"true_vertex_x", ana::SpillMultiVar(SPINEVAR_RT(vars::vertex_x, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_vertex_x", ana::SpillMultiVar(SPINEVAR_RR(vars::vertex_x, CUT, cuts::no_cut))});
    vars_selected.insert({"true_vertex_y", ana::SpillMultiVar(SPINEVAR_RT(vars::vertex_y, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_vertex_y", ana::SpillMultiVar(SPINEVAR_RR(vars::vertex_y, CUT, cuts::no_cut))});
    vars_selected.insert({"true_vertex_z", ana::SpillMultiVar(SPINEVAR_RT(vars::vertex_z, CUT, cuts::no_cut))});
    vars_selected.insert({"reco_vertex_z", ana::SpillMultiVar(SPINEVAR_RR(vars::vertex_z, CUT, cuts::no_cut))});
    vars_selected.insert({"muon_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_softmax, CUT, cuts::no_cut))});
    vars_selected.insert({"proton_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_proton_softmax, CUT, cuts::no_cut))});
    vars_selected.insert({"mip_softmax", ana::SpillMultiVar(SPINEVAR_RR(vars::leading_muon_mip_softmax, CUT, cuts::no_cut))});
    vars_selected.insert({"flash_time", ana::SpillMultiVar(SPINEVAR_RR(vars::flash_time, CUT, cuts::no_cut))});
    vars_selected.insert({"flash_total", ana::SpillMultiVar(SPINEVAR_RR(vars::flash_total_pe, CUT, cuts::no_cut))});
    vars_selected.insert({"flash_hypothesis", ana::SpillMultiVar(SPINEVAR_RR(vars::flash_hypothesis, CUT, cuts::no_cut))});

    analysis.AddTree("selectedNu", vars_selected, false);    

    /**
     * @brief Add a set of variables for signal interactions to the analysis.
     * @details This adds a set of variables to the analysis by creating a
     * map of variable names and SpillMultiVars that provide the functionality
     * to calculate the variables. These names are used in the TTree that is
     * created by the Tree class to store the results of the analysis.
     */
    #define SIGCUT cuts::muon2024::signal_1muNp
    std::map<std::string, ana::SpillMultiVar> vars_signal;
    vars_signal.insert({"nu_id", ana::SpillMultiVar(SPINEVAR_TT(vars::neutrino_id, SIGCUT))});
    vars_signal.insert({"baseline", ana::SpillMultiVar(SPINEVAR_TT(vars::true_neutrino_baseline, SIGCUT))});
    vars_signal.insert({"pdg", ana::SpillMultiVar(SPINEVAR_TT(vars::true_neutrino_pdg, SIGCUT))});
    vars_signal.insert({"cc", ana::SpillMultiVar(SPINEVAR_TT(vars::true_neutrino_cc, SIGCUT))});
    vars_signal.insert({"category", ana::SpillMultiVar(SPINEVAR_TT(vars::muon2024::category, SIGCUT))});
    vars_signal.insert({"interaction_mode", ana::SpillMultiVar(SPINEVAR_TT(vars::neutrino_interaction_mode, SIGCUT))});
    vars_signal.insert({"true_edep", ana::SpillMultiVar(SPINEVAR_TT(vars::true_neutrino_energy, SIGCUT))});
    vars_signal.insert({"true_muon_x", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_muon_end_x, SIGCUT))});
    vars_signal.insert({"true_muon_y", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_muon_end_y, SIGCUT))});
    vars_signal.insert({"true_muon_z", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_muon_end_z, SIGCUT))});
    vars_signal.insert({"true_proton_x", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_proton_end_x, SIGCUT))});
    vars_signal.insert({"true_proton_y", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_proton_end_y, SIGCUT))});
    vars_signal.insert({"true_proton_z", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_proton_end_z, SIGCUT))});
    vars_signal.insert({"true_tmuon", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_muon_ke, SIGCUT))});
    vars_signal.insert({"true_tproton", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_proton_ke, SIGCUT))});
    vars_signal.insert({"true_ptmuon", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_muon_pt, SIGCUT))});
    vars_signal.insert({"true_ptproton", ana::SpillMultiVar(SPINEVAR_TT(vars::leading_proton_pt, SIGCUT))});
    vars_signal.insert({"true_theta_mu", ana::SpillMultiVar(SPINEVAR_TT(vars::muon_polar_angle, SIGCUT))});
    vars_signal.insert({"true_phi_mu", ana::SpillMultiVar(SPINEVAR_TT(vars::muon_azimuthal_angle, SIGCUT))});
    vars_signal.insert({"true_opening_angle", ana::SpillMultiVar(SPINEVAR_TT(vars::muon2024::opening_angle, SIGCUT))});
    vars_signal.insert({"true_dpT", ana::SpillMultiVar(SPINEVAR_TT(vars::interaction_pt, SIGCUT))});
    vars_signal.insert({"true_dphiT", ana::SpillMultiVar(SPINEVAR_TT(vars::phiT, SIGCUT))});
    vars_signal.insert({"true_edalphaT", ana::SpillMultiVar(SPINEVAR_TT(vars::alphaT, SIGCUT))});
    vars_signal.insert({"true_vertex_x", ana::SpillMultiVar(SPINEVAR_TT(vars::vertex_x, SIGCUT))});
    vars_signal.insert({"true_vertex_y", ana::SpillMultiVar(SPINEVAR_TT(vars::vertex_y, SIGCUT))});
    vars_signal.insert({"true_vertex_z", ana::SpillMultiVar(SPINEVAR_TT(vars::vertex_z, SIGCUT))});
    
    analysis.AddTree("signalNu", vars_signal, true);

    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();
}
