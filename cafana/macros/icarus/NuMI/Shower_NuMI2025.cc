
/**
 * @file ccpi0AnaMC.C
 * @brief The main analysis macro for the ICARUS numu CC pi0 selection..
 * @details This macro drives the analysis by configuring the variables, cuts,
 * and samples to be used in the analysis. This is accomplished through the use
 * of the Analysis class, which containerizes the configuration of the analysis
 * and reduces the amount of boilerplate code needed to run the analysis.
 * @author lkashur@colostate.edu
*/
#define PLACEHOLDERVALUE std::numeric_limits<double>::quiet_NaN()
#define PIDFUNC pvars::pid
#define PRIMARYFUNC pvars::primary_classification
//#define PIDFUNC pvars::custom_pid
#define PROTON_BINDING_ENERGY 30.9 // MeV
#define BEAM_IS_NUMI true
#define WRITE_PURITY_TREES false

#include "include/mctruth.h"
#include "include/variables.h"
#include "include/cuts.h"

#include "include/NuMI2025/variables_nue.h"
#include "include/NuMI2025/cuts_nue.h"


#include "include/spinevar.h"
#include "include/analysis.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "TDirectory.h"
#include "TFile.h"

int main()
{
    // Output filename
    //ana::Analysis analysis("sbnd_ccpi0_mc_03_april_2025");
    //ana::Analysis analysis("sbnd_ccpi0_mc_03_april_2025");
    //ana::Analysis analysis("icarus_bnb_ccpi0_mc_08_april_2025");
    ana::Analysis analysis("testtest");

    ana::SpectrumLoader mc("/pnfs/icarus/persistent/users/dcarber/spine/combined_files/NuMI_CV_flat_cafs/NuMI_CV_flat_caf_combo*.root");
    analysis.AddLoader("mc", &mc, true);

    // SBND
    //ana::SpectrumLoader sbnd("/pnfs/icarus/persistent/users/mueller/sbnd/updated/flat/larcv_sbnd_bnb_cosmics_spine_updated.flat.root");
    //analysis.AddLoader("sbnd", &sbnd, true);

    //ana::SpectrumLoader intime("/pnfs/icarus/persistent/users/mueller/sbnd/larcv_sbnd_intime_spine.flat.root");
    //analysis.AddLoader("intime", &intime, true);

    /**
     * @brief Add variabls for selected interactions (in-phase) to the analysis.
     * @details This adds a set of variables to the analysis by creating a map
     * of variable names and SpillMultiVars that provide the functionality to
     * create the variables.  These names are used in the TTree that is created
     * by the Tree class to store the results of the analysis.
     */
    #define CUT cuts::has_shower
    #define TCUT cuts::no_cut
    std::map<std::string, ana::SpillMultiVar> vars_reco_shower;
    //vars_selected_nu_phase.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    //vars_selected_nu_phase.insert({"nu_E", SpineVar<TTYPE,RTYPE>(&vars::neutrino_E, &CUT, &TCUT)});
    //vars_selected_nu_phase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::nue::cut_type, &CUT, &TCUT)}); // GUNDAM
    //vars_selected_nu_phase.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::nue::is_signal_mc, &CUT, &TCUT)}); // GUNDAM
    //vars_selected_nu_phase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::nue::is_not_data, &CUT, &TCUT)}); // GUNDAM
    //vars_selected_nu_phase.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    //vars_selected_nu_phase.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::nue::category, &CUT, &TCUT)});
    vars_reco_shower.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::nue::category_topology, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_topology", SpineVar<RTYPE,RTYPE>(&vars::nue::reco_topology, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_topology", SpineVar<TTYPE,RTYPE>(&vars::nue::truth_topology, &CUT, &TCUT)});
    //vars_selected_nu_phase.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_electron_energy", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_ke, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_electron_energy", SpineVar<TTYPE,RTYPE>(&vars::nue::leading_electron_ke, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_electron_pT_mag", SpineVar<RTYPE,RTYPE>(&vars::nue::electron_transverse_momentum_mag, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_electron_pT_mag", SpineVar<TTYPE,RTYPE>(&vars::nue::electron_transverse_momentum_mag, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_electron_NuMI_angle", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_NuMI_angle, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_electron_NuMI_angle", SpineVar<TTYPE,RTYPE>(&vars::nue::leading_electron_NuMI_angle, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_electron_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_vertex_distance, &CUT, &TCUT)});   
    vars_reco_shower.insert({"reco_electron_axial_spread", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_axial_spread, &CUT, &TCUT)});   
    vars_reco_shower.insert({"reco_electron_dir_spread", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_directional_spread, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_electron_softmax", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_softmax, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_electron_primary_score", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_primary_score, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_electron_dedx", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_dedx, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_NuMI_polar", SpineVar<RTYPE,RTYPE>(&vars::nue::NuMI_polar_angle, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_NuMI_polar", SpineVar<TTYPE,RTYPE>(&vars::nue::NuMI_polar_angle, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_NuMI_azi", SpineVar<RTYPE,RTYPE>(&vars::nue::NuMI_azimuthal_angle, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_NuMI_azi", SpineVar<TTYPE,RTYPE>(&vars::nue::NuMI_azimuthal_angle, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_proton_energy", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_proton_ke, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_proton_energy", SpineVar<TTYPE,RTYPE>(&vars::nue::leading_proton_ke, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_proton_pT_mag", SpineVar<RTYPE,RTYPE>(&vars::nue::proton_transverse_momentum_mag, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_proton_pT_mag", SpineVar<TTYPE,RTYPE>(&vars::nue::proton_transverse_momentum_mag, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_proton_NuMI_angle", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_proton_NuMI_angle, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_proton_NuMI_angle", SpineVar<TTYPE,RTYPE>(&vars::nue::leading_proton_NuMI_angle, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_proton_softmax", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_proton_softmax, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_proton_muon_softmax", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_proton_muon_softmax, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_proton_pion_softmax", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_proton_pion_softmax, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_proton_primary_score", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_proton_primary_score, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_dpT", SpineVar<RTYPE,RTYPE>(&vars::nue::delta_pT, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_dpT", SpineVar<TTYPE,RTYPE>(&vars::nue::delta_pT, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_dalphaT", SpineVar<RTYPE,RTYPE>(&vars::nue::delta_alphaT, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_dalphaT", SpineVar<TTYPE,RTYPE>(&vars::nue::delta_alphaT, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_dphiT", SpineVar<RTYPE,RTYPE>(&vars::nue::delta_phiT, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_dphiT", SpineVar<TTYPE,RTYPE>(&vars::nue::delta_phiT, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_opening_angle", SpineVar<RTYPE,RTYPE>(&vars::nue::opening_angle, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_opening_angle", SpineVar<TTYPE,RTYPE>(&vars::nue::opening_angle, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_visible_energy", SpineVar<RTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_visible_energy", SpineVar<TTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_vertex_x", SpineVar<TTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_vertex_y", SpineVar<TTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_reco_shower.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    vars_reco_shower.insert({"true_vertex_z", SpineVar<TTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});

    analysis.AddTree("SelectedReco_shower", vars_reco_shower, false);

    #undef TCUT
    #define TCUT cuts::has_shower
    std::map<std::string, ana::SpillMultiVar> vars_true_shower;
    vars_true_shower.insert({"nu_id", SpineVar<TTYPE,RTYPE>(&vars::neutrino_id, &CUT, &TCUT)});
    vars_true_shower.insert({"nu_E", SpineVar<TTYPE,RTYPE>(&vars::neutrino_E, &CUT, &TCUT)});
    //vars_selected_cosmic_phase.insert({"CutType", SpineVar<RTYPE,RTYPE>(&vars::nue::cut_type, &CUT, &TCUT)}); // GUNDAM
    //vars_selected_cosmic_phase.insert({"IsSignal", SpineVar<TTYPE,RTYPE>(&vars::ccpi0ana_phase::is_signal_mc, &CUT, &TCUT)}); // GUNDAM
    //vars_selected_cosmic_phase.insert({"IsData", SpineVar<RTYPE,RTYPE>(&vars::ccpi0ana_phase::is_not_data, &CUT, &TCUT)}); // GUNDAM
    vars_true_shower.insert({"baseline", SpineVar<MCTRUTH,RTYPE>(&mctruth::true_neutrino_baseline, &CUT, &TCUT)});
    vars_true_shower.insert({"category", SpineVar<TTYPE,RTYPE>(&vars::nue::category, &CUT, &TCUT)});
    vars_true_shower.insert({"category_topology", SpineVar<TTYPE,RTYPE>(&vars::nue::category_topology, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_topology", SpineVar<RTYPE,RTYPE>(&vars::nue::reco_topology, &CUT, &TCUT)});
    vars_true_shower.insert({"true_topology", SpineVar<TTYPE,RTYPE>(&vars::nue::truth_topology, &CUT, &TCUT)});
    vars_true_shower.insert({"interaction_mode", SpineVar<MCTRUTH,RTYPE>(&mctruth::interaction_mode, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_electron_energy", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_ke, &CUT, &TCUT)});
    vars_true_shower.insert({"true_electron_energy", SpineVar<TTYPE,RTYPE>(&vars::nue::leading_electron_ke, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_electron_pT_mag", SpineVar<RTYPE,RTYPE>(&vars::nue::electron_transverse_momentum_mag, &CUT, &TCUT)});
    vars_true_shower.insert({"true_electron_pT_mag", SpineVar<TTYPE,RTYPE>(&vars::nue::electron_transverse_momentum_mag, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_electron_NuMI_angle", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_NuMI_angle, &CUT, &TCUT)});
    vars_true_shower.insert({"true_electron_NuMI_angle", SpineVar<TTYPE,RTYPE>(&vars::nue::leading_electron_NuMI_angle, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_electron_conv_dist", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_vertex_distance, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_electron_axial_spread", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_axial_spread, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_electron_dir_spread", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_directional_spread, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_electron_softmax", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_softmax, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_electron_primary_score", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_primary_score, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_electron_dedx", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_electron_dedx, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_NuMI_polar", SpineVar<RTYPE,RTYPE>(&vars::nue::NuMI_polar_angle, &CUT, &TCUT)});
    vars_true_shower.insert({"true_NuMI_polar", SpineVar<TTYPE,RTYPE>(&vars::nue::NuMI_polar_angle, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_NuMI_azi", SpineVar<RTYPE,RTYPE>(&vars::nue::NuMI_azimuthal_angle, &CUT, &TCUT)});
    vars_true_shower.insert({"true_NuMI_azi", SpineVar<TTYPE,RTYPE>(&vars::nue::NuMI_azimuthal_angle, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_proton_energy", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_proton_ke, &CUT, &TCUT)});
    vars_true_shower.insert({"true_proton_energy", SpineVar<TTYPE,RTYPE>(&vars::nue::leading_proton_ke, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_proton_pT_mag", SpineVar<RTYPE,RTYPE>(&vars::nue::proton_transverse_momentum_mag, &CUT, &TCUT)});
    vars_true_shower.insert({"true_proton_pT_mag", SpineVar<TTYPE,RTYPE>(&vars::nue::proton_transverse_momentum_mag, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_proton_NuMI_angle", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_proton_NuMI_angle, &CUT, &TCUT)});
    vars_true_shower.insert({"true_proton_NuMI_angle", SpineVar<TTYPE,RTYPE>(&vars::nue::leading_proton_NuMI_angle, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_proton_softmax", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_proton_softmax, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_proton_muon_softmax", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_proton_muon_softmax, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_proton_pion_softmax", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_proton_pion_softmax, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_proton_primary_score", SpineVar<RTYPE,RTYPE>(&vars::nue::leading_proton_primary_score, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_dpT", SpineVar<RTYPE,RTYPE>(&vars::nue::delta_pT, &CUT, &TCUT)});
    vars_true_shower.insert({"true_dpT", SpineVar<TTYPE,RTYPE>(&vars::nue::delta_pT, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_dalphaT", SpineVar<RTYPE,RTYPE>(&vars::nue::delta_alphaT, &CUT, &TCUT)});
    vars_true_shower.insert({"true_dalphaT", SpineVar<TTYPE,RTYPE>(&vars::nue::delta_alphaT, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_dphiT", SpineVar<RTYPE,RTYPE>(&vars::nue::delta_phiT, &CUT, &TCUT)});
    vars_true_shower.insert({"true_dphiT", SpineVar<TTYPE,RTYPE>(&vars::nue::delta_phiT, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_opening_angle", SpineVar<RTYPE,RTYPE>(&vars::nue::opening_angle, &CUT, &TCUT)});
    vars_true_shower.insert({"true_opening_angle", SpineVar<TTYPE,RTYPE>(&vars::nue::opening_angle, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_visible_energy", SpineVar<RTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_true_shower.insert({"true_visible_energy", SpineVar<TTYPE,RTYPE>(&vars::visible_energy, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_vertex_x", SpineVar<RTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_true_shower.insert({"true_vertex_x", SpineVar<TTYPE,RTYPE>(&vars::vertex_x, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_vertex_y", SpineVar<RTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_true_shower.insert({"true_vertex_y", SpineVar<TTYPE,RTYPE>(&vars::vertex_y, &CUT, &TCUT)});
    vars_true_shower.insert({"reco_vertex_z", SpineVar<RTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    vars_true_shower.insert({"true_vertex_z", SpineVar<TTYPE,RTYPE>(&vars::vertex_z, &CUT, &TCUT)});
    analysis.AddTree("SelectedTrue_shower", vars_true_shower, true);

    
    
    /**
     * @brief Run the analysis.
     * @details This runs the analysis on the samples specified by the
     * SpectrumLoaders and variables added to the Analysis class. It loops over
     * each sample (here only one), applies the cuts and variables to the data,
     * and stores the results in a TFile.
     */
    analysis.Go();
}
