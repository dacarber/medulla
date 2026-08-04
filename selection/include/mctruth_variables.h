/**
 * @file mctruth_variables.h
 * @brief Definitions of analysis variables which can extract information from
 * the SRTrueInteraction object.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from the SRTrueInteraction object. Each variable
 * is implemented as a function which takes an SRTrueInteraction object as an
 * argument and returns a double. The association of an SRInteractionTruthDLP
 * object to an SRTrueInteraction object is handled upstream in the SpineVar
 * functions.
 * @author mueller@fnal.gov
 * @author rvizarr@fnal.gov
 */
#ifndef MCTRUTH_VARIABLES_H
#define MCTRUTH_VARIABLES_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"
#include "sbnanaobj/StandardRecord/SRVector3D.h"

#include "framework.h"
#include "variables.h"

/**
 * @namespace mctruth
 * @brief Namespace for organizing variables which act on true interactions.
 * @details This namespace is intended to be used for organizing variables
 * which act on true interactions. Each variable is implemented as a function
 * which takes an SRTrueInteraction object as an argument and returns a double.
 */
namespace mctruth
{
    /**
     * @brief Variable for the true neutrino energy.
     * @details This variable is intended to provide the true energy of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino energy.
     */
    template<typename T>
        double neutrino_energy(const T & obj) { return obj.E; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, neutrino_energy, neutrino_energy);

   /**
     * @brief Variable for the true interaction energy transfer.
     * @details This variable is intended to provide the true energy
     * transfer from the neutrino to the hadronic system. This is
     * defined in the lab frame.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true energy transfer into the hadronic system
     * in the lab frame.
     */
    template<typename T>
        double energy_transfer(const T & obj) { return obj.q0_lab; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, energy_transfer, energy_transfer);
  
    /**
     * @brief Variable for the true neutrino baseline.
     * @details This variable is intended to provide the true baseline of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino baseline.
     */
    template<typename T>
        double baseline(const T & obj) { return obj.baseline; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, baseline, baseline);

    /**
     * @brief Variable for the true neutrino PDG code.
     * @details This variable is intended to provide the true PDG code of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino PDG code.
     */
    template<typename T>
        double pdg(const T & obj) { return obj.pdg; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pdg, pdg);

    /**
     * @brief Variable for the PDG code of the parent of the neutrino.
     * @details This variable is intended to provide the PDG code of the
     * parent of the neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the PDG code of the parent of the neutrino.
     */
    template<typename T>
        double parent_pdg(const T & obj) { return obj.parent_pdg; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, parent_pdg, parent_pdg);

    /**
     * @brief Variable for the true neutrino current value.
     * @details This variable is intended to provide the true current value of
     * the parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino current value.
     */
    template<typename T>
        double cc(const T & obj) { return obj.iscc; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, cc, cc);

    /**
     * @brief Variable for the interaction mode of the interaction.
     * @details This variable is intended to provide the interaction mode of the
     * interaction. This is based on the GENIE interaction mode enumeration 
     * defined in the LArSoft MCNeutrino class.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the interaction mode.
     */
    template<typename T>
        double interaction_mode(const T & obj) { return obj.genie_mode; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, interaction_mode, interaction_mode);

    /**
     * @brief Variable for the interaction type of the interaction.
     * @details This variable is intended to provide the interaction type of the
     * interaction. This is based on the GENIE interaction type enumeration 
     * defined in the LArSoft MCNeutrino class.
     * @param T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the interaction type.
     */
    template<typename T>
        double interaction_type(const T & obj) { return obj.genie_inttype; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, interaction_type, interaction_type);

    /**
     * @brief Variable for the true off-axis angle of the neutrino.
     * @details This variable is intended to provide the true off-axis angle of
     * the parent neutrino that produced the interaction. The off-axis angle is
     * calculated as the angle between the neutrino momentum vector and the
     * beam axis (defined as the z-axis in both SBND and ICARUS).
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true off-axis angle of the neutrino in degrees.
     */
    template<typename T>
    double off_axis_angle(const T & obj)
    {
        const auto & neutrino_momentum = obj.momentum;
        double mag = std::sqrt(
            neutrino_momentum.x * neutrino_momentum.x +
            neutrino_momentum.y * neutrino_momentum.y +
            neutrino_momentum.z * neutrino_momentum.z
        );
        return 180./3.141592653589793 * std::acos(
            neutrino_momentum.z / mag
        );
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, off_axis_angle, off_axis_angle);

    /**
     * @brief Variable for the true neutrino initial energy (alias for
     * neutrino_energy).
     * @details This variable provides the same information as @ref
     * neutrino_energy but under the name "energy_init" to match the branch
     * naming convention used in trees that also record particle-level
     * "energy_init" values. Use @c type = "mctruth" in the TOML branch list.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino energy (GeV in CAF units; MeV in SPINE units).
     */
    template<typename T>
    double energy_init(const T & obj) { return obj.E; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, energy_init, energy_init);

    // -------------------------------------------------------------------------
    // Interaction vertex / GENIE-native kinematics
    // -------------------------------------------------------------------------

    /**
     * @brief Variable for the x-coordinate of the true interaction vertex.
     * @details Uses obj.position, matching the coordinate already relied on
     * by mctruth::fiducial_cut and mctruth::vertex_z_veto_cut in
     * mctruth_cuts.h.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the x-coordinate of the true interaction vertex.
     */
    template<typename T>
    double vertex_x(const T & obj) { return obj.position.x; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, vertex_x, vertex_x);

    /**
     * @brief Variable for the y-coordinate of the true interaction vertex.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the y-coordinate of the true interaction vertex.
     */
    template<typename T>
    double vertex_y(const T & obj) { return obj.position.y; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, vertex_y, vertex_y);

    /**
     * @brief Variable for the z-coordinate of the true interaction vertex.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the z-coordinate of the true interaction vertex.
     */
    template<typename T>
    double vertex_z(const T & obj) { return obj.position.z; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, vertex_z, vertex_z);

    /**
     * @brief Variable for the GENIE-native hadronic invariant mass W.
     * @details Read directly from the GENIE record (obj.w) rather than
     * reconstructed from final-state kinematics, unlike vars::W.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true hadronic invariant mass W.
     */
    template<typename T>
    double W(const T & obj) { return obj.w; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, W, W);

    /**
     * @brief Variable for the GENIE-native four-momentum transfer squared.
     * @details Read directly from the GENIE record (obj.Q2) rather than
     * reconstructed from final-state kinematics, unlike vars::Qsquared.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true four-momentum transfer squared Q^2.
     */
    template<typename T>
    double Qsquared(const T & obj) { return obj.Q2; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, Qsquared, Qsquared);

    // -------------------------------------------------------------------------
    // Shared helpers for final-state-primary kinematics/multiplicities
    // -------------------------------------------------------------------------

    /**
     * @brief Map a true PDG code to the pvars::Particle_t enumeration.
     * @details Used so that the GENIE final-state primary list (obj.prim)
     * can be classified with the same species ordering used by the
     * reco/true pvars::pid() function, and so the same
     * pcuts::final_state_signal_thresholds vector (populated from the
     * "general.fsthresh" TOML field) can be reused unmodified.
     * @param pdg the (possibly signed) true PDG code of the primary.
     * @return the pvars::Particle_t value, or pvars::kUnknown if the PDG
     * code is not photon/electron/muon/charged-pion/proton.
     */
    inline int genie_pid(int pdg)
    {
        switch(std::abs(pdg))
        {
            case 22:   return pvars::kPhoton;
            case 11:   return pvars::kElectron;
            case 13:   return pvars::kMuon;
            case 211:  return pvars::kPion;
            case 2212: return pvars::kProton;
            default:   return pvars::kUnknown;
        }
    }

    /**
     * @brief Rest mass corresponding to a pvars::Particle_t value.
     * @param pid the pvars::Particle_t value (as returned by genie_pid).
     * @return the rest mass in GeV, or 0 for photons/unrecognized species.
     */
    inline double genie_pid_mass(int pid)
    {
        switch(pid)
        {
            case pvars::kElectron: return ELECTRON_MASS;
            case pvars::kMuon:     return MUON_MASS;
            case pvars::kPion:     return PION_MASS;
            case pvars::kProton:   return PROTON_MASS;
            default:               return 0.0; // photon, or unrecognized species.
        }
    }

    /**
     * @brief Check if a true GENIE final-state primary meets the "final
     * state signal" requirement used by the kinematic/multiplicity
     * variables below.
     * @details Mirrors pcuts::final_state_signal's reco/true definition
     * (primary_classification(p) && ke(p) > threshold): requires the
     * primary to be a direct product of the GENIE interaction
     * (p.start_process == 0 — the same field mctruth_cuts.h's
     * no_neutral_pions/at_least_one_pi0 use to mean "primary"), to be one
     * of the five tracked species, and to have kinetic energy above the
     * same pcuts::final_state_signal_thresholds vector used at the
     * reco/true level. This is a distinct, stricter definition than the
     * existing per-species cuts in mctruth_cuts.h (single_muon,
     * no_photons, etc.), which do not filter on start_process.
     * @tparam P the type of a single entry of obj.prim.
     * @param p the GENIE final-state primary to check.
     * @return true if the primary meets the final state signal requirement.
     */
    template<typename P>
    bool prim_final_state_signal(const P & p)
    {
        if(p.start_process != 0) return false;
        int pid = genie_pid(p.pdg);
        if(pid == pvars::kUnknown) return false;
        double ke = p.genE - genie_pid_mass(pid);
        return ke > pcuts::final_state_signal_thresholds[pid];
    }

    // -------------------------------------------------------------------------
    // Kinematic variables built from obj.prim
    // -------------------------------------------------------------------------

    /**
     * @brief Variable for the transverse momentum of the interaction
     * counting only true final-state signal primaries.
     * @details Generator-level analog of vars::dpT: sums the transverse
     * momentum of every prim_final_state_signal() primary in obj.prim,
     * using each primary's own start position (p.start) for the local
     * beam-direction estimate — matching vars::dpT, which uses each
     * particle's own start point rather than the shared interaction
     * vertex — and generation-level momentum (p.genp), which pairs with
     * the p.genE used elsewhere in this file.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true transverse momentum of the final-state primaries.
     */
    template<typename T>
    double dpT(const T & obj)
    {
        utilities::three_vector pt = {0, 0, 0};
        for(const auto & p : obj.prim)
        {
            if(prim_final_state_signal(p))
            {
                utilities::three_vector momentum = {p.genp.x, p.genp.y, p.genp.z};
                utilities::three_vector vtx = {p.start.x, p.start.y, p.start.z};
                utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
                pt = utilities::add(pt, this_pt);
            }
        }
        return utilities::magnitude(pt);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, dpT, dpT);

    /**
     * @brief Variable for the transverse momentum of the interaction
     * counting only the leading true charged lepton and proton.
     * @details Generator-level analog of vars::dpT_lp.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true transverse momentum of the leading charged lepton
     * and proton, or PLACEHOLDERVALUE if either is absent.
     */
    template<typename T>
    double dpT_lp(const T & obj)
    {
        utilities::three_vector l_pt = {0, 0, 0};
        utilities::three_vector p_pt = {0, 0, 0};
        double l_ke(0), p_ke(0);
        for(const auto & p : obj.prim)
        {
            if(!prim_final_state_signal(p)) continue;
            const int pid = genie_pid(p.pdg);
            const double ke = p.genE - genie_pid_mass(pid);
            utilities::three_vector momentum = {p.genp.x, p.genp.y, p.genp.z};
            utilities::three_vector vtx = {p.start.x, p.start.y, p.start.z};
            if((pid == pvars::kElectron || pid == pvars::kMuon) && ke > l_ke)
            {
                l_ke = ke;
                l_pt = utilities::transverse_momentum(momentum, vtx);
            }
            else if(pid == pvars::kProton && ke > p_ke)
            {
                p_ke = ke;
                p_pt = utilities::transverse_momentum(momentum, vtx);
            }
        }
        if(l_ke == 0 || p_ke == 0)
            return PLACEHOLDERVALUE;
        return utilities::magnitude(utilities::add(l_pt, p_pt));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, dpT_lp, dpT_lp);

    /**
     * @brief Variable for dphi_T of the interaction at the GENIE generator
     * level.
     * @details Generator-level analog of vars::dphiT: transverse angle
     * between the leading charged lepton and the (summed) hadronic system,
     * built from true final-state signal primaries.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true phi_T of the interaction.
     */
    template<typename T>
    double dphiT(const T & obj)
    {
        utilities::three_vector lepton_pt = {0, 0, 0};
        utilities::three_vector hadronic_pt = {0, 0, 0};
        for(const auto & p : obj.prim)
        {
            if(!prim_final_state_signal(p)) continue;
            const int pid = genie_pid(p.pdg);
            utilities::three_vector momentum = {p.genp.x, p.genp.y, p.genp.z};
            utilities::three_vector vtx = {p.start.x, p.start.y, p.start.z};
            utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
            if(pid == pvars::kElectron || pid == pvars::kMuon)
                lepton_pt = this_pt;
            else if(pid > pvars::kMuon) // pion or proton -> hadronic system.
                hadronic_pt = utilities::add(hadronic_pt, this_pt);
        }
        return std::acos(-1 * utilities::dot_product(lepton_pt, hadronic_pt) / (utilities::magnitude(lepton_pt) * utilities::magnitude(hadronic_pt)));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, dphiT, dphiT);

    /**
     * @brief Variable for dalpha_T of the interaction at the GENIE
     * generator level.
     * @details Generator-level analog of vars::dalphaT: transverse angle
     * between the total final-state transverse momentum and the leading
     * charged lepton's transverse momentum.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true alpha_T of the interaction.
     */
    template<typename T>
    double dalphaT(const T & obj)
    {
        utilities::three_vector lepton_pt = {0, 0, 0};
        utilities::three_vector total_pt = {0, 0, 0};
        for(const auto & p : obj.prim)
        {
            if(!prim_final_state_signal(p)) continue;
            const int pid = genie_pid(p.pdg);
            utilities::three_vector momentum = {p.genp.x, p.genp.y, p.genp.z};
            utilities::three_vector vtx = {p.start.x, p.start.y, p.start.z};
            utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
            if(pid == pvars::kElectron || pid == pvars::kMuon)
                lepton_pt = this_pt;
            total_pt = utilities::add(total_pt, this_pt);
        }
        return std::acos(-1 * utilities::dot_product(total_pt, lepton_pt) / (utilities::magnitude(total_pt) * utilities::magnitude(lepton_pt)));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, dalphaT, dalphaT);

    /**
     * @brief Variable for the opening angle between the leading true
     * electron and proton.
     * @details Generator-level analog of vars::opening_angle. The reco/true
     * version uses each particle's start_dir; there is no equivalent field
     * on a GENIE primary, so the normalized generation-level momentum
     * (p.genp) is used as the direction instead.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true opening angle between the leading electron and
     * proton, or kNoMatchValue if either is absent.
     */
    template<typename T>
    double opening_angle(const T & obj)
    {
        double e_ke(-1), p_ke(-1);
        bool has_e(false), has_p(false);
        utilities::three_vector e_dir = {0, 0, 0}, p_dir = {0, 0, 0};
        for(const auto & p : obj.prim)
        {
            if(!prim_final_state_signal(p)) continue;
            const int pid = genie_pid(p.pdg);
            const double ke = p.genE - genie_pid_mass(pid);
            utilities::three_vector momentum = {p.genp.x, p.genp.y, p.genp.z};
            if(pid == pvars::kElectron && ke > e_ke)
            {
                e_ke = ke;
                e_dir = utilities::normalize(momentum);
                has_e = true;
            }
            else if(pid == pvars::kProton && ke > p_ke)
            {
                p_ke = ke;
                p_dir = utilities::normalize(momentum);
                has_p = true;
            }
        }
        if(!has_e || !has_p) return kNoMatchValue;
        return std::acos(utilities::dot_product(e_dir, p_dir));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, opening_angle, opening_angle);

    /**
     * @brief Variable for the total true visible energy of the interaction.
     * @details Generator-level analog of vars::visible_energy: sums
     * p.genE (total energy) over true final-state signal primaries, with
     * the same proton binding-energy correction applied by vars::visible_energy
     * (subtracting mass - PROTON_BINDING_ENERGY for each proton, leaving
     * KE + binding energy for protons but full E for other species,
     * mirroring the reco/true definition exactly).
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true visible energy of the interaction.
     */
    template<typename T>
    double visible_energy(const T & obj)
    {
        double energy(0);
        for(const auto & p : obj.prim)
        {
            if(!prim_final_state_signal(p)) continue;
            const int pid = genie_pid(p.pdg);
            energy += p.genE;
            if(pid == pvars::kProton)
                energy -= (genie_pid_mass(pid) - PROTON_BINDING_ENERGY);
        }
        return energy;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, visible_energy, visible_energy);

    // -------------------------------------------------------------------------
    // Multiplicities built from obj.prim
    // -------------------------------------------------------------------------

    /**
     * @brief Count true final-state primaries of a given species.
     * @details Shared helper for the primary/secondary multiplicity
     * variables below. Mirrors variables.h's detail::count_multiplicities:
     * a single KE threshold (GeV) is applied uniformly, and primary/
     * secondary status is taken from p.start_process (0 == primary).
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to count primaries in.
     * @param pid the pvars::Particle_t species to count, or a negative
     * value to count every tracked species (used by the "total secondary"
     * multiplicity).
     * @param want_primary true to count primary (start_process == 0)
     * particles, false to count non-primary (secondary) particles.
     * @param threshold the kinetic-energy threshold in GeV.
     * @return the count of matching true final-state primaries.
     */
    template<typename T>
    double count_prim_multiplicity(const T & obj, int pid, bool want_primary, double threshold)
    {
        double count(0);
        for(const auto & p : obj.prim)
        {
            const bool is_primary = (p.start_process == 0);
            if(is_primary != want_primary) continue;
            const int this_pid = genie_pid(p.pdg);
            if(pid >= 0 && this_pid != pid) continue;
            if(pid < 0 && this_pid == pvars::kUnknown) continue;
            const double ke = p.genE - genie_pid_mass(this_pid);
            if(ke >= threshold) ++count;
        }
        return count;
    }

    template<typename T>
    double photon_multiplicity(const T & obj, std::vector<double> params={0.025,})
    { return count_prim_multiplicity(obj, pvars::kPhoton, true, params[0]); }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, photon_multiplicity, photon_multiplicity);

    template<typename T>
    double electron_multiplicity(const T & obj, std::vector<double> params={0.025,})
    { return count_prim_multiplicity(obj, pvars::kElectron, true, params[0]); }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, electron_multiplicity, electron_multiplicity);

    template<typename T>
    double muon_multiplicity(const T & obj, std::vector<double> params={0.025,})
    { return count_prim_multiplicity(obj, pvars::kMuon, true, params[0]); }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, muon_multiplicity, muon_multiplicity);

    template<typename T>
    double pion_multiplicity(const T & obj, std::vector<double> params={0.025,})
    { return count_prim_multiplicity(obj, pvars::kPion, true, params[0]); }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_multiplicity, pion_multiplicity);

    template<typename T>
    double proton_multiplicity(const T & obj, std::vector<double> params={0.025,})
    { return count_prim_multiplicity(obj, pvars::kProton, true, params[0]); }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, proton_multiplicity, proton_multiplicity);

    template<typename T>
    double secondary_photon_multiplicity(const T & obj, std::vector<double> params={0.0})
    { return count_prim_multiplicity(obj, pvars::kPhoton, false, params[0]); }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, secondary_photon_multiplicity, secondary_photon_multiplicity);

    template<typename T>
    double secondary_electron_multiplicity(const T & obj, std::vector<double> params={0.0})
    { return count_prim_multiplicity(obj, pvars::kElectron, false, params[0]); }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, secondary_electron_multiplicity, secondary_electron_multiplicity);

    template<typename T>
    double secondary_muon_multiplicity(const T & obj, std::vector<double> params={0.0})
    { return count_prim_multiplicity(obj, pvars::kMuon, false, params[0]); }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, secondary_muon_multiplicity, secondary_muon_multiplicity);

    template<typename T>
    double secondary_pion_multiplicity(const T & obj, std::vector<double> params={0.0})
    { return count_prim_multiplicity(obj, pvars::kPion, false, params[0]); }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, secondary_pion_multiplicity, secondary_pion_multiplicity);

    template<typename T>
    double secondary_proton_multiplicity(const T & obj, std::vector<double> params={0.0})
    { return count_prim_multiplicity(obj, pvars::kProton, false, params[0]); }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, secondary_proton_multiplicity, secondary_proton_multiplicity);

    /**
     * @brief Variable for the total non-primary true particle multiplicity.
     * @details Generator-level analog of vars::secondary_multiplicity:
     * counts all non-primary (p.start_process != 0) prim entries of any of
     * the five tracked species with kinetic energy above threshold. Unlike
     * the reco/true version, this does not attempt to count untracked
     * species (e.g. neutrons, nuclear fragments) since they have no
     * well-defined mass/KE mapping in this file.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @param params [0] KE threshold in GeV (default 0 GeV).
     * @return total count of secondary tracked-species particles above
     * threshold.
     */
    template<typename T>
    double secondary_multiplicity(const T & obj, std::vector<double> params={0.0})
    { return count_prim_multiplicity(obj, -1, false, params[0]); }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, secondary_multiplicity, secondary_multiplicity);

} // namespace mctruth
#endif
