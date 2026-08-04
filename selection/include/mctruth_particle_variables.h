/**
 * @file mctruth_particle_variables.h
 * @brief Definitions of selectors and analysis variables which act on a
 * single true final-state primary (an entry of SRTrueInteraction::prim).
 * @details This file is the MCTruth-particle-level counterpart to
 * selectors.h/particle_variables.h: selectors.h::leading_electron etc.
 * select an index into an interaction's SPINE particle list
 * (TType/RType::particles), while mctruth_selectors::leading_electron etc.
 * select an index into a GENIE interaction's true final-state primary list
 * (MCTruth::prim). Similarly, particle_variables.h's pvars namespace
 * extracts a value from a single SPINE particle, while this file's
 * mctruth_pvars namespace extracts a value from a single GENIE primary.
 * @author mueller@fnal.gov
 * @author rvizarr@fnal.gov
 */
#ifndef MCTRUTH_PARTICLE_VARIABLES_H
#define MCTRUTH_PARTICLE_VARIABLES_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

#include "framework.h"
#include "mctruth_variables.h"

/**
 * @namespace mctruth_selectors
 * @brief Namespace for organizing selectors which act on the GENIE-truth
 * MCTruth object, selecting an index into MCTruth::prim.
 * @details Mirrors selectors.h's leading_electron/leading_proton/
 * leading_secondary_electron/leading_secondary_proton, but restricted to
 * the four selectors actually used with type = "both_particle" in the
 * existing TOML configurations. Unlike selectors::leading_particle_index,
 * which applies no energy threshold, these also apply no threshold — only
 * a PDG match (and, for the "secondary" variants, a p.start_process != 0
 * filter) — for a faithful, unmodified port of the reco/true semantics.
 */
namespace mctruth_selectors
{
    /**
     * @brief Finds the index of the highest-KE MCTruth::prim entry of a
     * given species, optionally restricted to non-primary particles.
     * @details Shared helper for the four selectors below. Applies no
     * energy threshold, matching selectors::leading_particle_index /
     * leading_secondary_particle_index exactly.
     * @param obj the SRTrueInteraction to operate on.
     * @param pid the pvars::Particle_t species to select.
     * @param secondary_only if true, only consider primaries with
     * p.start_process != 0 (non-primary, mirroring
     * !pvars::primary_classification(p) at the reco/true level).
     * @return the index of the matching highest-KE primary, or kNoMatch.
     */
    inline size_t leading_prim_index(const MCTruth & obj, int pid, bool secondary_only)
    {
        double leading_ke(0);
        size_t index(kNoMatch);
        for(size_t i = 0; i < obj.prim.size(); ++i)
        {
            const auto & p = obj.prim[i];
            if(mctruth::genie_pid(p.pdg) != pid) continue;
            if(secondary_only && p.start_process == 0) continue;
            const double ke = p.genE - mctruth::genie_pid_mass(pid);
            if(ke > leading_ke)
            {
                leading_ke = ke;
                index = i;
            }
        }
        return index;
    }

    template<typename T>
    size_t leading_electron(const T & obj) { return leading_prim_index(obj, pvars::kElectron, false); }
    REGISTER_MCTRUTH_SELECTOR(leading_electron, leading_electron);

    template<typename T>
    size_t leading_proton(const T & obj) { return leading_prim_index(obj, pvars::kProton, false); }
    REGISTER_MCTRUTH_SELECTOR(leading_proton, leading_proton);

    template<typename T>
    size_t leading_secondary_electron(const T & obj) { return leading_prim_index(obj, pvars::kElectron, true); }
    REGISTER_MCTRUTH_SELECTOR(leading_secondary_electron, leading_secondary_electron);

    template<typename T>
    size_t leading_secondary_proton(const T & obj) { return leading_prim_index(obj, pvars::kProton, true); }
    REGISTER_MCTRUTH_SELECTOR(leading_secondary_proton, leading_secondary_proton);
}

/**
 * @namespace mctruth_pvars
 * @brief Namespace for organizing variables which act on a single true
 * final-state primary (an entry of MCTruth::prim).
 * @details Mirrors particle_variables.h's pvars namespace, ported to the
 * GENIE-truth primary type. Only the subset of pvars functions actually
 * used with type = "both_particle" in the existing TOML configurations are
 * provided; SPINE-reconstruction-specific quantities (softmax PID scores,
 * MCS/CSDA energy branching, cluster-derived spread/straightness, etc.)
 * have no generator-level analog and are not ported.
 */
namespace mctruth_pvars
{
    /**
     * @brief Variable for the x-coordinate of the primary's starting point.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the x-coordinate of the primary's starting point.
     */
    template<typename T>
    double start_x(const T & p) { return p.start.x; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, start_x, start_x);

    /**
     * @brief Variable for the y-coordinate of the primary's starting point.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the y-coordinate of the primary's starting point.
     */
    template<typename T>
    double start_y(const T & p) { return p.start.y; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, start_y, start_y);

    /**
     * @brief Variable for the z-coordinate of the primary's starting point.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the z-coordinate of the primary's starting point.
     */
    template<typename T>
    double start_z(const T & p) { return p.start.z; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, start_z, start_z);

    /**
     * @brief Variable for the x-coordinate of the primary's end point.
     * @details Guards against an infinite end point the same way
     * pvars::end_x does, though a GENIE primary's G4-tracked end point is
     * not expected to be infinite in practice.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the x-coordinate of the primary's end point.
     */
    template<typename T>
    double end_x(const T & p) { return std::isinf(p.end.x) ? PLACEHOLDERVALUE : (double)p.end.x; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, end_x, end_x);

    /**
     * @brief Variable for the y-coordinate of the primary's end point.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the y-coordinate of the primary's end point.
     */
    template<typename T>
    double end_y(const T & p) { return std::isinf(p.end.y) ? PLACEHOLDERVALUE : (double)p.end.y; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, end_y, end_y);

    /**
     * @brief Variable for the z-coordinate of the primary's end point.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the z-coordinate of the primary's end point.
     */
    template<typename T>
    double end_z(const T & p) { return std::isinf(p.end.z) ? PLACEHOLDERVALUE : (double)p.end.z; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, end_z, end_z);

    /**
     * @brief Variable for the true kinetic energy of the primary.
     * @details GENIE-native analog of pvars::ke, using p.genE and the
     * species rest mass (via mctruth::genie_pid/genie_pid_mass) rather than
     * SPINE's calorimetric/CSDA/MCS branching.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the true kinetic energy of the primary.
     */
    template<typename T>
    double ke(const T & p) { return p.genE - mctruth::genie_pid_mass(mctruth::genie_pid(p.pdg)); }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, ke, ke);

    /**
     * @brief Variable aliasing the true kinetic energy of the primary.
     * @details There is no generator-level concept of a separately
     * calorimetrically-reconstructed energy (pvars::calo_ke sums energy
     * deposited per SPINE spacepoint) — a GENIE primary carries only its
     * generator/G4-truth kinematics. This alias exists so branches
     * requesting "calo_ke" under type = "both_mctruth_particle" resolve to
     * the true kinetic energy rather than failing to register.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the true kinetic energy of the primary (identical to ke(p)).
     */
    template<typename T>
    double calo_ke(const T & p) { return ke(p); }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, calo_ke, calo_ke);

    /**
     * @brief Variable for the track length of the primary.
     * @details Read directly from the GENIE/G4-truth record (p.length),
     * unlike pvars::length, which reads a SPINE-reconstructed length.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the true track length of the primary.
     */
    template<typename T>
    double length(const T & p) { return p.length; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, length, length);

    /**
     * @brief Variable for the magnitude of the primary's momentum.
     * @details GENIE-native analog of pvars::p, computed from the
     * generation-level momentum vector (p.genp) rather than a
     * SPINE-reconstructed field.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the magnitude of the true momentum of the primary.
     */
    template<typename T>
    double p(const T & p)
    {
        utilities::three_vector momentum = {p.genp.x, p.genp.y, p.genp.z};
        return utilities::magnitude(momentum);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, p, p);

    /**
     * @brief Variable for the transverse momentum of the primary.
     * @details GENIE-native analog of pvars::dpT: transverse component of
     * the primary's momentum (p.genp) with respect to the assumed neutrino
     * direction, using the primary's own start position (p.start) for the
     * local beam-direction estimate — matching pvars::dpT, which uses each
     * particle's own start point.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the true transverse momentum of the primary.
     */
    template<typename T>
    double dpT(const T & p)
    {
        utilities::three_vector momentum = {p.genp.x, p.genp.y, p.genp.z};
        utilities::three_vector vtx = {p.start.x, p.start.y, p.start.z};
        utilities::three_vector pt = utilities::transverse_momentum(momentum, vtx);
        return utilities::magnitude(pt);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, dpT, dpT);

    /**
     * @brief Variable for the polar angle (w.r.t. the z-axis) of the
     * primary.
     * @details GENIE-native analog of pvars::polar_angle. There is no
     * MCTruth::prim analog of pvars::start_dir, so the direction is taken
     * as the normalized generation-level momentum (p.genp) instead.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the true polar angle of the primary.
     */
    template<typename T>
    double polar_angle(const T & p)
    {
        utilities::three_vector momentum = {p.genp.x, p.genp.y, p.genp.z};
        utilities::three_vector dir = utilities::normalize(momentum);
        return std::acos(std::get<2>(dir));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, polar_angle, polar_angle);

    /**
     * @brief Variable for the azimuthal angle of the primary.
     * @details GENIE-native analog of pvars::azimuthal_angle, using the
     * normalized generation-level momentum (p.genp) as the direction.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the true azimuthal angle of the primary.
     */
    template<typename T>
    double azimuthal_angle(const T & p)
    {
        utilities::three_vector momentum = {p.genp.x, p.genp.y, p.genp.z};
        utilities::three_vector dir = utilities::normalize(momentum);
        return std::atan2(std::get<1>(dir), std::get<0>(dir));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, azimuthal_angle, azimuthal_angle);

    /**
     * @brief Variable for the polar angle of the primary w.r.t. the NuMI
     * beam direction.
     * @details GENIE-native analog of pvars::NuMI_polar_angle. Uses the
     * primary's start position (p.start) and the normalized
     * generation-level momentum (p.genp) as the direction, in place of
     * pvars::start_point/start_dir.
     * @tparam T the type of the GENIE primary (an entry of MCTruth::prim).
     * @param p the primary to apply the variable on.
     * @return the true NuMI polar angle in radians (0 = along beam
     * direction).
     */
    template<typename T>
    double NuMI_polar_angle(const T & p)
    {
        utilities::three_vector momentum = {p.genp.x, p.genp.y, p.genp.z};
        utilities::three_vector dir = utilities::normalize(momentum);
        double bx = 31512.0380 + p.start.x;
        double by =  3364.4912 + p.start.y;
        double bz = 73363.2532 + p.start.z;
        double bmag = std::sqrt(bx*bx + by*by + bz*bz);
        double cos_theta = (std::get<0>(dir)*bx + std::get<1>(dir)*by + std::get<2>(dir)*bz) / bmag;
        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
        return std::acos(cos_theta);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruthParticle, NuMI_polar_angle, NuMI_polar_angle);
} // namespace mctruth_pvars
#endif
