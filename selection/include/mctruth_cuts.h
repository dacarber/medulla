/**
 * @file mctruth_cuts.h
 * @brief Definitions of analysis cuts applied at the GENIE generator level.
 * @details This file contains definitions of analysis cuts which operate
 * directly on the SRTrueInteraction object. Each cut is implemented as a
 * function which takes an SRTrueInteraction object as an argument and returns
 * a bool. The association of an SRInteractionTruthDLP object to an
 * SRTrueInteraction object is handled upstream in the framework.
 * @author mueller@fnal.gov
 * @author rvizarr@fnal.gov
 */
#ifndef MCTRUTH_CUTS_H
#define MCTRUTH_CUTS_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"
#include <cmath>
#include <limits>

#include "framework.h"
#include "variables.h"

/**
 * @namespace mctruth
 * @brief Namespace for organizing cuts which act on true interactions.
 * @details This namespace is intended to be used for organizing cuts
 * which act on true interactions. Each cut is implemented as a function
 * which takes an SRTrueInteraction object as an argument and returns a bool.
 */
namespace mctruth {
/**
 * @brief Cut for charged current interactions at the generator (GENIE) level.
 * @details Distinct from the SPINE truth-level iscc cut. Uses obj.iscc
 * directly from the MCTruth object.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @return true if the interaction is charged current.
 */
template <typename T> bool iscc(const T &obj) { return obj.iscc; }
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, iscc, iscc);

/**
 * @brief Cut for exactly one true final state muon above threshold.
 * @details Applied at the GENIE generator level using obj.prim. The
 * kinetic energy is computed from the GENIE genE field.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @param params KE threshold in MeV, defaults to 143.425 MeV.
 * @return true if exactly one muon above threshold.
 */
template <typename T>
bool single_muon(const T &obj, std::vector<double> params = {
                                   143.425,
                               }) {
  if (params.empty())
    params.push_back(143.425);
  double upper =
      params.size() > 1 ? params[1] : std::numeric_limits<double>::infinity();
  int num_muons(0);
  for (const auto &p : obj.prim) {
    if (abs(p.pdg) == 13) {
      double ke = (p.genE - (MUON_MASS));
      if (ke >= params[0] && ke <= upper)
        num_muons++;
    }
  }
  return num_muons == 1;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, single_muon, single_muon);

/**
 * @brief Cut for zero true final state charged pions above threshold.
 * @details Applied at the GENIE generator level using obj.prim. The
 * kinetic energy is computed from the GENIE genE field.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @param params KE threshold in MeV, defaults to 25 MeV.
 * @return true if no charged pions above threshold.
 */
template <typename T>
bool no_charged_pions(const T &obj, std::vector<double> params = {
                                        25.0,
                                    }) {
  if (params.empty())
    params.push_back(25.0);
  double upper =
      params.size() > 1 ? params[1] : std::numeric_limits<double>::infinity();
  for (const auto &p : obj.prim) {
    if (abs(p.pdg) == 211) {
      double ke = (p.genE - (PION_MASS));
      if (ke >= params[0] && ke <= upper)
        return false;
    }
  }
  return true;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_charged_pions,
                   no_charged_pions);

/**
 * @brief Cut for zero true final state neutral pions.
 * @details Applied at the GENIE generator level using obj.prim. Only
 * primary particles (start_process == 0) are considered.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @return true if no neutral pions are present.
 */
template <typename T> bool no_neutral_pions(const T &obj) {
  for (const auto &p : obj.prim) {
    if (p.start_process != 0)
      continue;
    if (p.pdg == 111)
      return false;
  }
  return true;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_neutral_pions,
                   no_neutral_pions);

/**
 * @brief Cut for zero true final state photons above threshold.
 * @details Applied at the GENIE generator level using obj.prim. The
 * energy is computed from the GENIE genE field.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @param params energy threshold in MeV, defaults to 25 MeV.
 * @return true if no photons above threshold.
 */
template <typename T>
bool no_photons(const T &obj, std::vector<double> params = {
                                  25.0,
                              }) {
  if (params.empty())
    params.push_back(25.0);
  double upper =
      params.size() > 1 ? params[1] : std::numeric_limits<double>::infinity();
  for (const auto &p : obj.prim) {
    if (p.pdg == 22) {
      double energy = p.genE;
      if (energy >= params[0] && energy <= upper)
        return false;
    }
  }
  return true;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_photons, no_photons);

/**
 * @brief Cut for zero true final state electrons above threshold.
 * @details Applied at the GENIE generator level using obj.prim. The
 * kinetic energy is computed from the GENIE genE field.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @param params KE threshold in MeV, defaults to 25 MeV.
 * @return true if no electrons above threshold.
 */
template <typename T>
bool no_electrons(const T &obj, std::vector<double> params = {
                                    25.0,
                                }) {
  if (params.empty())
    params.push_back(25.0);
  double upper =
      params.size() > 1 ? params[1] : std::numeric_limits<double>::infinity();
  for (const auto &p : obj.prim) {
    if (abs(p.pdg) == 11) {
      double ke = (p.genE - (ELECTRON_MASS));
      if (ke >= params[0] && ke <= upper)
        return false;
    }
  }
  return true;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_electrons, no_electrons);

/**
 * @brief Cut for exactly one true final state proton above threshold.
 * @details Applied at the GENIE generator level using obj.prim. The
 * kinetic energy is computed from the GENIE genE field.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @param params KE threshold in MeV, defaults to 50 MeV.
 * @return true if exactly one proton above threshold.
 */
template <typename T>
bool single_proton(const T &obj, std::vector<double> params = {
                                     50.0,
                                 }) {
  if (params.empty())
    params.push_back(50.0);
  double upper =
      params.size() > 1 ? params[1] : std::numeric_limits<double>::infinity();
  int count(0);
  for (const auto &p : obj.prim) {
    if (p.pdg == 2212) {
      double ke = (p.genE - (PROTON_MASS));
      if (ke >= params[0] && ke <= upper)
        count++;
    }
  }
  return count == 1;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, single_proton, single_proton);

/**
 * @brief Apply a cut for a maximum energy transfer
 * @details This function applies a cut to select interactions below
 * a maximum energy transfer (aka, omega or q0).This cut is on the
 * true energy transfer from the neutrino to the hadronic system.
 * This is defined in the lab frame.
 * @tparam T the type of interaction (true or reco).
 * @param obj the interaction to select on.
 * @param params the parameters for the cut. In this case, this is a
 * double representing an energy transfer in GeV
 * @return true if the energy transfer is below the given value
 */
template <class T>
bool is_below_max_e_transfer(const T &obj, std::vector<double> params = {
                                               0.1,
                                           }) {
  return obj.q0_lab < params[0];
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, is_below_max_e_transfer,
                   is_below_max_e_transfer);

/**
 * @brief Cut to select neutrino interactions at the GENIE level.
 * @details At the MCTruth scope the framework only invokes this cut when
 * the SPINE truth interaction has a valid GENIE match (nu_id >= 0), so
 * all objects reaching this cut are already neutrinos.  Returns true
 * unconditionally; exists to mirror the "true"-scope neutrino cut so
 * category definitions can use the same names in both scopes.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @return always true.
 */
template <typename T> bool neutrino(const T &obj) { return true; }
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, neutrino, neutrino);

/**
 * @brief Fiducial volume cut at the GENIE generator level.
 * @details Checks whether the neutrino interaction vertex (obj.position)
 * lies within the ICARUS fiducial volume.  The bounds mirror those used
 * by the SPINE truth-DLP fiducial cut, adapted for the SRTrueInteraction
 * position coordinate system:
 *   |x| in (10, 190) cm  — 10 cm from cathode, 190 cm from wire planes
 *   z in (10, 450) cm    — 10 cm from upstream/downstream ends
 *   y in (-190, 120) cm  — vertical extent of the active volume
 * An additional veto removes the cathode-corner region that is excluded
 * from the SPINE is_fiducial definition.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @return true if the vertex is in the fiducial volume.
 */
template <typename T> bool fiducial_cut(const T &obj) {
  const double x = obj.position.x;
  const double y = obj.position.y;
  const double z = obj.position.z;
  bool in_fv = (std::abs(x) > 10.0 && std::abs(x) < 190.0) &&
               (y > -190.0 && y < 120.0) && (z > 10.0 && z < 450.0);
  bool cathode_corner = (x > 210.215 && y > 60.0 && z > 290.0 && z < 390.0);
  return in_fv && !cathode_corner;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, fiducial_cut, fiducial_cut);

/**
 * @brief Veto true interactions whose vertex lies within a z-range.
 * @details Mirrors the reco/true-level vertex_z_veto_cut, but operates on
 * the SRTrueInteraction position coordinate (obj.position.z) rather than
 * obj.vertex[2].
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @param params the [low, high] bounds of the z-range to veto.
 * @return false if the vertex z is within [params[0], params[1]], true
 * otherwise.
 */
template <typename T>
bool vertex_z_veto_cut(const T &obj, std::vector<double> params = {}) {
  if (params.size() == 2 && obj.position.z >= params[0] &&
      obj.position.z <= params[1])
    return false;
  else
    return true;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, vertex_z_veto_cut,
                   vertex_z_veto_cut);

/**
 * @brief Cut for exactly one true final state electron above threshold.
 * @details Applied at the GENIE generator level using obj.prim. The
 * kinetic energy is computed from the GENIE genE field.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @param params KE threshold in MeV, defaults to 25 MeV.
 * @return true if exactly one electron above threshold.
 */
template <typename T>
bool single_electron(const T &obj, std::vector<double> params = {
                                       25.0,
                                   }) {
  if (params.empty())
    params.push_back(25.0);
  double upper =
      params.size() > 1 ? params[1] : std::numeric_limits<double>::infinity();
  int count(0);
  for (const auto &p : obj.prim) {
    if (std::abs(p.pdg) == 11) {
      double ke = (p.genE - (ELECTRON_MASS));
      if (ke >= params[0] && ke <= upper)
        count++;
    }
  }
  return count == 1;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, single_electron,
                   single_electron);

/**
 * @brief Cut for zero true final state muons above threshold.
 * @details Applied at the GENIE generator level using obj.prim. The
 * kinetic energy is computed from the GENIE genE field.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @param params KE threshold in MeV, defaults to 143.425 MeV.
 * @return true if no muons above threshold.
 */
template <typename T>
bool no_muons(const T &obj, std::vector<double> params = {
                                143.425,
                            }) {
  if (params.empty())
    params.push_back(143.425);
  double upper =
      params.size() > 1 ? params[1] : std::numeric_limits<double>::infinity();
  for (const auto &p : obj.prim) {
    if (std::abs(p.pdg) == 13) {
      double ke = (p.genE - (MUON_MASS));
      if (ke >= params[0] && ke <= upper)
        return false;
    }
  }
  return true;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_muons, no_muons);

/**
 * @brief Cut for zero true final state protons within kinetic energy range.
 * @details Applied at the GENIE generator level using obj.prim. The
 * kinetic energy is computed from the GENIE genE field.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @param params KE lower threshold in MeV (defaults to 50 MeV), and optional
 * upper threshold in MeV.
 * @return true if no protons within kinetic energy range.
 */
template <typename T>
bool no_protons(const T &obj, std::vector<double> params = {
                                  50.0,
                              }) {
  if (params.empty())
    params.push_back(50.0);
  double upper =
      params.size() > 1 ? params[1] : std::numeric_limits<double>::infinity();
  for (const auto &p : obj.prim) {
    if (p.pdg == 2212) {
      double ke = (p.genE - (PROTON_MASS));
      if (ke >= params[0] && ke <= upper)
        return false;
    }
  }
  return true;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_protons, no_protons);

/**
 * @brief Cut to select interactions with more than one proton within kinetic
 * energy range.
 * @details Applied at the GENIE generator level using obj.prim. The
 * kinetic energy is computed from the GENIE genE field.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @param params KE lower threshold in MeV (defaults to 50 MeV), and optional
 * upper threshold in MeV.
 * @return true if more than one proton within kinetic energy range.
 */
template <typename T>
bool multiproton(const T &obj, std::vector<double> params = {
                                   50.0,
                               }) {
  if (params.empty())
    params.push_back(50.0);
  double upper =
      params.size() > 1 ? params[1] : std::numeric_limits<double>::infinity();
  int count(0);
  for (const auto &p : obj.prim) {
    if (p.pdg == 2212) {
      double ke = (p.genE - (PROTON_MASS));
      if (ke >= params[0] && ke <= upper)
        count++;
    }
  }
  return count > 1;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, multiproton, multiproton);

/**
 * @brief Cut to select interactions with at least one neutral pion.
 * @details Applied at the GENIE generator level using obj.prim. Only
 * primary particles (start_process == 0) are considered.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @param params energy threshold in MeV, defaults to 0 MeV (no threshold).
 * @return true if at least one pi0 is present.
 */
template <typename T>
bool at_least_one_pi0(const T &obj, std::vector<double> params = {
                                        0.0,
                                    }) {
  for (const auto &p : obj.prim) {
    if (p.start_process != 0)
      continue;
    if (p.pdg == 111) {
      double energy = p.genE;
      if (energy >= params[0])
        return true;
    }
  }
  return false;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, at_least_one_pi0,
                   at_least_one_pi0);

/**
 * @brief Containment cut restricted to true track-like final state
 * primaries, at the GENIE generator level.
 * @details Generator-level analog of cuts::track_containment_cut
 * (cuts.h), which requires every SPINE-reconstructed particle classified
 * as a track (semantic_type == 1) to be contained. There is no semantic
 * classification on a GENIE primary, so track-like species are identified
 * by PDG code instead (muon: 13, charged pion: 211, proton: 2212).
 * Applied to every matching primary in obj.prim regardless of energy or
 * primary/secondary status, mirroring how the reco/true cut applies to
 * every track-classified particle rather than only final-state-signal
 * ones. Containment is read from p.cont_tpc (containment within the TPC
 * that would read out the particle's hits, matching the "track" wording
 * of the cut name); p.contained (full-cryostat containment) is also
 * available on the primary if a looser/stricter definition is desired.
 * @tparam T the type of the object to apply the cut on.
 * @param obj the SRTrueInteraction to apply the cut on.
 * @return true if every true track-like primary is contained.
 */
template <typename T> bool track_containment_cut(const T &obj) {
  for (const auto &p : obj.prim) {
    const int pdg = std::abs(p.pdg);
    if ((pdg == 13 || pdg == 211 || pdg == 2212) && !p.cont_tpc)
      return false;
  }
  return true;
}
REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, track_containment_cut,
                   track_containment_cut);

} // namespace mctruth
#endif
