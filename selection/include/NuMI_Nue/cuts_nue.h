/**
 * @file cuts_eon2024.h
 * @brief Header file for definitions of analysis cuts specific to the eon2024
 * analysis.
 * @details This file contains definitions of analysis cuts which can be used
 * to select interactions specific to the eon2024 analysis. The cuts are
 * intended to be used in conjunction with the generic cuts defined in cuts.h.
 * @author eeller@fnal.gov
*/
#ifndef CUTS_MUON2024_H
#define CUTS_MUON2024_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "include/utilities.h"
#include "include/selectors.h"

/**
 * @namespace cuts::eon2024
 * @brief Namespace for organizing cuts specific to the eon2024 analysis.
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions specific to the eon2024 analysis. Each cut is implemented as
 * a function which takes an interaction object as an argument and returns a
 * boolean. The function should be templated on the type of interaction object if
 * the cut is intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the cuts
 * namespace, which is used for organizing generic cuts which act on interactions.
 */
namespace cuts::nue2024
{
    std::vector<double> flash_times;
    static const bool initialized = [] {
        if constexpr(!BEAM_IS_NUMI)
            flash_times = {0.0,1.6};
        else
            flash_times = {0.0,9.6};
        return true;
    }();

    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), and 1muNp
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), and 1muNp topological cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 1muNp topological cut.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
    bool quality_cuts(const T & obj) { 
            size_t i(selectors::leading_particle_index(obj, 1));
            size_t p(selectors::leading_particle_index(obj, 4));

            return obj.particles[i].axial_spread >0.02 && obj.particles[i].directional_spread < 0.25 && obj.particles[i].vertex_distance <7.5 && obj.particles[p].pid_scores[4] >0.6 && obj.particles[p].pid_scores[2] <0.04;
        }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, quality_cuts, quality_cuts);
    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), and 1muNp
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), and 1muNp topological cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 1muNp topological cut.
     * @note This cut is intended to be used for the muon2024 analysis.
     */
    template<class T>
    bool axial_spread_cuts(const T & obj) { 
            size_t i(selectors::leading_particle_index(obj, 1));

            return obj.particles[i].axial_spread >0.02;
        }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, axial_spread_cuts, axial_spread_cuts);

    template<class T>
    bool directional_spread_cuts(const T & obj) { 
            size_t i(selectors::leading_particle_index(obj, 1));

            return obj.particles[i].directional_spread < 0.25;
        }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, directional_spread_cuts, directional_spread_cuts);

    template<class T>
    bool vertex_distance_cuts(const T & obj) { 
            size_t i(selectors::leading_particle_index(obj, 1));

            return obj.particles[i].vertex_distance <7.5;
        }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, vertex_distance_cuts, vertex_distance_cuts);

    template<class T>
    bool proton_pid_cuts(const T & obj) { 
            size_t p(selectors::leading_particle_index(obj, 4));

            return obj.particles[p].pid_scores[4] >0.6;
        }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, proton_pid_cuts, proton_pid_cuts);
    
    template<class T>
    bool muon_pid_cuts(const T & obj) { 
            size_t p(selectors::leading_particle_index(obj, 4));

            return obj.particles[p].pid_scores[2] <0.04;
        }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, muon_pid_cuts, muon_pid_cuts);
    /**
     * @brief Apply a 1e1p topological (final state) cut.
     * @details The interaction est have a topology matching 1e1p as defined by
     * the conditions in the @ref utilities::count_primaries() function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a 1e1p topology.
     * @note This cut is intended to be used for the eon2024 analysis.
     */
    template<class T>
    bool topological_1e1p_cut(const T & obj)
    {
        std::vector<uint32_t> c(utilities::count_primaries(obj));
        return c[0] == 0 && c[1] == 1 && c[2] == 0 && c[3] == 0 && c[4] == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, topological_1e1p_cut, topological_1e1p_cut);

    /**
     * @brief Apply a 1eNp topological (final state) cut.
     * @details The interaction est have a topology matching 1eNp as defined by
     * the conditions in the @ref utilities::count_primaries() function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a 1eNp topology.
     * @note This cut is intended to be used for the eon2024 analysis.
     */
    template<class T>
    bool topological_1eNp_cut(const T & obj)
    {
        std::vector<uint32_t> c(utilities::count_primaries(obj));
        return c[0] == 0 && c[1] == 1 && c[2] == 0 && c[3] == 0 && c[4] >= 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, topological_1eNp_cut, topological_1eNp_cut);
    
    /**
     * @brief Apply a 1eX topological (final state) cut.
     * @details The interaction est have a topology matching 1eX as defined by
     * the conditions in the @ref utilities::count_primaries() function.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a 1eX topology.
     * @note This cut is intended to be used for the eon2024 analysis.
     */
    template<class T>
    bool topological_1eX_cut(const T & obj)
    {
        std::vector<uint32_t> c(utilities::count_primaries(obj));
        return c[1] == 1 && (c[0] > 0 || c[2] > 0 || c[3] > 0 || c[4] > 0);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, topological_1eX_cut, topological_1eX_cut);

    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), and 1e1p
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), and 1e1p topological cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 1e1p topological cut.
     * @note This cut is intended to be used for the eon2024 analysis.
     */
    template<class T>
    bool all_1e1p_cut(const T & obj) { return fiducial_cut<T>(obj) && track_containment_cut<T>(obj) && flash_cut<T>(obj, flash_times) && topological_1e1p_cut<T>(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, all_1e1p_cut, all_1e1p_cut);

    /**
     * @brief Apply a fiducial volume, flash time (BNB), and 1e1p
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, flash time (BNB), and
     * 1e1p topological cut on the interaction using the logical "and" of each
     * previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, flash time,
     * and 1e1p topological cut.
     * @note This cut is intended to be used for the eon2024 analysis.
     */
    template<class T>
    bool all_1e1p_no_containment_cut(const T & obj) { return fiducial_cut<T>(obj) && flash_cut<T>(obj, flash_times) && topological_1e1p_cut<T>(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, all_1e1p_no_containment_cut, all_1e1p_no_containment_cut);

    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), and 1eNp
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), and 1eNp topological cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 1eNp topological cut.
     * @note This cut is intended to be used for the eon2024 analysis.
     */
    template<class T>
    bool all_1eNp_cut(const T & obj) { return fiducial_cut<T>(obj) && track_containment_cut<T>(obj) && flash_cut<T>(obj, flash_times) && topological_1eNp_cut<T>(obj); }//track_
    REGISTER_CUT_SCOPE(RegistrationScope::Both, all_1eNp_cut, all_1eNp_cut);

    /**
     * @brief Apply a fiducial volume, flash time (BNB), and 1eNp
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, flash time (BNB), and
     * 1eNp topological cut on the interaction using the logical "and" of each
     * previously defined cut.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, flash time,
     * and 1eNp topological cut.
     * @note This cut is intended to be used for the eon2024 analysis.
     */
    template<class T>
    bool all_1eNp_no_containment_cut(const T & obj) { return fiducial_cut<T>(obj) && flash_cut<T>(obj, flash_times) && topological_1eNp_cut<T>(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, all_1eNp_no_containment_cut, all_1eNp_no_containment_cut);

    /**
     * @brief Apply a fiducial volume, containment, flash time (BNB), and 1eX
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, containment, flash time
     * (BNB), and 1eX topological cut on the interaction using the logical "and"
     * of each previously defined cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * flash time, and 1eX topological cut.
     * @note This cut is intended to be used for the eon2024 analysis.
     */
    template<class T>
    bool all_1eX_cut(const T & obj) { return fiducial_cut<T>(obj) && track_containment_cut<T>(obj) && flash_cut<T>(obj, flash_times) && topological_1eX_cut<T>(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, all_1eX_cut, all_1eX_cut);

    /**
     * @brief Apply a fiducial volume, flash time (BNB), and 1eX
     * topological cut (logical "and" of each).
     * @details This function applies a fiducial volume, flash time (BNB), and
     * 1eX topological cut on the interaction using the logical "and" of each
     * previously defined cut.
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, flash time,
     * and 1eX topological cut.
     * @note This cut is intended to be used for the eon2024 analysis.
     */
    template<class T>
    bool all_1eX_no_containment_cut(const T & obj) { return fiducial_cut<T>(obj) && flash_cut<T>(obj, flash_times) && topological_1eX_cut<T>(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, all_1eX_no_containment_cut, all_1eX_no_containment_cut);

    /**
     * @brief Apply a cut to select the 1e1p signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1e1p signal.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1e1p topological cut.
     * @note This cut is intended to be used for the eon2024 analysis for
     * defining the signal.
     */
    template<class T>
    bool signal_1e1p(const T & obj) { return neutrino(obj) && fiducial_cut(obj) && containment_cut(obj) && topological_1e1p_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, signal_1e1p, signal_1e1p);

    /**
     * @brief Apply a cut to select the 1e1p signal.
     * @details This function applies a cut on the final state, and fiducial
     * volume, of the interaction. This is the "true" 1e1p signal.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, and 1e1p
     * topological cut.
     * @note This cut is intended to be used for the eon2024 analysis for
     * defining the signal.
     */
    template<class T>
    bool signal_1e1p_no_containment(const T & obj) { return neutrino(obj) && fiducial_cut(obj) && topological_1e1p_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, signal_1e1p_no_containment, signal_1e1p_no_containment);

    /**
     * @brief Apply a cut to select the 1e1p non-signal.
     * @details This function applies a cut on the final state, and fiducial volume,
     * and containment of the interaction. This is the "true" 1e1p non-signal 
     * (1e1p topology, but not signal).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1e1p topological cut.
     * @note This cut is intended to be used for the eon2024 analysis for
     * defining a complement to the signal.
     */
    template<class T>
    bool nonsignal_1e1p(const T & obj) { return neutrino(obj) && !(fiducial_cut(obj) && containment_cut(obj)) && topological_1e1p_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, nonsignal_1e1p, nonsignal_1e1p);

    /**
     * @brief Apply a cut to select the 1e1p non-signal.
     * @details This function applies a cut on the final state, and fiducial
     * volume, of the interaction. This is the "true" 1e1p non-signal 
     * (1e1p topology, but not signal).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, and 1e1p
     * topological cut.
     * @note This cut is intended to be used for the eon2024 analysis for
     * defining a complement to the signal.
     */
    template<class T>
    bool nonsignal_1e1p_no_containment(const T & obj) { return neutrino(obj) && !fiducial_cut(obj) && topological_1e1p_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, nonsignal_1e1p_no_containment, nonsignal_1e1p_no_containment);

    /**
     * @brief Apply a cut to select the 1eNp signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1eNp signal.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1eNp topological cut.
     * @note This cut is intended to be used for the eon2024 analysis for
     * defining the signal.
     */
    template<class T>
    bool signal_1eNp(const T & obj) { return neutrino(obj) && fiducial_cut(obj) && track_containment_cut(obj) && topological_1eNp_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, signal_1eNp, signal_1eNp);

    /**
     * @brief Apply a cut to select the 1eNp signal.
     * @details This function applies a cut on the final state, and fiducial
     * volume, of the interaction. This is the "true" 1eNp signal.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, and 1eNp
     * topological cut.
     * @note This cut is intended to be used for the eon2024 analysis for
     * defining the signal.
     */
    template<class T>
    bool signal_1eNp_no_containment(const T & obj) { return neutrino(obj) && fiducial_cut(obj) && topological_1eNp_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, signal_1eNp_no_containment, signal_1eNp_no_containment);

    /**
     * @brief Apply a cut to select the 1eNp non-signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1eNp non-signal
     * (1eNp topology, but not signal).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1eNp topological cut.
     * @note This cut is intended to be used for the eon2024 analysis for
     * defining a complement to the signal.
     */
    template<class T>
    bool nonsignal_1eNp(const T & obj) { return neutrino(obj) && !(fiducial_cut(obj) && containment_cut(obj)) && topological_1eNp_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, nonsignal_1eNp, nonsignal_1eNp);

    /**
     * @brief Apply a cut to select the 1eNp non-signal.
     * @details This function applies a cut on the final state, and fiducial
     * volume, of the interaction. This is the "true" 1eNp non-signal 
     * (1eNp topology, but not signal).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, and 1eNp
     * topological cut.
     * @note This cut is intended to be used for the eon2024 analysis for
     * defining a complement to the signal.
     */
    template<class T>
    bool nonsignal_1eNp_no_containment(const T & obj) { return neutrino(obj) && !fiducial_cut(obj) && topological_1eNp_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, nonsignal_1eNp_no_containment, nonsignal_1eNp_no_containment);

    /**
     * @brief Apply a cut to select the 1eX signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1eX signal
     * definition.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1eX topological cut.
     * @note This cut is intended to be used for the eon2024 analysis for
     * defining the signal.
     */
    template<class T>
    bool signal_1eX(const T & obj) { return neutrino(obj) && fiducial_cut(obj) && containment_cut(obj) && topological_1eX_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, signal_1eX, signal_1eX);

    /**
     * @brief Apply a cut to select the 1eX signal.
     * @details This function applies a cut on the final state, and fiducial
     * volume, of the interaction. This is the "true" 1eX signal definition.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, and 1eX
     * topological cut.
     * @note This cut is intended to be used for the eon2024 analysis for
     * defining the signal.
     */
    template<class T>
    bool signal_1eX_no_containment(const T & obj) { return neutrino(obj) && fiducial_cut(obj) && topological_1eX_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, signal_1eX_no_containment, signal_1eX_no_containment);

    /**
     * @brief Apply a cut to select the 1eX non-signal.
     * @details This function applies a cut on the final state, fiducial volume,
     * and containment of the interaction. This is the "true" 1eX non-signal
     * (1eX topology, but not signal).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, containment,
     * and 1eX topological cut.
     * @note This cut is intended to be used for the eon2024 analysis for
     * defining a complement to the signal.
     */
    template<class T>
    bool nonsignal_1eX(const T & obj) { return neutrino(obj) && !(fiducial_cut(obj) && containment_cut(obj)) && topological_1eX_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, nonsignal_1eX, nonsignal_1eX);

    /**
     * @brief Apply a cut to select the 1eX non-signal.
     * @details This function applies a cut on the final state, and fiducial
     * volume, of the interaction. This is the "true" 1eX non-signal 
     * (1eX topology, but not signal).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction passes the fiducial volume, and 1eX
     * topological cut.
     * @note This cut is intended to be used for the eon2024 analysis for
     * defining a complement to the signal.
     */
    template<class T>
    bool nonsignal_1eX_no_containment(const T & obj) { return neutrino(obj) && !fiducial_cut(obj) && topological_1eX_cut(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, nonsignal_1eX_no_containment, nonsignal_1eX_no_containment);
}
#endif // CUTS_MUON2024_H
