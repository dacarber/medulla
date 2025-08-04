/**
 * @file vars_muon2024.h
 * @brief Header file for definitions of analysis variables specific to the
 * muon2024 analysis.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from objs specific to the muon2024
 * analysis. Each variable is implemented as a function which takes an
 * obj object as an argument and returns a double. These are the
 * building blocks for producing high-level plots of the selected objs.
 * @author mueller@fnal.gov
 */
#ifndef VARS_MUON2024_H
#define VARS_MUON2024_H

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"

#include "include/selectors.h"
#include "include/framework.h"
#include "include/cuts.h"
#include "include/muon2024/cuts_muon2024.h"

/**
 * @namespace vars::muon2024
 * @brief Namespace for organizing variables specific to the muon2024 analysis.
 * @details This namespace is intended to be used for organizing variables which
 * act on objs specific to the muon2024 analysis. Each variable is
 * implemented as a function which takes an obj object as an argument
 * and returns a double. The function should be templated on the type of
 * obj object if the variable is intended to be used on both true and
 * reconstructed objs.
 * @note The namespace is intended to be used in conjunction with the vars
 * namespace, which is used for organizing generic variables which act on
 * objs.
 */
namespace vars::nue
{
    /**
     * @brief Variable for enumerating obj categories.
     * @details This variable provides a basic categorization of objs
     * using only signal, neutrino background, and cosmic background as the
     * three categories.
     * 0: 1mu1p (contained and fiducial)
     * 1: 1mu1p (not contained or not fiducial)
     * 2: 1muNp (N > 1, contained and fiducial)
     * 3: 1muNp (N > 1, not contained or fiducial)
     * 4: 1muX (not 1muNp, contained and fiducial)
     * 5: 1muX (not 1muNp, not contained or fiducial)
     * 6: Other CC nu
     * 7: Other NC nu
     * 8: Cosmic
     * @tparam T the type of obj (true or reco).
     * @param obj The obj to apply the variable on.
     * @return the enumerated category of the obj.
    */
    template<class T>
    double category(const T & obj)
    {   
        double cat(11);
        if(obj.nu_id >= 0){
            std::vector<uint32_t> counts(utilities::count_primaries(obj));
            if(counts[1] == 1 && counts[2] == 0 )
                {
                    if(counts[0] == 0 && counts[3] == 0 && counts[4] == 1 && cuts::track_containment_cut(obj) && obj.is_fiducial) cat = 0; 
                    else if(counts[0] == 0 && counts[3] == 0 && counts[4] == 1 && cuts::track_containment_cut(obj)) cat = 8;
                    else if(counts[0] == 0 && counts[3] == 0 && counts[4] == 1) cat = 9;
                    else if(counts[0] == 0 && counts[3] == 0 && counts[4] > 1 && cuts::track_containment_cut(obj) && obj.is_fiducial) cat = 1;
                    else if(counts[0] == 0 && counts[3] == 0 && counts[4] > 1 && cuts::track_containment_cut(obj)) cat = 8;
                    else if(counts[0] == 0 && counts[3] == 0 && counts[4] > 1) cat = 9;
                    else if(counts[0] == 0 && counts[3] == 0 && counts[4] == 0) cat = 2;
                    else if(counts[0] == 0 && counts[3] == 0 && counts[4] == 0) cat = 2;
                    else if(counts[0] == 0 && counts[3] == 1 && counts[4] > 0) cat = 3;
                    else if(obj.current_type == 0) cat = 4;
                    else if(obj.current_type == 1) cat = 5;
                    else cat =10;
                }
                else if(obj.current_type == 0 && counts[2] == 1) cat = 7;
                else if(obj.current_type == 1) cat = 5;
                else if(obj.current_type == 0 && obj.pdg_code == 12) cat = 6;
                else if(obj.current_type == 0 && obj.pdg_code == 14) cat = 7;
                else if(obj.current_type == 1) cat = 5;
                else cat = 10;

        }
        return cat;
    
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, category, category);

    /**
     * @brief 
     * @details 
     * @tparam T the type of obj (true or reco).
     * @param obj The obj to apply the variable on.
     * @return 
    */
    template<class T>
    double topology(const T & obj)
    {   
        double topo(0);
            std::vector<uint32_t> counts(utilities::count_primaries(obj));
            topo = counts[0]*1e8 + counts[1]*1e6 + counts[2]*1e4 + counts[3]*1e2 + counts[4];

        return topo;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, topology, topology);

    template<class T>
        double neutrino_energy(const T & obj)
        {
            double nu_energy(0);
            for(auto &p : obj.particles)
            {
                if (p.is_primary == false) continue;
                if (p.pid == 0){
                    nu_energy+=pvars::ke(p);
                } 
                else if (p.pid == 1){
                    nu_energy+=pvars::ke(p) + ELECTRON_MASS;
                }
                else if (p.pid == 2){
                    nu_energy+=pvars::ke(p) + MUON_MASS;
                }
                else if (p.pid == 3){
                    nu_energy+=pvars::ke(p) + PION_MASS;
                }
                else if (p.pid == 4){
                    nu_energy+=pvars::ke(p) + 40;
                }
            }
            return nu_energy;
        }
    '''
    REGISTER_VAR_SCOPE(RegistrationScope::Both, neutrino_energy, neutrino_energy);
    template<class T>
        double Qsquared(const T & obj)
        {
            double nu_energy = neutrino_energy(obj);
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>){
                std::vector<uint32_t> counts(utilities::count_primaries(obj));
                double nu_energy = obj.energy_init + 40*counts[4];
            }
            
            size_t i(selectors::leading_particle_index(obj, 1));
            double electron_energy = pvars::ke(obj.particles[i]) + ELECTRON_MASS;
            TVector3 p;
            p.SetX(obj.particles[i].momentum[0]);
            p.SetY(obj.particles[i].momentum[1]);
            p.SetZ(obj.particles[i].momentum[2]);
            TVector3 beamdir(0.39431672, 0.04210058, 0.91800973);
            double cos_theta = p.Dot(beamdir)/(p.Mag()*beamdir.Mag());
            double Q = 2*nu_energy*(electron_energy - p.Mag()*cos_theta) - std::pow(ELECTRON_MASS,2);
            return Q;

        }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, Qsquared, Qsquared);
    template<class T>
        double W(const T & obj)
        {
            //if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>){
            //    proton_num = utilities::count_primaries()
            //    nu_energy = obj.energy_init;
            //}
            double MN =37147.393;
            double nu_energy = neutrino_energy(obj);
            size_t i(utilities::leading_particle_index(obj, 1));
            double electron_energy = pvars::ke(obj.particles[i]) + ELECTRON_MASS;
            double Q = Qsquared(obj);
            double W = std::sqrt(std::pow(MN,2)+2*MN*(nu_energy-electron_energy) - Q);
            
            return W;

        }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, W, W);
    template<class T>
        double cos_theta(const T & particle)
        {


            TVector3 p;
            p.SetX(particle.momentum[0]);
            p.SetY(particle.momentum[1]);
            p.SetZ(particle.momentum[2]);
            TVector3 beamdir(0.39431672, 0.04210058, 0.91800973);
            double cos_theta = p.Dot(beamdir)/(p.Mag()*beamdir.Mag());
            
            return cos_theta;

        }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, cos_theta, cos_theta);

    template<class T>
        double leading_proton_softmax(const T & obj)
        {
            size_t i(selectors::leading_proton(obj));
            return obj.particles[i].pid_scores[4];
        }
    template<class T>
        double leading_proton_ke(const T & obj)
        {
            size_t i(selectors::leading_proton(obj));
            return pvars::ke(obj.particles[i]);
        }
    template<class T>
        double leading_proton_muon_softmax(const T & obj)
        {
            size_t i(selectors::leading_proton(obj));
            return obj.particles[i].pid_scores[2];
        }
    template<class T>
        double leading_proton_pion_softmax(const T & obj)
        {
            size_t i(selectors::leading_proton(obj));
            return obj.particles[i].pid_scores[3];
        }
    template<class T>
        double leading_proton_primary_score(const T & obj)
        {
            size_t i(selectors::leading_proton(obj));
            return obj.particles[i].primary_scores[1];
        }

    template<class T>
        double leading_electron_ke(const T & obj)
        {
            size_t i(selectors::leading_electron(obj));
            return pvars::ke(obj.particles[i]);
        }
    template<class T>
        double leading_electron_vertex_distance(const T & obj)
        {   
            
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
            {
               return -5; 
            }
            size_t i(selectors::leading_electron(obj));
            return obj.particles[i].vertex_distance;
            
        }
    template<class T>
        double leading_electron_axial_spread(const T & obj)
        {
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
            {
               return -5;
            }
            size_t i(selectors::leading_electron(obj));
            return obj.particles[i].axial_spread;
        }
    template<class T>
        double leading_electron_directional_spread(const T & obj)
        {
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
            {
               return -5;
            }
            size_t i(selectors::leading_electron(obj));
            return obj.particles[i].directional_spread;
        }
    template<class T>
        double leading_electron_start_straightness(const T & obj)
        {
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
            {
               return -5;
            }
            size_t i(selectors::leading_electron(obj));
            return obj.particles[i].start_straightness;
        }
    template<class T>
        double leading_electron_dedx(const T & obj)
        {  
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
            {
               return -5;
            }
            size_t i(selectors::leading_electron(obj));
            return obj.particles[i].start_dedx;
        }
    template<class T>
        double leading_electron_primary_score(const T & obj)
        {
            size_t i(selectors::leading_electron(obj));
            return obj.particles[i].primary_scores[1];
        }
    template<class T>
        double leading_electron_softmax(const T & obj)
        {
            size_t i(selectors::leading_electron(obj));
            return obj.particles[i].pid_scores[0];
        }
    template<class T>
        double leading_electron_NuMI_polar_angle(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 1));
            double angle(cos_theta(obj.particles[i]));
            return angle;
        }
    template<class T>
        double leading_proton_NuMI_polar_angle(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 4));
            double angle(cos_theta(obj.particles[i]));
            return angle;
        }
    template<class T>
        double leading_electron_NuMI_azimuthal_angle(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 1));
            double angle(NuMI_azimuthal_angle(obj.particles[i]));
            return angle;
        }
    template<class T>
        double leading_proton_NuMI_azimuthal_angle(const T & obj)
        {
            size_t i(utilities::leading_particle_index(obj, 4));
            double angle(NuMI_azimuthal_angle(obj.particles[i]));
            return angle;
        }
    '''
    /**
     * @brief Variable for the opening angle between leading muon and proton.
     * @details The leading muon and proton are defined as the particles with the
     * highest kinetic energy. The opening angle is defined as the arccosine of
     * the dot product of the momentum vectors of the leading muon and proton.
     * @tparam T the type of obj (true or reco).
     * @param obj the obj to apply the variable on.
     * @return the opening angle between the leading muon and
     * proton.
     */
    template<class T>
    double opening_angle(const T & obj)
    {
        size_t mi = selectors::leading_electron(obj);
        size_t pi = selectors::leading_proton(obj);
        if(mi == kNoMatch || pi == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
        else
        {
            auto & m(obj.particles[mi]);
            auto & p(obj.particles[pi]);
            return std::acos(m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, opening_angle, opening_angle);
}
#endif // VARS_MUON2024_H