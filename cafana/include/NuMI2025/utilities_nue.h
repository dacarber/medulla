
/**
 * @file utilities_ccccpi0ana_phase.h
 * @brief Header file for definitions of utility functions for supporting
 * analysis variables and cuts.
 * @details This file contains definitions of utility functions which are used
 * to support the implementation of analysis variables and cuts. These functions
 * are intended to be used to simplify the implementation of variables and cuts
 * by providing common functionality which can be reused across multiple
 * variables and cuts.
 * @author lkashur@colostate.edu
 */
#ifndef UTILITIES_NUE_H
#define UTILITIES_NUE_H
#include <iostream>
#include <vector>
#include <TVector3.h>
#include "include/cuts.h"
//#include "include/beaminfo.h"

#define MIN_PHOTON_ENERGY 25
#define MIN_ELECTRON_ENERGY 70
#define MIN_MUON_ENERGY 25
#define MIN_PION_ENERGY 25
#define MIN_PROTON_ENERGY 40
#define MIN_KAON_ENERGY 25


struct truth_inter {
  int num_primary_electrons;
  int num_primary_electrons_thresh;
  int num_primary_protons;
  int num_primary_protons_thresh;
  int num_primary_photons;
  int num_primary_photons_thresh;
  int num_primary_muons;
  int num_primary_muons_thresh;
  int num_primary_pions;
  int num_primary_pions_thresh;
  double transverse_momentum_mag;
  bool is_fiducial;
  bool has_contained_tracks;
  bool is_neutrino;
  bool is_cc;
  double leading_electron_energy;
  double leading_electron_conv_dist;
  double leading_electron_axial_spread;
  double leading_electron_dir_spread;
  double leading_electron_softmax;
  double leading_electron_primary;
  double leading_proton_energy;
  double leading_proton_softmax;
  double leading_proton_muon_softmax;
  double leading_proton_pion_softmax;
  double leading_proton_primary;

};

struct reco_inter {
  double transverse_momentum_mag;
  double leading_electron_energy;
  double leading_electron_conv_dist;
  double leading_electron_axial_spread;
  double leading_electron_dir_spread;
  double leading_electron_softmax;
  double leading_electron_primary;
  double leading_proton_energy;
  double leading_proton_softmax;
  double leading_proton_muon_softmax;
  double leading_proton_pion_softmax;
  double leading_proton_primary;
};


/**
 * @namespace utilities_ccpi0ana_phase
 * @brief Namespace for organizing utility functions for supporting analysis
 * variables and cuts.
 * @details This namespace is intended to be used for organizing utility
 * functions which are used to support the implementation of analysis variables
 * and cuts. These functions are intended to be used to simplify the
 * implementation of variables and cuts by providing common functionality which
 * can be reused across multiple variables and cuts.
 * @note The namespace is intended to be used in conjunction with the
 * vars and cuts namespaces, which are used for organizing variables and cuts
 * which act on interactions.
 */
namespace utilities_nue
{
    /**
     * @brief Check if the particle meets final state signal requirements.
     * @details Particle must be primary and have an energy above threshold.
     * Muons must have a length of at least 50 cm (143.425 MeV), pions
     * an energy of at least 25 MeV, and photons an energy of at least 25 MeV.
     * @tparam T the type of particle (true or reco).
     * @param p the particle to check.
     * @return true if the particle is a final state signal particle.
     */
    template<class T>
        bool final_state_signal(const T & p)
        {
	  bool passes(false);
	  if(p.is_primary)
	  {  
            if(PIDFUNC(p) == 0 && p.ke >= MIN_PHOTON_ENERGY) passes = true; // Photons
            if(PIDFUNC(p) == 1 && p.ke >= MIN_ELECTRON_ENERGY) passes = true; // Electrons
            if(PIDFUNC(p) == 2 && p.ke >= MIN_MUON_ENERGY) passes = true; // Muons
	    			if(PIDFUNC(p) == 3 && p.ke >= MIN_PION_ENERGY) passes = true; // Pions
            if(PIDFUNC(p) == 4 && p.ke >= MIN_PROTON_ENERGY) passes = true; // Protons
            if(PIDFUNC(p) == 5 && p.ke >= MIN_KAON_ENERGY) passes = true; // Kaons 
	  }
          return passes;
	}

    /**
     * @brief Count the primaries of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to find the topology of.
     * @return the count of primaries of each particle type within the ineraction.
     */
    template<class T>
        std::vector<uint32_t> count_primaries(const T & obj)
        {
	  std::vector<uint32_t> counts(6, 0);
	  for(auto &p : obj.particles)
          {

            if(final_state_signal(p))
	    {
	        ++counts[PIDFUNC(p)];
	    }
          }
	  return counts;
	}

    /**
     * @brief Structure for true interaction information.
     * @details This structure is needed because true neutral
     * pions are not kept track of in SPINE labeling.  Instead,
     * photons belonging to neutral pions are grouped in pairs.
     * @tparam T the type of interaction (true).
     * @param obj the interaction to select on.
     * @return a truth_inter structure.
     * @note This structure is intented to be used for the ccpi0ana analysis. 
     */
    template<class T> 
      truth_inter truth_interaction_info(const T & obj)
      {
	// Initialize struct
	truth_inter s;
	  
	// Initialize relevant TVector3s
	TVector3 vertex(obj.vertex[0], obj.vertex[1], obj.vertex[2]);
	TVector3 beamdir;
	if constexpr(!BEAM_IS_NUMI){
	    beamdir.SetX(0);
	    beamdir.SetY(0);
	    beamdir.SetZ(1);
	  }
	else{
	  beamdir.SetX(315.120380 + vertex[0]);
	  beamdir.SetY(33.644912 + vertex[1]);
	  beamdir.SetZ(733.632532 + vertex[2]);
	}
	//TVector3 beamdir(BEAMDIR);
	  
	// Initialize output variables
	int primary_electron_count(0);
  int primary_electron_count_thresh(0);
  int primary_proton_count(0);
  int primary_proton_count_thresh(0);
  int primary_photon_count(0);
  int primary_photon_count_thresh(0);
  int primary_muon_count(0);
  int primary_muon_count_thresh(0);
  int primary_pion_count(0);
  int primary_pion_count_thresh(0);
  bool is_fiducial(false);
  bool has_contained_tracks(false);
  bool is_neutrino(false);
  bool is_cc(false);
  double transverse_momentum_mag;
  double leading_electron_energy;
  double leading_electron_conv_dist;
  double leading_electron_axial_spread;
  double leading_electron_dir_spread;
  double leading_electron_softmax;
  double leading_electron_primary;
  double leading_proton_energy;
  double leading_proton_softmax;
  double leading_proton_muon_softmax;
  double leading_proton_pion_softmax;
  double leading_proton_primary;
  double pT0(0), pT1(0), pT2(0);
	// Particle loop
	size_t leading_electron_index(0);
	size_t leading_proton_index(0);
	double max_electron_ke(-99999);
	for(size_t i(0); i < obj.particles.size(); ++i)
	{
	  const auto & p = obj.particles[i];
	  // Primaries
	  if(p.is_primary)
	  {
	    // Transverse momentum calculation
	    TVector3 _p(p.momentum[0], p.momentum[1], p.momentum[2]);
	    TVector3 pL = _p.Dot(beamdir) * beamdir;
	    TVector3 pT = _p - pL;
	    pT0 += pT[0];
	    pT1 += pT[1];
	    pT2 += pT[2];
		      
	    // Electrons
	    if(PIDFUNC(p) == 1)
	    {
	      primary_electron_count++;
	      if(p.ke >= MIN_ELECTRON_ENERGY) primary_electron_count_thresh++;
	      if(p.ke > max_electron_ke)
	      {
					max_electron_ke = p.ke;
					leading_electron_index = i;
	      }
	    }
	    // Protons
	    if(PIDFUNC(p) == 4)
	    {
	      primary_proton_count++;
	      if(p.ke >= MIN_PROTON_ENERGY) primary_proton_count_thresh++;
	    }
	    // Protons
	    if(PIDFUNC(p) == 0)
	    {
	      primary_photon_count++;
	      if(p.ke >= MIN_PHOTON_ENERGY) primary_photon_count_thresh++;
	    }
	    if(PIDFUNC(p) == 2)
	    {
	      primary_muon_count++;
	      if(p.ke >= MIN_MUON_ENERGY) primary_muon_count_thresh++;
	    }
	    if(PIDFUNC(p) == 3)
	    {
	      primary_pion_count++;
	      if(p.ke >= MIN_PION_ENERGY) primary_pion_count_thresh++;
	    }
	  } // end primary loop
	} // end particle loop

	
	
	// Obtain info about signal particles, if they exist
	if(primary_electron_count_thresh == 1 && primary_proton_count_thresh == 1 && obj.current_type == 0 && obj.is_fiducial)
	{      
	  // Get leading muon info
	  const auto & electron = obj.particles[leading_electron_index];
	  const auto & proton = obj.particles[leading_proton_index];

	  TVector3 electron_momentum(electron.momentum[0], electron.momentum[1], electron.momentum[2]);

	  s.leading_electron_energy = leading_electron_energy;
	  s.leading_proton_energy = (proton.ke);
	  s.leading_electron_conv_dist = (electron.vertex_distance);
	  s.leading_electron_axial_spread =  (electron.axial_spread);
	  s.leading_electron_dir_spread =  (electron.directional_spread);
	  s.leading_electron_softmax=  (electron.pid_scores[0]);
	  s.leading_electron_primary =  (electron.primary_scores[1]);
	  s.leading_proton_softmax =  (proton.pid_scores[4]);
	  s.leading_proton_muon_softmax =  (proton.pid_scores[2]);
	  s.leading_proton_pion_softmax =  (proton.pid_scores[3]);
	  s.leading_proton_primary =  (proton.primary_scores[1]);
	        
	} // end signal 
	else
	{
		s.leading_electron_energy = -5;
	  s.leading_proton_energy = -5;
	  s.leading_electron_conv_dist = -5;
	  s.leading_electron_axial_spread =  -5;
	  s.leading_electron_dir_spread =  -5;
	  s.leading_electron_softmax=  -5;
	  s.leading_electron_primary =  -5;
	  s.leading_proton_softmax =  -5;
	  s.leading_proton_muon_softmax =  -5;
	  s.leading_proton_pion_softmax =  -5;
	  s.leading_proton_primary = -5;
	  
	}
	
	s.num_primary_electrons= primary_electron_count;
	s.num_primary_electrons_thresh = primary_electron_count_thresh;
	s.num_primary_protons = primary_proton_count;
	s.num_primary_protons_thresh = primary_proton_count_thresh;
	s.num_primary_photons= primary_photon_count;
	s.num_primary_photons_thresh = primary_photon_count_thresh;
	s.num_primary_muons = primary_muon_count;
	s.num_primary_muons_thresh = primary_muon_count_thresh;
	s.num_primary_pions = primary_pion_count;
	s.num_primary_pions_thresh = primary_pion_count_thresh;
	s.transverse_momentum_mag = sqrt(pow(pT0, 2) + pow(pT1, 2) + pow(pT2, 2));
	s.is_fiducial = obj.is_fiducial;
	s.has_contained_tracks = cuts::track_containment_cut<T>(obj);
	if(obj.nu_id > -1) is_neutrino = true;
	s.is_neutrino = is_neutrino;
	if(obj.current_type == 0) is_cc = true;
	s.is_cc = is_cc;

	return s;
      }

    /**
     * @brief Structure for reco numu cc pi0 interaction information.
     * @details This structure stores information about reconstructed
     * muon and neutral pion.
     * @tparam T the type of interaction (reco).
     * @param obj the interaction to select on.
     * @return a reco_pi0 structure.
     * @note This structure is intented to be used for the ccpi0ana analysis.
     */
    template<class T> 
      reco_inter reco_interaction_info(const T & obj)
      {
	// Initialize structure
      	int primary_electron_count(0);
  int primary_electron_count_thresh(0);
  int primary_proton_count(0);
  int primary_proton_count_thresh(0);
	reco_inter s;
	  
	// Initialize relevant TVector3s
	TVector3 vertex(obj.vertex[0], obj.vertex[1], obj.vertex[2]);
        TVector3 beamdir;
        if constexpr(!BEAM_IS_NUMI){
            beamdir.SetX(0);
            beamdir.SetY(0);
            beamdir.SetZ(1);
          }
        else{
          beamdir.SetX(315.120380 + vertex[0]);
          beamdir.SetY(33.644912 + vertex[1]);
          beamdir.SetZ(733.632532 + vertex[2]);
        }
	//TVector3 beamdir(BEAMDIR);

	// Initialize output variables
	double pT0(0), pT1(0), pT2(0);
	size_t leading_electron_index(0);
	size_t leading_proton_index(0);
	double max_electron_ke(-99999);
	double max_proton_ke(-99999);

	// Loop over particles
	for(size_t i(0); i < obj.particles.size(); ++i)
	{
	    const auto & p = obj.particles[i];
	    if(!p.is_primary) continue; // Primaries

	    // Transverse momentum calculation
	    TVector3 _p(p.momentum[0], p.momentum[1], p.momentum[2]);
	    TVector3 pL = _p.Dot(beamdir) * beamdir;
	    TVector3 pT = _p - pL;
	    pT0 += pT[0];
	    pT1 += pT[1];
	    pT2 += pT[2];
	    // Muons
	    if(PIDFUNC(p) == 1)
	    {
		if(p.ke > max_electron_ke)
		{
		    max_electron_ke = p.ke;
		    leading_electron_index = i;
		}
	    }

	    // Photons
	    if(PIDFUNC(p) == 4)
	    {
				if(p.ke > max_proton_ke)
				{
				    max_proton_ke = p.ke;
				    leading_proton_index = i;
				}
			} // end particle loop
	 const auto & electron = obj.particles[leading_electron_index];
	 const auto & proton = obj.particles[leading_proton_index];

	// Fill struct
	s.transverse_momentum_mag = sqrt(pow(pT0, 2) + pow(pT1, 2) + pow(pT2, 2));
	s.leading_electron_energy = max_electron_ke;
  s.leading_proton_energy = (max_proton_ke);
  s.leading_electron_conv_dist = (electron.vertex_distance);
  s.leading_electron_axial_spread =  (electron.axial_spread);
  s.leading_electron_dir_spread =  (electron.directional_spread);
  s.leading_electron_softmax=  (electron.pid_scores[0]);
  s.leading_electron_primary =  (electron.primary_scores[1]);
  s.leading_proton_softmax =  (proton.pid_scores[4]);
  s.leading_proton_muon_softmax =  (proton.pid_scores[2]);
  s.leading_proton_pion_softmax =  (proton.pid_scores[3]);
  s.leading_proton_primary =  (proton.primary_scores[1]);



	return s;
      }

	}
}
#endif // UTILITIES_NUE_H
