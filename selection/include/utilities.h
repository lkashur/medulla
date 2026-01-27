/**
 * @file utilities.h
 * @brief Header file for definitions of utility functions acting on
 * interactions.
 * @details This file contains definitions of utility functions which are used
 * to support the implementation of analysis variables and cuts. These functions
 * are intended to be used to simplify the implementation of variables and cuts
 * by providing common functionality which can be reused across multiple
 * variables and cuts.
 * @author mueller@fnal.gov
 */
#ifndef UTILITIES_H
#define UTILITIES_H
#include <vector>
#include <array>

#define PI0_MASS 134.9768

#include "framework.h"
#include "include/particle_variables.h"
#include "include/particle_cuts.h"
#include "include/particle_utilities.h"

struct pi0 {
  double_t vertex_x;
  double_t vertex_y;
  double_t vertex_z;
  size_t leading_shower_index;
  double leading_shower_ke;
  double leading_shower_conv_dist;
  size_t subleading_shower_index;
  double subleading_shower_ke;
  double subleading_shower_conv_dist;
  double showers_costheta;
  double mass;
  double energy;
  double momentum;
  double beam_costheta;
};


/**
 * @namespace utilities
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
namespace utilities
{

    /**
     * @brief Count the primaries of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to find the topology of.
     * @return the count of primaries of each particle type within the
     * interaction.
     */
    template<class T>
        std::vector<uint32_t> count_primaries(const T & obj)
        {
            std::vector<uint32_t> counts(5, 0);
            for(auto &p : obj.particles)
            {
                if(pcuts::final_state_signal(p))
                    ++counts[pvars::pid(p)];
            }
            return counts;
        }

    /**
     * @brief List of good runs for ICARUS Run 2.
     * @details This list contains the run numbers that are considered good
     * for analysis in ICARUS Run 2. The list is used to filter out runs that
     * are known to have issues or are not suitable for analysis.
     * @note This is intended to be removed once we have a proper version-
     * controlled run list committed to sbnana.
     */
    constexpr std::array<unsigned int, 229>
    icarus_good_runs_run2{  9301,  9302,  9303,  9307,  9308,  9309,  9310,  9311,  9312,  9313,
                            9314,  9316,  9317,  9318,  9327,  9328,  9329,  9330,  9332,  9333,
                            9335,  9337,  9338,  9339,  9340,  9341,  9342,  9343,  9344,  9346,
                            9347,  9353,  9354,  9356,  9357,  9358,  9359,  9360,  9361,  9362,
                            9363,  9364,  9365,  9366,  9367,  9380,  9383,  9384,  9385,  9386,
                            9387,  9388,  9389,  9390,  9391,  9392,  9393,  9394,  9409,  9412,
                            9415,/*9435,*/9436,  9437,  9438,  9439,  9441,  9445,  9448,  9450,
                            9458,  9460,  9472,  9473,  9474,  9477,  9478,  9481,  9482,  9499,
                            9504,  9518,  9563,  9565,  9569,  9570,  9583,  9587,  9588,  9589,
                            9590,  9593,  9594,  9595,  9597,  9598,  9599,  9602,  9610,  9626,
                            9627,  9631,  9647,  9648,  9649,  9658,  9672,  9675,  9688,  9689,
                            9690,  9691,  9692,  9693,  9694,  9695,  9696,  9697,  9698,  9699,
                            9700,  9703,  9704,  9705,  9714,  9715,  9716,  9717,  9721,  9723,
                            9724,  9725,  9726,  9728,  9729,  9730,  9731,  9732,  9733,  9734,
                            9735,  9743,  9744,  9745,  9746,  9747,  9750,  9752,  9753,  9755,
                            9758,  9762,  9763,  9764,  9765,  9783,  9788,  9791,  9792,  9793,
                            9794,  9795,  9796,  9807,  9834,  9835,  9837,  9838,  9840,  9841,
                            9844,  9847,  9849,  9851,  9854,  9855,  9860,  9862,  9868,  9870,
                            9892,  9894,  9896,  9897,  9914,  9917,  9919,  9921,  9922,  9924,
                            9925,  9926,  9929,  9932,  9935,  9940,  9941,  9942,  9944,  9945,
                            9946,  9949,  9950,  9951,  9953,  9954,  9956,  9959,  9960,  9961,
                            9970,  9971,  9974,  9977,  9979,  9981,  9982,  9986, 10054, 10059,
                           10061, 10062, 10064, 10065, 10066, 10067, 10084, 10085, 10096, 10097};
    
    /**
     * @brief Boolean function to check if a run is in the list of good runs
     * for ICARUS Run 2.
     * @param run the run number to check.
     * @return true if the run is in the list of good runs, false otherwise.
     */
    bool is_icarus_good_run(unsigned int run)
    {
        return std::find(icarus_good_runs_run2.begin(), icarus_good_runs_run2.end(), run) != icarus_good_runs_run2.end();
    }

    /**
     * @brief Utility function to find the index of the first optical flash
     * @details This function returns the index of the first optical flash in
     * the event, shifted by a given amount. It is intended to be used to find
     * the flash nearest to the trigger of the event. The optional shift is
     * used to allow for tuning to correct for the natural offset of the
     * reconstructed flash time from the trigger time, which is not zero
     * despite all systems being referenced to the trigger. This function uses
     * the `firsttime` field of the optical flash.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @param shift the amount to shift the time of the flash by, default is 0.
     * @return size_t the index of the first optical flash in the event, or
     * kNoMatch if no flash is found.
     */
    template<typename T>
    size_t first_opflash_firsttime(const T & sr, double shift=0.0)
    {
        // This variable returns the index of the first optical flash in the
        // event, shifted by a given amount.
        if(sr.opflashes.empty()) return kNoMatch;
        size_t first_flash_index = 0;
        double min_time_diff = std::abs(sr.opflashes[0].firsttime - shift);
        for(size_t i = 1; i < sr.opflashes.size(); ++i)
        {
            double time_diff = std::abs(sr.opflashes[i].firsttime - shift);
            if(time_diff < min_time_diff)
            {
                min_time_diff = time_diff;
                first_flash_index = i;
            }
        }
        return first_flash_index;
    }

    /**
     * @brief Utility function to find the index of the first optical flash
     * @details This function returns the index of the first optical flash in
     * the event, shifted by a given amount. It is intended to be used to find
     * the flash nearest to the trigger of the event. The optional shift is
     * used to allow for tuning to correct for the natural offset of the
     * reconstructed flash time from the trigger time, which is not zero
     * despite all systems being referenced to the trigger. This function uses
     * the `time` field of the optical flash.
     * @tparam T the top-level record.
     * @param sr the StandardRecord to apply the variable on.
     * @param shift the amount to shift the time of the flash by, default is 0.
     * @return size_t the index of the first optical flash in the event, or
     * kNoMatch if no flash is found.
     */
    template<typename T>
    size_t first_opflash_rawtime(const T & sr, double shift=0.0)
    {
        // This variable returns the index of the first optical flash in the
        // event, shifted by a given amount.
        if(sr.opflashes.empty()) return kNoMatch;
        size_t first_flash_index = 0;
        double min_time_diff = std::abs(sr.opflashes[0].time - shift);
        for(size_t i = 1; i < sr.opflashes.size(); ++i)
        {
            double time_diff = std::abs(sr.opflashes[i].time - shift);
            if(time_diff < min_time_diff)
            {
                min_time_diff = time_diff;
                first_flash_index = i;
            }
        }
        return first_flash_index;
    }

    static std::vector<double> reco_pi0_shower_thresholds = {25.0, 25.0};

    /**
     * @brief Group true pi0 daughters by parent ID
     * @param obj the interaction from which to count
     * @return an unordered map with parent ID as key and daughter particle indices as value
     */
    std::unordered_map< int, std::vector<size_t> > get_true_pi0s(const caf::SRInteractionTruthDLPProxy & obj, bool primaries = true, std::vector<double> params = {0.0,})
    {
        std::unordered_map< int, std::vector<size_t> > true_pi0s;
	
	/**'
	 * Loop over particles, grouping photons
	 * and electrons that share a parent pi0.
	 */
	for(size_t i(0); i < obj.particles.size(); ++i)
	{
	    const auto & p = obj.particles[i];
	    

	    if(primaries == true and !p.is_primary) continue;
	    if(primaries == false and p.is_primary) continue;

	    if(p.parent_pdg_code == 111 && (p.pdg_code == 22 || p.pdg_code == 11 || p.pdg_code == -11))
	    {
	        true_pi0s[p.parent_track_id].push_back(i);
	    }
	}

	/**
	 * Remove pi0 candidates that are either
	 * subthreshold or have fewer than two daughters.
	 */
	std::vector<size_t> bad_pi0_ids;
        for(auto const & _pi0 : true_pi0s)
        {

	  int num_pi0_daughters(0);
	  three_vector pi0_momentum(0,0,0);
	  double pi0_ke(0);
	  for(auto & _pidx : _pi0.second)
	    {
	      const auto & _p = obj.particles[_pidx];
	      three_vector _p_momentum(_p.momentum[0], _p.momentum[1], _p.momentum[2]);
	      std::get<0>(pi0_momentum) += std::get<0>(_p_momentum);
          std::get<1>(pi0_momentum) += std::get<1>(_p_momentum);
          std::get<2>(pi0_momentum) += std::get<2>(_p_momentum);
	      num_pi0_daughters++;
	    }
        double x = std::get<0>(pi0_momentum);
        double y = std::get<1>(pi0_momentum);
        double z = std::get<2>(pi0_momentum);
        double pi0_p = std::sqrt(x*x + y*y + z*z);
	    pi0_ke = std::sqrt(std::pow(PI0_MASS, 2) + std::pow(pi0_p, 2)) - PI0_MASS;
	    
	    // TEST: RELAX SIGNAL DEFINITION
	    //if (pi0_ke < params[0]) bad_pi0_ids.push_back(_pi0.first);
            if(num_pi0_daughters < 2 || pi0_ke < params[0]) bad_pi0_ids.push_back(_pi0.first);
	}
	for(size_t i=0; i<bad_pi0_ids.size(); i++)
	{
            true_pi0s.erase(bad_pi0_ids[i]);
	}
	return true_pi0s;
    }

    /**
     * @brief Count the number of true primary neutral pions in interaction.
     * @param obj the interaction from which to count.
     * @param params the neutral pion kinetic energy threshold to consider.
     * @return the number of true primary neutral pions in interaction.
     */
    double true_primary_pi0_multiplicity(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={0.0,})
    {
        std::unordered_map< int, std::vector<size_t> > true_primary_pi0s = get_true_pi0s(obj, true, params);
	return true_primary_pi0s.size();
    }

    /**
     * @brief Count the number of true nonprimary neutral pions in interaction.
     * @param obj the interaction from which to count.
     * @param params the neutral pion kinetic energy threshold to consider.
     * @return the number of true nonprimary neutral pions in interaction.
     */
    double true_nonprimary_pi0_multiplicity(const caf::SRInteractionTruthDLPProxy & obj, std::vector<double> params={0.0,})
    {
        std::unordered_map< int, std::vector<size_t> > true_nonprimary_pi0s = get_true_pi0s(obj, false, params);
	return true_nonprimary_pi0s.size();
    }

    /**
     * @brief Extract information about neutral pion.
     * @details Extract information about neutral pion from its daughter particles, 
     * assuming only one neutral pion is present in interaction.
     * @param obj the interaction from which to count.
     * @param params the neutral pion kinetic energy threshold to consider.
     * @return the number of true nonprimary neutral pions in interaction.
     */
    template<class T>
    pi0 pi0_info(const T & obj)
    {
        /**
	 * Initialize output structure.
	 */
        pi0 s;
	three_vector momentum(0,0,0);
	double beam_costheta(-5);
	double vertex_x(-5);
	double vertex_y(-5);
	double vertex_z(-5);
	size_t leading_shower_index(kNoMatch);
	double leading_shower_ke(-5);
	double leading_shower_conv_dist(-5);
        size_t subleading_shower_index(kNoMatch);
        double subleading_shower_ke(-5);
        double subleading_shower_conv_dist(-5);
	double showers_costheta(-5);
	double mass(-5);
	double energy(-5);

	/**
	 * Neutrino beam direction and neutrino interaction vertex,
	 * to be used in later calculations.
	 */
	three_vector beamdir(0,0,1);
	three_vector vertex(obj.vertex[0], obj.vertex[1], obj.vertex[2]);

	/**
	 * True neutral pions
	 */
	if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
		       {
			   /**
			    * Retrieve map with true neutral pion daughter indices.
			    */
			   std::unordered_map< int, std::vector<size_t> > true_primary_pi0s = get_true_pi0s(obj, true);
			   
			   /**
			    * Calculate neutral pion momentum from daugher particles.
			    * Ideally the true neutral pion momentum would be stored
			    * in SPINE, but this is not currently the case.
			    */
			   int num_pi0_daughters(0);
			   int num_pi0_photon_daughters(0);
			   std::vector<size_t> pi0_daughter_indices;
			   for(auto const & _pi0 : true_primary_pi0s)
			   {
			       for(auto & _pidx : _pi0.second)
			       {
				   const auto & _p = obj.particles[_pidx];
				   three_vector _p_momentum(_p.momentum[0], _p.momentum[1], _p.momentum[2]);
                   std::get<0>(momentum) += std::get<0>(_p_momentum);
                   std::get<1>(momentum) += std::get<1>(_p_momentum);
                   std::get<2>(momentum) += std::get<2>(_p_momentum);
				   pi0_daughter_indices.push_back(_pidx);
				   num_pi0_daughters++;
				   if(_p.pid == 0) num_pi0_photon_daughters++;
			       }
			   }
			   
			   /**
			    * Calculate neutral pion angle with respect to neutrino beam.
			    */
               double raw_dot = utilities::dot_product(momentum, beamdir);
               double p_mag = utilities::magnitude(momentum);
               double result = (p_mag > 0) ? (raw_dot / p_mag) : 0.0;
			   beam_costheta = result;

			   /**
			    * If neutral pion decays to two photons (the usual case),
			    * retrieve information about leading and subleading photons.
			    */
			   if(num_pi0_photon_daughters == 2)
			   {
			       const auto & pi0_daughter0 = obj.particles[pi0_daughter_indices[0]];
			       const auto & pi0_daughter1 = obj.particles[pi0_daughter_indices[1]];

			       size_t _leading_shower_index;
			       size_t _subleading_shower_index;
			       if(pi0_daughter0.ke > pi0_daughter1.ke)
			       {
				   _leading_shower_index = pi0_daughter_indices[0];
				   _subleading_shower_index = pi0_daughter_indices[1];
			       }
			       else
			       {
				   _leading_shower_index = pi0_daughter_indices[1];
				   _subleading_shower_index = pi0_daughter_indices[0];
			       }
			       
			       /**
				* Leading photon information
				*/
			       const auto & leading_shower = obj.particles[_leading_shower_index];
			       leading_shower_index = _leading_shower_index;
			       leading_shower_ke = leading_shower.ke;
			       three_vector leading_shower_start_point(leading_shower.start_point[0], leading_shower.start_point[1], leading_shower.start_point[2]);
			       three_vector leading_shower_dir(leading_shower.momentum[0], leading_shower.momentum[1], leading_shower.momentum[2]);
			       leading_shower_dir = utilities::normalize(leading_shower_dir);
			       leading_shower_conv_dist = utilities::magnitude(utilities::subtract(vertex, leading_shower_start_point));
			       three_vector leading_shower_momentum(leading_shower.momentum[0], leading_shower.momentum[1], leading_shower.momentum[2]);
			       
			       /**
				* Subleading photon information
				*/
			       const auto & subleading_shower = obj.particles[_subleading_shower_index];
                               subleading_shower_index = _subleading_shower_index;
                               subleading_shower_ke = subleading_shower.ke;
                               three_vector subleading_shower_start_point(subleading_shower.start_point[0], subleading_shower.start_point[1], subleading_shower.start_point[2]);
                               three_vector subleading_shower_dir(subleading_shower.momentum[0], subleading_shower.momentum[1], subleading_shower.momentum[2]);
                               subleading_shower_dir = utilities::normalize(subleading_shower_dir);
                               subleading_shower_conv_dist = utilities::magnitude(utilities::subtract(vertex, subleading_shower_start_point));
                               three_vector subleading_shower_momentum(subleading_shower.momentum[0], subleading_shower.momentum[1], subleading_shower.momentum[2]);

			       /**
				* Information derived from leading and subleading photons.
				*/
			       showers_costheta = utilities::dot_product(leading_shower_dir,subleading_shower_dir);
			       mass = sqrt(2*leading_shower_ke*subleading_shower_ke*(1-showers_costheta));
			       energy = leading_shower_ke + subleading_shower_ke;
			       vertex_x = leading_shower.parent_position[0];
			       vertex_y = leading_shower.parent_position[1];
			       vertex_z = leading_shower.parent_position[2];
			   }
			   
			   /**
			    * Store output.
			    */
			   s.vertex_x = vertex_x;
			   s.vertex_y = vertex_y;
			   s.vertex_z = vertex_z;
			   s.momentum = utilities::magnitude(momentum)/1000.0;
			   // Override if only one pi0 daughter exists
			   if (num_pi0_daughters < 2)
			     s.momentum = -5;
			   s.beam_costheta = beam_costheta;
			   s.leading_shower_index = leading_shower_index;
			   s.leading_shower_ke = leading_shower_ke/1000.0;
			   s.leading_shower_conv_dist = leading_shower_conv_dist;
			   s.subleading_shower_index = subleading_shower_index;
                           s.subleading_shower_ke = subleading_shower_ke/1000.0;
                           s.subleading_shower_conv_dist = subleading_shower_conv_dist;
			   s.showers_costheta = showers_costheta;
			   s.mass = mass;
			   s.energy = energy/1000.0;
		       }
	
	/**
	 * Reco neutral pions
	 */
	else
	{
	    /**
	     * Loop over every possible pair of photons.
	     * Pair that has best agreement with true neutral pion mass
	     * is chosen as belonging to candidate neutral pion.
	     */
	    std::vector<std::pair< std::pair<size_t, size_t>, double> > showers_metric;
	    for(size_t i(0); i < obj.particles.size(); ++i)
	    {
	        /**
		 * Photon 0
		 */
	        const auto & p = obj.particles[i];
		if(!(p.is_primary && p.pid == 0)) continue;
		
		three_vector sh0_start(p.start_point[0], p.start_point[1], p.start_point[2]);
		three_vector sh0_start_dir = utilities::normalize(utilities::subtract(sh0_start, vertex));
		three_vector sh0_dir(p.start_dir[0], p.start_dir[1], p.start_dir[2]);

		for(size_t j(0); j < obj.particles.size(); ++j)
		{
		  if(j == i) continue;
		  
		  /**
		   * Photon 1
		   */
		  const auto & q = obj.particles[j];
		  if(!(q.is_primary && q.pid == 0)) continue;

		  three_vector sh1_start(q.start_point[0], q.start_point[1], q.start_point[2]);
		  three_vector sh1_start_dir = utilities::normalize(utilities::subtract(sh1_start, vertex));
		  three_vector sh1_dir(q.start_dir[0], q.start_dir[1], q.start_dir[2]);

		  double _leading_shower_ke(0);
		  double _subleading_shower_ke(0);

		  if(pvars::calo_ke(p) > pvars::calo_ke(q))
		    {
		      _leading_shower_ke = pvars::calo_ke(p);
		      _subleading_shower_ke = pvars::calo_ke(q);
		    }
		  else
		    {
		      _leading_shower_ke = pvars::calo_ke(q);
		      _subleading_shower_ke = pvars::calo_ke(p);
		    }
		  if(_leading_shower_ke < reco_pi0_shower_thresholds[1] || _subleading_shower_ke < reco_pi0_shower_thresholds[0]) continue;

		  /**
		   * Calculate invariant mass of photon pair.
		   */
		  double _showers_costheta = utilities::dot_product(sh0_start_dir,sh1_start_dir);
		  double _mass = sqrt(2*_leading_shower_ke*_subleading_shower_ke*(1-_showers_costheta));
		  showers_metric.push_back(std::make_pair(std::make_pair(i, j), mass));
		  //std::sort(showers_metric.begin(), showers_metric.end(), [](const std::pair<std::pair<size_t, size_t>, double> &a, const std::pair<std::pair<size_t, size_t>, double> &b)
		  //{ return abs(a.second - 134.9768) < abs(b.second - 134.9768);});
		}
	    }
	    
	    /**
	     * Sort combinations of photon pairs by agreement with true neutral pion mass
	     * and identify leading and subleading photons.
	     */
	    if(!showers_metric.empty())
	    {
	        std::sort(showers_metric.begin(), showers_metric.end(), [](const std::pair<std::pair<size_t, size_t>, double> &a, const std::pair<std::pair<size_t, size_t>, double> &b)
			  { return abs(a.second - PI0_MASS) < abs(b.second - PI0_MASS);});
	      
	        std::pair<size_t, size_t> sh_pair_ids = showers_metric[0].first; // best metric agreement
		if(pvars::calo_ke(obj.particles[sh_pair_ids.first]) > pvars::calo_ke(obj.particles[sh_pair_ids.second]))
		{
		    leading_shower_index = sh_pair_ids.first;
		    subleading_shower_index = sh_pair_ids.second;
		}
		else
		{
		    leading_shower_index = sh_pair_ids.second;
		    subleading_shower_index = sh_pair_ids.first;
		}

		/**
		 * Leading photon information
		 */
		const auto & leading_shower = obj.particles[leading_shower_index];
		leading_shower_ke = pvars::calo_ke(leading_shower);
		three_vector leading_shower_start_point(leading_shower.start_point[0], leading_shower.start_point[1], leading_shower.start_point[2]);
		three_vector leading_shower_dir = utilities::subtract(leading_shower_start_point, vertex);
		leading_shower_dir = utilities::normalize(leading_shower_dir);
		leading_shower_conv_dist = utilities::magnitude(utilities::subtract(vertex, leading_shower_start_point));
		three_vector leading_shower_momentum = utilities::scale(leading_shower_dir, leading_shower_ke);

		/**
		 * Subleading photon information
		 */
		const auto & subleading_shower = obj.particles[subleading_shower_index];
                subleading_shower_ke = pvars::calo_ke(subleading_shower);
                three_vector subleading_shower_start_point(subleading_shower.start_point[0], subleading_shower.start_point[1], subleading_shower.start_point[2]);
                three_vector subleading_shower_dir = utilities::subtract(subleading_shower_start_point, vertex);
                subleading_shower_dir = utilities::normalize(subleading_shower_dir);
                subleading_shower_conv_dist = utilities::magnitude(utilities::subtract(vertex, subleading_shower_start_point));
		three_vector subleading_shower_momentum =  utilities::scale(subleading_shower_dir, subleading_shower_ke);

		/**
		 * Information derived from leading and subleading photons
		 */
		showers_costheta = utilities::dot_product(leading_shower_dir,subleading_shower_dir);
		mass = sqrt(2*leading_shower_ke*subleading_shower_ke*(1-showers_costheta));
		energy = leading_shower_ke + subleading_shower_ke;
		momentum = utilities::add(leading_shower_momentum, subleading_shower_momentum);
		beam_costheta = utilities::dot_product(utilities::normalize(momentum),beamdir);

		/**
		 * Store output.
		 */
		s.vertex_x = std::get<0>(vertex);
		s.vertex_y = std::get<1>(vertex);
		s.vertex_z = std::get<2>(vertex);
		s.momentum = utilities::magnitude(momentum)/1000.0;
		s.beam_costheta = beam_costheta;
		s.leading_shower_index = leading_shower_index;
		s.leading_shower_ke = leading_shower_ke/1000.0;
		s.leading_shower_conv_dist = leading_shower_conv_dist;
		s.subleading_shower_index = subleading_shower_index;
		s.subleading_shower_ke = subleading_shower_ke/1000.0;
		s.subleading_shower_conv_dist = subleading_shower_conv_dist;
		s.showers_costheta = showers_costheta;
		s.mass = mass;
		s.energy = energy/1000.0;
	    }
	    else
	    {
	        s.vertex_x = -5;
		s.vertex_y = -5;
		s.vertex_z = -5;
	        s.momentum = -5;
		s.beam_costheta = -5;
		s.leading_shower_index = kNoMatch;
		s.leading_shower_ke = -5;
		s.leading_shower_conv_dist = -5;
		s.subleading_shower_index = kNoMatch;
		s.subleading_shower_ke = -5;
		s.subleading_shower_conv_dist = -5;
		s.showers_costheta = -5;
		s.mass = -5;
		s.energy = -5;
	    }
	}

	return s;
    }
}
#endif // UTILITIES_H
