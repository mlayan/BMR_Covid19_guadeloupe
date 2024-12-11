#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>
#include <map>
#include "Period.h"
#include "Episode.h"
#include "utils.h"

using namespace std;

//----------------------------------------------------------------------
// Period class - Methods
//----------------------------------------------------------------------
// Constructor
Period::Period() : m_size(0), m_nDays(0), m_periodID(0) {

	m_episodes.resize(0);
	m_subsectors.clear();
	m_idsPerSubsector.resize(0);
	m_colonizedPerSubsector.resize(0);
	m_colonizedPrev.resize(0);
	m_patientPerSubsector.resize(0);
	//m_patientsTot.resize(0);
	m_covidPrev.resize(0);
	m_nbeds.resize(0);
	m_positives.resize(0);
	m_negatives.resize(0);
	m_acquisitions.resize(0);
	m_subsectorsInInteraction.resize(0);
}


// Method to create a new Period object in a loop
void Period::newPeriod() {
	m_size = 0;
	m_nDays = 0;  
	m_periodID = 0;

	m_episodes.clear();
	m_subsectors.clear();
	m_idsPerSubsector.clear();
	m_colonizedPerSubsector.clear();
	m_colonizedPrev.clear();
	m_patientPerSubsector.clear();
	//m_patientsTot.clear();
	m_covidPrev.clear();
	m_nbeds.clear();
	m_positives.clear();
	m_negatives.clear();
	m_acquisitions.clear();
	m_subsectorsInInteraction.clear();
}


// Add one Episode object to a Period object 
void Period::addEpisode(Episode episode, int period) {

	if (m_size == 0) {
		m_periodID = period;
	} else {
		if (m_periodID != period) cout << "Error in the loading of period" << endl;
	}

	if (m_nDays <= episode.discharge()) m_nDays = episode.discharge()+1;

	m_size += 1; 
	m_episodes.push_back(episode);

	if (episode.type() == "acq") m_acquisitions.push_back(episode.id());
	if (episode.type() == "pos") m_positives.push_back(episode.id());
	if (episode.type() == "neg") m_negatives.push_back(episode.id());
}


// Add prevalence matrix
void Period::addCovidPrev(std::vector<std::vector<double>> prevalenceMatrix, int display) {
	m_covidPrev = prevalenceMatrix;

	if (display == 1) {
		cout << endl << "Covid prevalence:" << endl;
		for (size_t s=0; s<prevalenceMatrix.size(); s++) {
			for (size_t d=0; d<prevalenceMatrix[s].size(); d++) {
				cout << prevalenceMatrix[s][d] << " ";
			}
			cout << endl;
		}
	}
}

void Period::addNBeds(std::vector<std::vector<double>> demographicsMatrix, int display) {
	m_nbeds = demographicsMatrix;

	if (display == 1) {
	cout << endl << "Beds occupied:" << endl;
	for (size_t s=0; s<demographicsMatrix.size(); s++) {
		for (size_t d=0; d<demographicsMatrix[s].size(); d++) {
			cout << demographicsMatrix[s][d] << " ";
		}
			cout << endl;
		}
	}
}

// Initialize colonization times
void Period::initialColonizationTimes(std::mt19937_64& gen, int display) {

	if (display == 1) {
		cout << "\nInitial augmented colonization times for period " << m_periodID+1 << ":" << endl;
		cout << "Number of acquisition episodes: " << m_acquisitions.size() << endl;
	}

	for (const auto &index : m_acquisitions) {
		int newColTime = m_episodes[index].newColTime(gen);
		m_episodes[index].setColTime(newColTime);

		if (display ==1) 
			cout << index << " " << m_episodes[index].lastNegativeTest() << " " << m_episodes[index].firstPositiveTest() << " " << newColTime << " " << m_episodes[index].colTime() << endl;
	}
}


// Initial computation of prevalences 
void Period::compute_prevalences(int display) {

	// Create m_subsectors
	std::vector<size_t> allSubsectors;
	for (size_t i=0; i<m_episodes.size(); i++) {
		if (find(allSubsectors.begin(), allSubsectors.end(), m_episodes[i].subsector()) == allSubsectors.end()) 
			allSubsectors.push_back(m_episodes[i].subsector());
	}
	sort(allSubsectors.begin(), allSubsectors.end());
	for (size_t i=0; i<allSubsectors.size(); i++) {
		m_subsectors.insert({allSubsectors[i], i});
		m_idsPerSubsector.push_back({});
	}

	if (display == 1) {
		cout << "\nInitial colonization prevalence for period " << m_periodID+1 << ":" << endl;
		for (auto i: m_subsectors) cout << i.first << " " << i.second << endl;
	}

	// Create m_idsPerSubsector
	for (size_t i=0; i<m_episodes.size();i++) {
		size_t s = m_subsectors[m_episodes[i].subsector()];
		m_idsPerSubsector[s].push_back(m_episodes[i].id());
	}
	
	// Update m_colonizedPerSubsector and m_colonizedPrev 
	std::vector<double> colonized(m_nDays);
	std::vector<double> allpatients(m_nDays);
	for (size_t d=0; d<m_nDays; d++) {
		m_colonizedPrev.push_back(0.0);
	}

	for (size_t d=0; d<m_nDays; d++) {
		double npos(0);
		double ntot(0);
		
		for (auto i : m_subsectors) {
			for (const auto &id : m_idsPerSubsector[i.second]) {
				if (m_episodes[id].getTransmissionStatus(d) == 1) npos +=1;
			}
			ntot += m_nbeds[i.second][d];
		}

		if (ntot != 0) {
			m_colonizedPrev[d] = npos/ntot;
			colonized[d] = npos;
			allpatients[d] = ntot;
		}
	}

	// Load prevalence data in matrices
	for (auto i : m_subsectors) {
		for (size_t d=0; d<m_nDays; d++) {
			double npos(0);
			//double ntot(0);
			for (const auto &id : m_idsPerSubsector[i.second]) {
				if (m_episodes[id].getTransmissionStatus(d) == 1) npos +=1;
			}
			colonized[d] = npos;
		}
		m_colonizedPerSubsector.push_back(colonized);
	}


	if (display == 1) {
		cout << endl;
		for (size_t d=0; d<m_nDays; d++) cout << m_colonizedPrev[d] << " ";
		cout << endl; 

		for (size_t i=0; i<m_colonizedPerSubsector.size();i++) {
			for (size_t j=0; j<m_colonizedPerSubsector[i].size(); j++) {
				cout << m_colonizedPerSubsector[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
}



// Update prevalences 
void Period::update_prevalences(int epi, int display) {

	// Update matrix of colonization prevalence
	size_t sub = m_subsectors[m_episodes[epi].subsector()];
	size_t start = m_episodes[epi].admissionTime();
	size_t end = m_episodes[epi].discharge();

	for (size_t d=start; d<end+1; d++) {
		m_colonizedPrev[d] = 0;
		double npos_sub(0.);

		for (const auto &index : m_idsPerSubsector[sub]) {
			if (m_episodes[index].getTransmissionStatus(int(d)) == 1) npos_sub+=1;
		}
		m_colonizedPerSubsector[sub][d] = npos_sub;

		// Update covid prevalence taking into account all subsectors
		double npos_tot(0.);
		double ntot(0.);
		for (auto i: m_subsectors) {
			npos_tot += m_colonizedPerSubsector[i.second][d];
			ntot += m_nbeds[i.second][d];
		}

		if (ntot != 0) {
			m_colonizedPrev[d] = npos_tot / ntot;
		} else {
			m_colonizedPrev[d] = 0.;
		}
	}

	if (display == 1) {
			cout << endl;
			for (size_t d=0; d<m_nDays; d++) cout << m_colonizedPrev[d] << " ";
			cout << endl; 

			for (size_t i=0; i<m_colonizedPerSubsector.size();i++) {
				for (size_t j=0; j<m_colonizedPerSubsector[i].size(); j++) {
					cout << m_colonizedPerSubsector[i][j] << " ";
				}
				cout << endl;
			}
			cout << endl;
		}

}


void Period::display_prevalence_matrix() {

	// Colonization prevalence per day (vector) 
	for (size_t d=0; d<m_nDays; d++) cout << m_colonizedPrev[d] << " ";
	cout << endl; 
	
	// Number of colonized individuals per subsector per day (matrix)
	for (size_t i=0; i<m_colonizedPerSubsector.size();i++) {
		for (size_t j=0; j<m_colonizedPerSubsector[i].size(); j++) {
			cout << m_colonizedPerSubsector[i][j] << " ";
		}
		cout << endl;
	}

	// When individuals are positives
	cout << "Positive inviduals:" << endl;
	for (const auto &i: m_positives) {
		cout << i << " " << m_subsectors[m_episodes[i].subsector()] << " " << m_episodes[i].admissionTime() << " " << m_episodes[i].discharge() << endl;
	}

	cout << "Acquisitions:" << endl;
	for (const auto &i: m_acquisitions) {
		cout << i << " " << m_subsectors[m_episodes[i].subsector()] << " " << m_episodes[i].admissionTime() << " " << m_episodes[i].colTime() << " " << m_episodes[i].discharge() << endl;
	}

	cout << endl;
}



// Propose new value for colonization time
int Period::newColTime(int index, std::mt19937_64& gen) {
	int colTime = m_episodes[index].newColTime(gen);
	return colTime;
}


// Set new value for colonization time
void Period::setColTime(int index, int colTime) {
	m_episodes[index].setColTime(colTime);
}


// Display data
void Period::display_data() {

	// Show "raw" data 
	// for (size_t i=0; i<m_size; i++) {
	// 	for (size_t d=0; d < m_episodes[i].periodLength(); d++) {
	// 		cout << m_periodID << " " << m_episodes[i].id() << " " << m_episodes[i].subsector() << " " << m_episodes[i].type() << " " << m_episodes[i].getStatus(d) << endl;
	// 	}
	// }

	// Show episode attributes  
	for (size_t i=0; i<m_size; i++) {
		cout << m_periodID << " " << m_episodes[i].id() << " " << m_episodes[i].subsector() << " " << m_episodes[i].type() << " " << m_episodes[i].admissionTime() << " " << m_episodes[i].firstPositiveTest() << " " << m_episodes[i].discharge() << endl;
	}
	cout << endl;
}


// Compute loglikelihood
double Period::compute_log_lik(
  std::vector<double> parameter, 
  std::vector<int> selectedParam, 
  int display
  ) {

  double LL = 0.0;

  // Log Likelihood of negative individuals
  for (const auto &index : m_negatives) { 
  	LL += m_episodes[index].compute_log_lik_neg(
  		parameter, 
  		selectedParam, 
  		m_periodID,
  		m_subsectors,
  		m_subsectorsInInteraction,
  		m_colonizedPrev,
  		m_covidPrev,
  		display
  		);
  }

	// Log Likelihood of individuals that acquired a resistant bacteria
  for (const auto &index : m_acquisitions) { 
  	LL += m_episodes[index].compute_log_lik_acq(
  		parameter, 
  		selectedParam, 
  		m_periodID,
  		m_subsectors,
  		m_subsectorsInInteraction,
  		m_colonizedPrev,
  		m_covidPrev,
  		display
  		);
  }

  if (display) cout << "LL total: " << LL << "\n\n";

  return LL;
}



