#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>
#include "Episode.h"
#include "utils.h"

using namespace std;

//----------------------------------------------------------------------
// Episode class - Methods
//----------------------------------------------------------------------
// Constructor
Episode::Episode() {
  m_id = 0; 
  m_subsector = 0;
  m_admissionTime = -1;
  m_discharge = -1;
  m_firstPositiveTest = -1;
  m_colTime = -1;
  m_type = "";
}

// Constructor during data loading
void Episode::newEpisode(size_t id, size_t subsector, int admission, int discharge, int detectionDay, 
  int lastNegative, std::string type, std::string data_augmentaton) {
  m_id = id;
  m_subsector = subsector;
  m_admissionTime=admission;
  m_discharge=discharge;
  if (data_augmentaton == "yes") {
    m_firstPositiveTest=detectionDay;
    m_lastNegativeTest=lastNegative;
  } else {
    m_colTime = detectionDay-1;
  }
  m_type = type;
} 

// Accessors
int Episode::getTransmissionStatus(int day) {
  int output = 0;

  if (day >= m_admissionTime && day <= m_discharge) {
    if (m_type == "pos") output = 1;

    // Individuals that are colonized during their stay cannot transmit the day they are colonized
    if (m_type == "acq" && day > m_colTime) output = 1;
  }

  return output;
}


int Episode::getStay(int day) {
  int output = 0;
  if (m_admissionTime <= day && day <= m_discharge) {
    output = 1;
  }
  return output;
}


// Data augmentation
int Episode::newColTime(std::mt19937_64& gen) {
  int output = m_lastNegativeTest + rdiscrete(gen, m_firstPositiveTest-m_lastNegativeTest);
  return output;
}

void Episode::setColTime(int colTime) {
  m_colTime = colTime;
}

// Compute log likelihood for negative individuals
double Episode::compute_log_lik_neg(
  std::vector<double> parameters, 
  std::vector<int> selectedParam,
  int period,
  std::map<size_t, size_t> all_subsectors,
  std::vector<size_t> subsectorsInInteraction,
  std::vector<double> colonizedPrev,
  std::vector<std::vector<double>> covidPrev,
  int display
  ) {

  double LL=0.0;
  double loglambda = 0.0;
  size_t subsector = all_subsectors[m_subsector];

  for (int day=m_admissionTime; day<m_discharge+1; day++) {
    // npos = 0;
    // ntot = 0;

    loglambda = parameters[4] + parameters[7] * covidPrev[subsector][day];
    if (period == 0) loglambda += parameters[0];
    if (period == 1) loglambda += parameters[1] * covidPrev[subsector][day];
    if (period == 2) loglambda += parameters[2] * covidPrev[subsector][day] + parameters[5];
    if (period == 3) loglambda += parameters[3] * covidPrev[subsector][day] + parameters[6];

    // for (const auto &subsector : m_subsectorsInInteraction) {
    //   npos += colonizedPerSubsector[subsector][day];
    //   ntot += patientPerSubsector[subsector][day];
    // }

    LL += -exp(loglambda) * colonizedPrev[day];
    if (display == 1) {
      cout << "LL ind " << m_id << " - day " << day << " - type " << m_type << ": " << loglambda << " " << exp(loglambda) << " " << colonizedPrev[day] << endl; 
      if (isinf(exp(loglambda))) cout << m_subsector << " " << subsector << " " << covidPrev.size() << " " << covidPrev[subsector].size() << " " << day << " " << covidPrev[subsector][day] << endl;
    }


  }
  if (display == 1) cout << "LL ind " << m_id << " - total: " << LL << endl; 


  return LL;
}



// Compute log likelihood for individuals that acquired a resistant bacteria
double Episode::compute_log_lik_acq(
  std::vector<double> parameters, 
  std::vector<int> selectedParam,
  int period,
  std::map<size_t, size_t> all_subsectors,
  std::vector<size_t> subsectorsInInteraction,
  std::vector<double> colonizedPrev,
  std::vector<std::vector<double>> covidPrev,
  int display
  ) {

  double LL=0.0;
  double loglambda = 0.0;
  size_t subsector = all_subsectors[m_subsector];
  
  // Period of survival
  for (int day=m_admissionTime; day<m_colTime; day++) {
    // npos = 0;
    // ntot = 0;

    loglambda = parameters[4] + parameters[7] * covidPrev[subsector][day];
    if (period == 0) loglambda += parameters[0];
    if (period == 1) loglambda += parameters[1] * covidPrev[subsector][day];
    if (period == 2) loglambda += parameters[2] * covidPrev[subsector][day] + parameters[5];
    if (period == 3) loglambda += parameters[3] * covidPrev[subsector][day] + parameters[6];

    // for (const auto &subsector : m_subsectorsInInteraction) {
    //   npos += colonizedPerSubsector[subsector][day];
    //   ntot += patientPerSubsector[subsector][day];
    // }
    
    if (display == 1) {
      cout << "LL ind " << m_id << " - day " << day << " - type " << m_type << ": " << loglambda << " " << exp(loglambda) << " " << colonizedPrev[day] << endl; 
    if (isinf(exp(loglambda))) cout << m_subsector << " " << subsector << " " << covidPrev.size() << " " << covidPrev[subsector].size() << " " << day << " " << covidPrev[subsector][m_colTime] << endl;
    }
    LL += -exp(loglambda) * colonizedPrev[day];
  }

  // Colonization time
  // npos = 0;
  // ntot = 0;

  loglambda = parameters[4] + parameters[7] * covidPrev[subsector][m_colTime];
  if (period == 0) loglambda += parameters[0];
  if (period == 1) loglambda += parameters[1] * covidPrev[subsector][m_colTime];
  if (period == 2) loglambda += parameters[2] * covidPrev[subsector][m_colTime] + parameters[5];
  if (period == 3) loglambda += parameters[3] * covidPrev[subsector][m_colTime] + parameters[6];

  // for (const auto &subsector : m_subsectorsInInteraction) {
  //   npos += colonizedPerSubsector[subsector][m_colTime];
  //   ntot += patientPerSubsector[subsector][m_colTime];
  // }
  
  LL += log( 1 - exp(-exp(loglambda) * colonizedPrev[m_colTime]) );

  if (display == 1) {
    cout << "LL ind " << m_id << " - day " << m_colTime << " - type " << m_type << ": " << loglambda << " " << exp(loglambda) << " " << colonizedPrev[m_colTime] << endl; 
    cout << "LL ind " << m_id << " - total: " << LL << endl; 
    if (isinf(exp(loglambda))) cout << m_subsector << " " << subsector << " " << covidPrev.size() << " " << covidPrev[subsector].size() << " " << m_colTime << " " << covidPrev[subsector][m_colTime] << endl;
  }

  return LL;
}