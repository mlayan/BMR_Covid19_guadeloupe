#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
// #include <boost/random.hpp>
#include "McmcObject.h"
#include "Period.h"
#include "utils.h"

using namespace std;

//----------------------------------------------------------------------
// MCMC class
//----------------------------------------------------------------------
// Constructor
McmcObject::McmcObject() {
    m_gen.seed(std::random_device{}());
    m_globalLogLik = 0.0;
    m_iterations = 0;
    m_iterTimeInfection = 0;
    m_acceptedMoveData = 0;
	m_proposedMoveData = 0;
    m_rateForRandomWalk.resize(0);
    m_numberOfMoveAccepted.resize(0);
    m_numberOfMoveProposed.resize(0);
    m_accept.resize(0);
    m_parameter.resize(0);
    m_selectedParam.resize(0);
    m_periodLogLik.resize(0);
    std::vector<Period> m_data;
    m_data.resize(0);
    m_nPeriod = 0;
}

McmcObject::McmcObject(
                       int nIterations,
                       std::vector<Period> data,
                       std::vector<double> parameter,
                       std::vector<int> selectedParameter,
                       std::vector<double> rateForRandomWalk,
                       int nIterTimeInfection
                       ) {
    // Default values
    m_globalLogLik= 0.0;
    m_numberOfMoveAccepted.resize(0);
    m_numberOfMoveProposed.resize(0);
    m_accept.resize(parameter.size());
    m_acceptedMoveData = 0;
	m_proposedMoveData = 0;

    // Copy input vectors
    //m_gen.seed(seed);
    m_gen.seed(std::random_device{}());
    cout << "=============MCMC===============" << endl;
    cout << "Seed: " << m_gen() << endl;

    m_periodLogLik.resize(data.size());
    m_iterations = nIterations;
    m_iterTimeInfection = nIterTimeInfection;
    m_rateForRandomWalk = rateForRandomWalk;
    m_parameter = parameter;
    m_selectedParam = selectedParameter;
    m_data = data;
    m_nPeriod = data.size();
}

// Initialize parameter values
void McmcObject::initial_param_values(std::vector<std::string> paramNames, int display) {
    
    for (size_t parID=0; parID < m_parameter.size(); parID++) {
        if ( m_selectedParam[parID] == 1 ) {
            m_parameter[parID] = rnorm(m_gen, 0.0, 10.0);
            if (display == 1) 
                cout << paramNames[parID] << ": " << m_parameter[parID] << endl;
        }
    }
}

// Modification of attributes
void McmcObject::resetMoves() {
    m_numberOfMoveProposed.clear();
    m_numberOfMoveProposed.resize(m_parameter.size());
    m_numberOfMoveAccepted.clear();
    m_numberOfMoveAccepted.resize(m_parameter.size());
    m_acceptedMoveData = 0;
	m_proposedMoveData = 0;
}

// Initialize infection time
void McmcObject::initialize_augmented_data(std::string data_augmentation, int display) {
    
    size_t period;
    for (period=0; period<m_nPeriod; period++) {

        // Initialize infection time
        if (data_augmentation == "yes") m_data[period].initialColonizationTimes(m_gen, display);

        // Compute prevalences
        m_data[period].compute_prevalences(display);
    } 
}

// Initialize infection time
void McmcObject::initial_log_lik(int display) {

    m_globalLogLik = 0;
    
    for (size_t period=0; period < m_nPeriod; period++) {
        m_periodLogLik[period] = m_data[period].compute_log_lik(m_parameter, m_selectedParam, display);
        m_globalLogLik += m_periodLogLik[period];
        //m_data[period].display_prevalence_matrix();
    }
}


// Update parameter value in the mcmc chain
void McmcObject::update_parameter(int parID, double step) {

    // Random walk
	double oldValue = m_parameter[parID];
	double newValue = oldValue + rnorm(m_gen, 0.0, m_rateForRandomWalk[parID]);
    // cout << "Step " << step << ": " << parID << " - " << oldValue << " - " << newValue << endl; 

    // Log ratio of priors 
    double logRatioPrior = 0.5 * ( pow(oldValue, 2) - pow(newValue, 2) ) / pow(10.0,2);
    //double logRatioProposal = 0.0;

    // cout << logRatioPrior << endl ;

	// Compute new log likelihood
    m_parameter[parID] = newValue;
    std::vector<double> newLogLikPeriod(m_nPeriod, 0.0);
    double newLogLikGlobal(0.0);

    for (size_t period=0; period < m_nPeriod; period++) {
        int display = 0;

        if ( 
            (parID == 1 && period != 1) || 
            (parID == 2 && period != 2) || (parID == 5 && period != 2) ||
            (parID == 3 && period != 3) || (parID == 6 && period != 3) ||
            (parID == 7 && period == 0) || (parID == 0 && period != 0)
             ) {
            // No update of the likelihood when the parameter is period specific but does not match the period
            newLogLikPeriod[period] = m_periodLogLik[period];
            
        } else {
            newLogLikPeriod[period] = m_data[period].compute_log_lik(m_parameter, m_selectedParam, display);
        }     
        newLogLikGlobal += newLogLikPeriod[period];

/*        cout << endl << "Step " << step << " - Param " << parID << endl; 
        m_data[period].display_prevalence_matrix();*/
        if (isinf(newLogLikPeriod[period])) {
            cout << step << " " << parID << " " << period << " " << newLogLikPeriod[period] << " " << m_periodLogLik[period] << endl;
            double test = m_data[period].compute_log_lik(m_parameter, m_selectedParam, 1);
        }

    }

	// Update log likelihood
    // newLogLikGlobal = 0.0;
    // m_globalLogLik = 0.0;
	double Q = newLogLikGlobal - m_globalLogLik + logRatioPrior; // + logRatioProposal
    double random_number = log(runif(m_gen));
    // cout << random_number << " " << Q << endl;

	if (  random_number < Q )
	{
		m_globalLogLik = newLogLikGlobal;
		m_periodLogLik = newLogLikPeriod;
        m_numberOfMoveAccepted[parID]++;
        m_accept[parID] = (1 + (step - 1) * m_accept[parID]) / step;
	}
	else {
		m_parameter[parID] = oldValue;
        m_accept[parID] = (0 + (step - 1) * m_accept[parID]) / step;
	}
	
    m_numberOfMoveProposed[parID]++;
}



// Update augmented infection times
void McmcObject::update_augmented_acquisitions_times() {

    // Loop over houses
    size_t period;
    for (period = 0; period<m_data.size(); ++period) {
        std::vector<size_t> acquisitions = m_data[period].getAcquisitions();

        // Loop over episodes of acquisitions
        for (const auto &acq : acquisitions) {

            int display = 0;

            // Independent sampler from the incubation period
            double oldValue = m_data[period].getColTime(acq);
            double newValue = m_data[period].newColTime(acq, m_gen);

            // Update infection time and person-to-person transmission rates
            m_data[period].setColTime(acq, newValue);
            m_data[period].update_prevalences(acq, display);   
            
/*            if (period == 1) {
                cout << endl << "Period " << period << " - Ind " << acq << endl; 
                cout << oldValue << " " << newValue << " " << m_data[period].getSubsector(acq) << " " << m_data[period].getDischarge(acq) << endl;
                m_data[period].display_prevalence_matrix();
            }*/
            
            // Log likelihood
            display = 0;
            double newLogLikOfPeriod = m_data[period].compute_log_lik(m_parameter, m_selectedParam, display);
            double currentLogLik = m_periodLogLik[period];
            double differenceLogLik = newLogLikOfPeriod - currentLogLik;

            if( log(runif(m_gen)) < differenceLogLik ) {
                m_globalLogLik += differenceLogLik;
                m_periodLogLik[period] = newLogLikOfPeriod;
                m_acceptedMoveData++;

            } else {                
                m_data[period].setColTime(acq, oldValue);
                m_data[period].update_prevalences(acq, display);
                
            }

            m_proposedMoveData++;
        }
    }
}
