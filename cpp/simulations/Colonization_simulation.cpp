 #include <Rcpp.h>
#include "Colonization_simulation_rate.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
DataFrame simulateCol(
	DataFrame occupancy,
  NumericMatrix covid_data,
  NumericMatrix nBeds,
  int period,
  double intercept,
  double intercept2,
  double intercept3,
  double intercept4,
  double pcov2,
  double pcov3,
  double pcov4,
  double pcov
  )
  {

    int display = 0;


    // Display input parameters 
    if (display == 1) {
      Rcout << "###############################" << endl;
      Rcout << "### Start of the simulation ###" << endl;
      Rcout << "###############################" << "\n\n";

      Rcout << "-------- Parameter values --------" << endl;
      Rcout << "Baseline transmission rate: " << intercept << endl; 
      Rcout << "Period 2: " << intercept2 << endl;
      Rcout << "Period 3: " << intercept3 << endl;
      Rcout << "Period 3: " << intercept4 << endl;
      Rcout << "COVID-19 (period 2): " << pcov2 << endl; 
      Rcout << "COVID-19 (period 3): " << pcov3 << endl; 
      Rcout << "COVID-19 (period 4): " << pcov4 << endl; 
      Rcout << "COVID-19 (all periods): " << pcov << "\n\n"; 

    }


    // Initialize transmission vectors
    int nDays = nBeds.ncol();
    int nStays = occupancy.nrow();

    IntegerVector subsectors = occupancy["subsector"];
    IntegerVector subsectors_unique = sort_unique(subsectors);
    int nSubsectors = subsectors_unique.size();
    IntegerVector subsectors_recoded(nStays);
    for (int i=0; i<nStays; i++) {
      for (int s=0; s<nSubsectors; s++) {
        if (subsectors[i] == subsectors_unique[s]) subsectors_recoded[i] = s;
      }
    }

    IntegerVector ids = occupancy["id"];
    IntegerVector admissions = occupancy["admission"];
    IntegerVector discharges = occupancy["discharge"];

    IntegerVector colonizationDates (nStays, -1);
    IntegerVector transmissionDates (nStays, -1);
    IntegerMatrix nTransmitters(nSubsectors, nDays);
    IntegerVector allColonized;

    CharacterVector colType(nStays, "neg");


    // Number of individuals already colonized at admission
    // We assume that 6% of the individuals are already colonized at admission
    double p_positive;
    if (period == 0) p_positive = 0.11;
    if (period == 1) p_positive = 0.04;
    if (period == 2) p_positive = 0.12;
    int nPositives = max(int(rbinom(1, nStays, p_positive)[0]), 1);
    IntegerVector positives_ids = sample(ids, nPositives, false);
    for (const auto &p: positives_ids) {
      colType[p] = "pos";
      allColonized.push_back(p);
      colonizationDates[p] = admissions[p];
      transmissionDates[p] = admissions[p];

      for (int d=admissions[p]; d<discharges[p]+1; d++) {
        nTransmitters(subsectors_recoded[p], d)++;
      }
      if (display == 1) Rcout << subsectors[p] << " ";
    }

    // Print information
    IntegerVector positives_0 = allPositives(ids, admissions, colType, 0);
  
    if (display == 1) {
      Rcout << "-------- Initialize positive status --------" << endl;
      Rcout << "N positive patients at admission: " << nPositives << endl;
      Rcout << "N positive patients at time 0: " << positives_0.size() << endl;
      Rcout << "Ids of positive patients at time 0: "; 
      for (const auto &i: positives_0) Rcout << i << " ";
      Rcout << endl << "Subsectors of colonized individuals: " ;
      for (const auto &i: positives_0) Rcout << subsectors[i] << " ";

      Rcout << endl << "Update of colonization matrix:\n";
      for (int s=0; s<nSubsectors; s++) {
        for (int d=0; d<nDays; d++) Rcout << nTransmitters(s,d) << " ";
        Rcout << "\n";
      }
      Rcout << "\n";
    }

    // Simulate Epidemic within ICU ward
    // Simulations stop at nDays-1 because if a patient is infected 
    // at nDays, it will be tested positive only at nDays+1
    if (display==1) Rcout << "-------- Transmission of colonization --------" << endl;

    for (int curr_time = 0; curr_time < nDays; ++curr_time) { 

      if (display == 1) Rcout << "-------- Time " << curr_time << " --------" << endl;

      for (int s=0; s<nSubsectors; s++) {
        double pInfection = infectionRate(
          curr_time,
          period,
          s,
          subsectors_unique, 
          covid_data,
          nBeds, 
          nTransmitters, 
          intercept,
          intercept2,
          intercept3,
          intercept4,
          pcov2,
          pcov3,
          pcov4,
          pcov,
          display 
        ); 

        if (display == 1) Rcout << "Subsector " << subsectors_unique[s] << ": " << pInfection << endl;

        IntegerVector susceptibles = allSusceptiblesInSubsector(
          subsectors, 
          ids, 
          admissions, 
          discharges, 
          colonizationDates, 
          curr_time, 
          subsectors_unique[s]
        );

        for (const auto &ind : susceptibles) {

          double rn = runif(1)[0];
          if (display == 1) Rcout << ind << " " << rn << " " << pInfection << endl; 

          if (rn < pInfection) {

            // Rcout << ind << " " << subsectors[ind] << endl;
            allColonized.push_back(ind);
            colonizationDates[ind] = curr_time;
            transmissionDates[ind] = curr_time+1;
            colType[ind] = "acq";

            int lastTime = min(discharges[ind]+1, nDays-1);
            for (int d=curr_time+1; d<lastTime; d++) {
              nTransmitters(s, d) ++;
            }

          }
        }
      }

      // Rcout << endl;

    }

    if (display == 1) {
      Rcout << "-------- End of simulation --------" << "\n\n";
      Rcout << "Ids of all colonized individuals:" << endl;
      for (int i=0; i<allColonized.size(); i++) Rcout << allColonized[i] << " ";  
    }
    
    // Output
    return DataFrame::create(
    	_("subsector") = occupancy["subsector"],
    	_("id") = occupancy["id"],
    	_("admission") = occupancy["admission"],
      _("discharge") = occupancy["discharge"],
      _("detection_date") = transmissionDates,
      _("colType") = colType   
    	);
}




