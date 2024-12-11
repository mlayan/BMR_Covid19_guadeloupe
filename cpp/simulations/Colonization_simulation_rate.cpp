#include <Rcpp.h>
#include "Colonization_simulation_rate.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]


//////////////////////////////////////////////
// [[Rcpp::export]]
IntegerVector allSusceptibles(
	IntegerVector ids,
	IntegerVector admissions,
	IntegerVector colonizationDates,
	IntegerVector transmissionDates,
	int d
	) {

	IntegerVector out;

	for (const auto &ind : ids) {
		if (admissions[ind] == d) {
			out.push_back(ind);
		}
	}

	return out;
}

//////////////////////////////////////////////
// [[Rcpp::export]]
IntegerVector allPositives(
	IntegerVector ids, 
	IntegerVector admissions, 
	CharacterVector colType,
	int d
	) {

	IntegerVector out;

	for (const auto &ind : ids) {
		if (admissions[ind] == 0 && colType[ind] == "pos") {
			out.push_back(ind);
		}
	}

	return out;

}


//////////////////////////////////////////////
// [[Rcpp::export]]
IntegerVector allSusceptiblesInSubsector(
	IntegerVector subsectors,
	IntegerVector ids,
	IntegerVector admissions,
	IntegerVector discharges,
	IntegerVector colonizationDates,
	int d, 
	int subsector
	) {


	IntegerVector out;

	for (const auto &ind : ids) {
		if (subsectors[ind] == subsector) {
			if (colonizationDates[ind] < 0 && admissions[ind] <= d && d <= discharges[ind]) {
				out.push_back(ind);
			}
		}
	}

	return out;
}

//////////////////////////////////////////////
// [[Rcpp::export]]
double infectionRate(
	int d,    
	int period,
	int subsector,
	IntegerVector subsectors_unique, // To use when using information on sectors that are in contact
	NumericMatrix pcovid,
	NumericMatrix nBeds, 
	IntegerMatrix nTransmitters, 
	double intercept,
	double intercept2,
	double intercept3,
	double intercept4,
	double pcov2,
	double pcov3,
	double pcov4,
	double pcov,
	int display
  ) {
	
	// Lambda 
	double log_lambda_j = intercept + pcov * pcovid( subsector, d );
	if (period == 1) log_lambda_j += intercept2 + pcov2 * pcovid( subsector, d ); 
	if (period == 2) log_lambda_j += intercept3 + pcov3 * pcovid( subsector, d ); 
	if (period == 3) log_lambda_j += intercept4 + pcov4 * pcovid( subsector, d ); 

	// Number of colonized patients in subsector j and in all subsectors in contact with j
	// and their total number of patients 
	double colonized_j_and_others = 0.0;
	double total_j_and_others = 0.0;
	for (int i=0; i<subsectors_unique.size(); i++) {
		colonized_j_and_others += nTransmitters(i, d); 
		total_j_and_others += nBeds(i,d);
	}

	// Probability of infection at day d
	double out = 1 - exp( - exp(log_lambda_j) * colonized_j_and_others  / total_j_and_others );

	if (display == 1) {
		Rcout << "Subsector " << subsector << endl;
		Rcout << intercept << " " << intercept2 << " " << intercept3 << " " << intercept4 << " " << pcov2 << " " << pcov3 << " " << pcov4 << " " << pcov << endl;
		Rcout << log_lambda_j << endl;
		Rcout << pcovid( subsector, d ) << endl;
		Rcout << nTransmitters(subsector, d) << endl;
		Rcout << colonized_j_and_others << " " << total_j_and_others << endl;
	}

	return out;
  
}
