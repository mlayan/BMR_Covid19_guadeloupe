#ifndef COLONIZATION_SIMULATION_RATE__H
#define COLONIZATION_SIMULATION_RATE__H

#include <Rcpp.h>
using namespace Rcpp;

RcppExport IntegerVector allSusceptibles(
	IntegerVector ids,
	IntegerVector admissions,
	IntegerVector colonizationDates,
	IntegerVector transmissionDates,
	int d
	);

RcppExport IntegerVector allPositives(
	IntegerVector ids, 
	IntegerVector admissions, 
	CharacterVector colType,
	int d
	);

RcppExport IntegerVector allSusceptiblesInSubsector(
	IntegerVector subsectors,
	IntegerVector ids,
	IntegerVector admissions,
	IntegerVector discharges,
	IntegerVector colonizationDates,
	int d, 
	int subsector
	);

RcppExport double infectionRate(
	int d,    
	int period,
	int subsector,
	IntegerVector subsectors_unique,
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
	int display=0
    );

#endif
