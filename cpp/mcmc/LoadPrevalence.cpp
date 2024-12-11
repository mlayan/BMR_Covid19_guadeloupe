#include<vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "LoadPrevalence.h"

using namespace std;

std::vector<std::vector<std::vector<double>>> getPrevalenceMatrix(
    std::string fileName,
    int numberOfPeriodsColonization,
    int display
    ) {

    std::vector<std::vector<std::vector<double>>> output;

	std::ifstream infile(fileName.c_str());
	if(infile)
    {
         // Variables to store data
        int period(-1), lastperiod(-1), lastsubsector(-1), subsector(-1), t(-1);
        double val(-1.0);
        size_t subsector_index = 0;

        // Variables to verify data loading
        int numberOfPeriods=0;
        int numberOfSubsectors=0;

        // std::vector<std::vector<double>> currMatrix(10, std::vector<double>(0));
        std::vector<std::vector<double>> currMatrix;

        for ( std::string line; std::getline(infile, line); )
        {
            // Create a stringstream of the current line
            std::istringstream in(line);
            
            // Store information in variables 
            in >> period >> subsector >> t >> val;

            // Create prevalence matrix
            if (lastperiod != period && numberOfPeriods == 0) // First period first subsector
            {

                std::vector<double> newSubsector0 = {val};
                currMatrix.push_back(newSubsector0);
                numberOfPeriods++;
                numberOfSubsectors++;
                lastperiod = period;
                lastsubsector = subsector;

            } else if (lastperiod != period && numberOfPeriods > 0) {
                output.push_back(currMatrix);
                if (display == 1) {
                    cout << "Prevalence matrix for period " << numberOfPeriods << ":" << endl; 
                    for (size_t i=0; i < currMatrix.size(); i++) {
                        for (size_t j=0; j < currMatrix[i].size(); j++) {
                            cout << numberOfPeriods-1 << " " << i << " " << j << " " << currMatrix[i][j] << endl; 
                        }
                    }
                }
                
                subsector_index = 0;
                numberOfPeriods++;
                lastperiod = period;
                lastsubsector = subsector;
                std::vector<double> newSubsector1 = {val};
                currMatrix.resize(0);
                currMatrix.push_back(newSubsector1);

            } else if (lastperiod == period && lastsubsector != subsector) {

                std::vector<double> newSubsector2 = {val};
                currMatrix.push_back(newSubsector2);
                if (period == 0) numberOfSubsectors++;
                lastsubsector = subsector;
                subsector_index++;


            } else { // Add daily prevalence for subsector during period
                currMatrix[subsector_index].push_back(val);
            }
        }


        // Add last period-specific prevalence matrix 
        output.push_back(currMatrix);
        if (display == 1) {
            cout << "Prevalence matrix for period " << numberOfPeriods << ":" << endl; 
            for (size_t i=0; i < currMatrix.size(); i++) {
                for (size_t j=0; j < currMatrix[i].size(); j++) {
                    cout << numberOfPeriods-1 << " " << i << " " << j << " " << currMatrix[i][j] << endl; 
                }
            }
        }

        if (numberOfPeriods != numberOfPeriodsColonization) {
            cout << "ERROR: The number of periods is not in agreement between the colonization data and the prevalence data" << endl;
            cout << numberOfPeriodsColonization << " " << numberOfPeriods << endl;
        }

    }

    return(output);
}
