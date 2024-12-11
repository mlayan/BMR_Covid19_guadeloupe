#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <random>
#include <chrono>
#include <vector>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <experimental/filesystem>
#include "LoadPrevalence.h"
#include "McmcObject.h"
#include "Period.h"
#include "Episode.h"

using namespace std;
namespace fs = std::experimental::filesystem;

//----------------------------------------------------------------------
// Load data
//----------------------------------------------------------------------
std::vector<Period> buildData(
    std::string episodeDataFile, 
    std::string covidPrevalenceDataFile, 
    std::string nbedsDataFile,
    std::string data_augmentation
    )
{


    int display=0;

    // Variable with all the data
    std::vector<Period> output;
    output.resize(0);

    // Read Episode data file 
    /*
    The data file should be a space-separated table
    */
    std::ifstream infile(episodeDataFile.c_str());
    if(infile)
    {
        cout << "=============DATA===============" << endl;
        cout << "Read data file: " << episodeDataFile << endl;

        // Variables to store data
        int lastperiod(-1), admission(0), discharge(0), period(0), id(0), detectionDay(0), lastNegative(0);
        size_t subsector(0);
        std::string colType;
        std::vector<int> t_col(0);

        // Variables to verify data loading
        int numberOfPeriods(0);
        int numberOfEpisodes(0);
        int max_discharge(0);
        std::vector<int> numberOfDaysPerPeriod(0);
        std::vector<int> numberOfEpisodesPerPeriod(0);


        Episode currEpisode;
        Period currPeriod;

        for ( std::string line; std::getline(infile, line); )
        {
            // Create a stringstream of the current line
            std::istringstream in(line);
            
            // Store information in variables 
            if (data_augmentation == "yes") in >> period >> id >> subsector >> admission >> discharge >> detectionDay >> lastNegative >> colType;
            if (data_augmentation == "no") in >> period >> id >> subsector >> admission >> discharge >> detectionDay >> colType;

            // Create Episodes and Periods
            if (numberOfPeriods == 0)
                {
                    numberOfPeriods++;
                    numberOfEpisodes++;
                    lastperiod = period;
                    currEpisode.newEpisode(id, subsector, admission, discharge, detectionDay, lastNegative, colType, data_augmentation);
                    max_discharge = discharge;

                } else if (numberOfPeriods > 0 && lastperiod != period) { // Move from one period to another
                    
                    // Initialize the first episode of the first period
                    numberOfEpisodesPerPeriod.push_back(numberOfEpisodes);
                    numberOfDaysPerPeriod.push_back(max_discharge);
                    currPeriod.addEpisode(currEpisode, lastperiod);
                    
                    output.push_back(currPeriod);

                    numberOfPeriods++;
                    currPeriod.newPeriod();
                    lastperiod = period;
                    numberOfEpisodes=1;
                    max_discharge = discharge;
                    currEpisode.newEpisode(id, subsector, admission, discharge, detectionDay, lastNegative, colType, data_augmentation);

                } else { // Move from one episode to another in the same period

                    numberOfEpisodes++;
                    currPeriod.addEpisode(currEpisode, period);
                    currEpisode.newEpisode(id, subsector, admission, discharge, detectionDay, lastNegative, colType, data_augmentation);
                    if (discharge > max_discharge) max_discharge = discharge;

                } 
        }

        // Add last household
        currPeriod.addEpisode(currEpisode, lastperiod);
        output.push_back(currPeriod);
        numberOfDaysPerPeriod.push_back(max_discharge);
        numberOfEpisodesPerPeriod.push_back(numberOfEpisodes);

        // General information
        cout << "Number of periods: " << numberOfPeriods << endl;
        cout << "Output size: " << output.size() << "\n\n";

        cout << "Number of episodes per period:" << endl;
        for (size_t i = 0; i < numberOfEpisodesPerPeriod.size(); i++)
            cout << "Period " << i+1 << ": " << numberOfEpisodesPerPeriod[i] << " " << output[i].getSize() << endl;

        cout << endl << "Number of days per period:" << endl;
        for (size_t i = 0; i < numberOfEpisodesPerPeriod.size(); i++)
            cout << "Period " << i+1 << ": " << numberOfDaysPerPeriod[i]+1 << endl;

        cout << endl << "Acquisitions IDs per period: " << endl;
        for (size_t i = 0; i < numberOfDaysPerPeriod.size(); i++) {
            cout << "Period " << i+1 << ": ";
            for (size_t j = 0; j<output[i].getAcquisitions().size(); j++) cout << output[i].getAcquisitions()[j] << " ";
            cout << endl;
        }
        
        cout << endl << "Positive IDs per period: " << endl;
        for (size_t i = 0; i < numberOfDaysPerPeriod.size(); i++) {
            cout << "Period " << i+1 << ": ";
            for (size_t j = 0; j<output[i].getPositives().size(); j++) cout << output[i].getPositives()[j] << " ";
            cout << endl;
        }

        cout << endl << "Negative IDs per period: " << endl;
        for (size_t i = 0; i < numberOfDaysPerPeriod.size(); i++) {
            cout << "Period " << i+1 << ": ";
            for (size_t j = 0; j<output[i].getNegatives().size(); j++) cout << output[i].getNegatives()[j] << " ";
            cout << endl;
        }
        cout << endl;

        // Display data
        // Pay attention here, when the "raw" data are printed, "acq" episodes are transformed
        // into "pos" episodes
        if (display ==1) {
            for (size_t i=0; i<output.size(); i++) output[i].display_data();
        }

        // Add prevalence data
        std::vector<std::vector<std::vector<double>>> covid = getPrevalenceMatrix(covidPrevalenceDataFile, numberOfPeriods, display);
        cout << "Covid prevalence matrix loaded" << endl;
        std::vector<std::vector<std::vector<double>>> nbeds = getPrevalenceMatrix(nbedsDataFile, numberOfPeriods, display);
        cout << "Number of beds occupied per day loaded" << "\n\n";
        
        for (size_t i =0; i<output.size(); i++) {
            output[i].addCovidPrev(covid[i], 1);
            output[i].addNBeds(nbeds[i], 1);
        }

        return(output);

    } else {
        cout << "ERROR: Cannot open the file." << endl;
        return(output);
    }
}


//----------------------------------------------------------------------
// Run mcmc
//----------------------------------------------------------------------
void runMCMC(
    McmcObject mcmc,
    std::string outputFile,
    std::string augmentedDataFile,
    int pas,
    std::vector<int> idOfSelectedParameter,
    std::vector<std::string> paramNames,
    std::string data_augmentation
)
{
	ofstream output(outputFile.c_str());

	int iteration, iter, parID;
	int numberOfIteration = int(mcmc.iteration() / pas);
	int nIterTimeInfection = mcmc.getNIterTimeInf();

    // Initial state
    mcmc.initial_param_values(paramNames, 1);
	mcmc.initialize_augmented_data(data_augmentation, 1); // Initialize colonization times
	mcmc.initial_log_lik(0); 
    cout << "Initial log likelihood (start 1): " << mcmc.globalLogLik() << endl;
    int max_n_start = 1;
    while ( isinf(mcmc.globalLogLik()) ) {
        max_n_start++;
        mcmc.initial_param_values(paramNames, 1);
        if (data_augmentation == "yes") mcmc.initialize_augmented_data(data_augmentation, 0);
        mcmc.initial_log_lik(0);
        cout << "Initial log likelihood (start " << max_n_start <<"): " << mcmc.globalLogLik() << endl;
        if (max_n_start == 10) break;
    }

    if (! isinf(mcmc.globalLogLik()) ) {

        // Column names
        std::string colNames="iteration logLik ";
        for (const auto &p : paramNames) colNames += p + " " + p + "_p " + p + "_a ";
        colNames += "data_p data_a";
        output << colNames << endl;

        // Likelihood initial state
        output << "0 " << mcmc.globalLogLik() << " "; // Log likelihood
        for (size_t i = 0; i < mcmc.nParameters(); i++)
            output << mcmc.parameter(i) << " 0 0 ";
        output << "0 0" << endl;


        // Chain
        for (iteration = 0; iteration < numberOfIteration; iteration++)
        {
            mcmc.resetMoves();

            for (iter = 0; iter < pas; iter++)
            {
                for (size_t selectedParameter = 0; selectedParameter < idOfSelectedParameter.size(); selectedParameter++)
                {
                    parID = idOfSelectedParameter[selectedParameter];
                    mcmc.update_parameter(parID, iteration*pas+iter+1);
                }

                // Augmented data
                if (data_augmentation == "yes") {
                    for (int i=0; i < nIterTimeInfection; i++) 
                    {
                        mcmc.update_augmented_acquisitions_times(); 
                    }
                }

            }

            // Write log likelihood, parameter values, number of proposed/accepted move per parameter in the output file
            output << (iteration + 1) * pas << " " << mcmc.globalLogLik() << " "; // Log likelihood

            for (size_t i = 0; i < mcmc.nParameters(); i++)
                output << mcmc.parameter(i) << " " << mcmc.proposedMove(i) << " " << mcmc.acceptedMove(i) << " ";

            output << mcmc.proposedMoveData() << " " << mcmc.acceptedMoveData() << endl;
        }

        output.close();
    } else {
        cout << "Diverging likelihood after 10 starts" << endl;
    }

}


//----------------------------------------------------------------------
// Get list of parameters in the model
//----------------------------------------------------------------------
std::vector<int> getParamModel(std::string model, int numberOfParameters) {
    std::vector<int> params(numberOfParameters, 0);

    // intercept1, pcov2, pcov3, pcov4, intercept2, intercept3, intercept4, pcov
    //     0        1      2      3         4         5           6           7  

    if (model == "model1") {
        params[0] = 1;
        params[4] = 1;
        params[5] = 1;      
    }

    if (model=="model2") {
        params[4] = 1;
        params[1] = 1;
        params[2] = 1;
    }

    return params;
}



//----------------------------------------------------------------------
// Main function
//----------------------------------------------------------------------
int main(int argc, char **argv)
{


    cout << "================================" << endl;
    cout << "BMR TRANSMISSION MODEL" << endl;
    cout << "================================" << "\n\n";


    // Arguments passed to main
    std::string chainID = argv[1];
    std::string dataType = argv[2];
    std::string nSim = argv[3];
    std::string data_augmentation = "yes";
    std::string model;
    if (argc > 4) {
        model = argv[4];
        data_augmentation = argv[5];
    }

    //==========Model parameters==========
    // Initial values
    std::vector<std::string> parameterNames = {"intercept1", "pcov2", "pcov3", "pcov4", "intercept2", "intercept3", "intercept4", "pcov"};
    //                                            0            1        2        3         4          5              6              7      
    
    int numberOfParameters = parameterNames.size();
    std::vector<double> parameter(numberOfParameters);

    // Parameters to infer according to model
    std::vector<int> selectedParameter = getParamModel(model, numberOfParameters);

	int parameterNumber;
	std::vector<int> idOfSelectedParameter(0);
    for(parameterNumber=0; parameterNumber<numberOfParameters; parameterNumber++)
    {
        if(selectedParameter[parameterNumber]==1) idOfSelectedParameter.push_back(parameterNumber);
    }

    cout << "ID of selected parameters: "; 
    for (auto i = idOfSelectedParameter.begin(); i != idOfSelectedParameter.end(); ++i)
    	std::cout << parameterNames[*i] << ' ';
    cout << endl;


    //==========MCMC parameters==========
    int pas = 100; 
    int numberOfIteration = 100000;
    int numberOfIterationTimeInfection = 1;

    // Variance of random walk
    std::vector<double> rateForRandomWalk(numberOfParameters);
    //                                         intercept1 pcov2  pcov3   pcov4    intercept2   intercept3   intercept4      pcov
    if (model == "model1") rateForRandomWalk = { 1.8,     0.,    0.,     0.,      0.7,         1.3,         0.,             0.};
    if (model == "model2") rateForRandomWalk = { 0.,      1.,    6.0,    0.,      0.8,         0.,          0.,             0. };


    //==========Output files==========
    //Paths
    std::string dataDir="/pasteur/appa/homes/maylayan/MMMICovid/BMR_Covid19_guadeloupe/data/";
    std::string outputDir="/pasteur/appa/homes/maylayan/MMMICovid/BMR_Covid19_guadeloupe/results/";
    fs::path pathOutput = outputDir;
    if (fs::is_directory(pathOutput) == 0) fs::create_directories(pathOutput);  

    // Files
    std::string episodeDataFile, covidPrevalenceDataFile, nbedsDataFile;
    std::string outputFile, augmentedDataFile; 
    
    covidPrevalenceDataFile=dataDir+"covid_data_sim.txt";        
    nbedsDataFile=dataDir+"nbeds_data_sim.txt";  
    
    if (dataType == "simulation") {
        if (data_augmentation=="no") {
            episodeDataFile=dataDir+"simulated_fixed/"+model+"/episode_data_sim"+nSim+".txt";
            outputDir += "simulated_fixed/"+model;
        }
        
        if (data_augmentation=="yes") {
            episodeDataFile=dataDir+"simulated/"+model+"/episode_data_sim"+nSim+".txt";
            outputDir += "simulated/"+model;
        }

        fs::path pathFile = outputDir;
        if (fs::is_directory(pathFile) == 0) {
            fs::create_directories(pathFile);
        }
        outputFile= outputDir+"/mcmc_sim" + nSim + "_" + chainID + ".txt";
        augmentedDataFile=outputDir+"/augmented_data_sim" + nSim + "_" + chainID + ".txt"; 

    } else {
        episodeDataFile=dataDir+"empirical/episode_data.txt";            
        outputFile= outputDir + "empirical/mcmc_" + chainID + ".txt";
        augmentedDataFile=outputDir + "empirical/augmented_data_" + chainID + ".txt"; 
    }

    cout << "Input files: " << endl;
    cout << episodeDataFile << endl;
    cout << covidPrevalenceDataFile << endl;
    cout << nbedsDataFile << "\n\n";
    cout << "Output file: " << endl << outputFile << "\n\n";
    
    
    //==========Build data==========
    // Load data
    std::vector<Period> periodData = buildData(
        episodeDataFile, 
        covidPrevalenceDataFile, 
        nbedsDataFile, 
        data_augmentation
        );

    // Initialize MCMC object
    McmcObject mcmc(
        numberOfIteration, 
        periodData, 
        parameter, 
        selectedParameter, 
        rateForRandomWalk, 
        numberOfIterationTimeInfection
        );

    //==========MCMC==========
    runMCMC(
        mcmc,
        outputFile,
        augmentedDataFile,
        pas,
        idOfSelectedParameter,
        parameterNames,
        data_augmentation
        );

    return 0;
}
