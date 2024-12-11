#ifndef MCMCOBJECT__H
#define MCMCOBJECT__H
#include <random>
#include <vector>
#include "Period.h"


// Class
class McmcObject
{
public:
	McmcObject();
	McmcObject(
            int nIterations,
            std::vector<Period> data,
            std::vector<double> parameter,
            std::vector<int> selectedParameter,
            std::vector<double> rateForRandomWalk,
            int nIterTimeInfection = 1
            );
	~McmcObject() {};

	size_t getNumbPeriod() const { return m_nPeriod; };
	int iteration() const { return m_iterations; };
	int getNIterTimeInf() const { return m_iterTimeInfection; };
	int proposedMove(int index) const { return m_numberOfMoveProposed[index]; };
	int acceptedMove(int index) const { return m_numberOfMoveAccepted[index]; };
	int proposedMoveData() const { return m_proposedMoveData; };
	int acceptedMoveData() const { return m_acceptedMoveData; };
	double globalLogLik() const { return m_globalLogLik; };
	double periodLogLik(int index) const { return m_periodLogLik[index]; };
	std::vector<double> hhLogLik() const { return m_periodLogLik; };
	double rateRandomWalk(int index) const { return m_rateForRandomWalk[index]; };
	double parameter(int index) const { return m_parameter[index]; };
	size_t nParameters() const { return m_parameter.size(); };

	void resetMoves();
	void initial_param_values(std::vector<std::string> paramNames, int display = 0);
	void initialize_augmented_data(std::string data_augmentation, int display = 0);
	void initial_log_lik(int display = 0);
	void update_parameter(int parID, double step);
	void update_augmented_acquisitions_times();

private:
	int m_iterations;
	int m_iterTimeInfection;
	std::mt19937_64 m_gen;
	double m_globalLogLik;
	std::vector<double> m_periodLogLik;
	int m_acceptedMoveData;
	int m_proposedMoveData;
	std::vector<int> m_numberOfMoveAccepted;
	std::vector<int> m_numberOfMoveProposed;
	std::vector<double> m_accept;
	std::vector<double> m_parameter;
	std::vector<int>  m_selectedParam;
	std::vector<double> m_rateForRandomWalk;
	std::vector<Period> m_data;
	size_t m_nPeriod;
};

#endif
