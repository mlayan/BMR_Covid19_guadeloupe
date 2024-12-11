#ifndef DEF_PERIOD__H
#define DEF_PERIOD__H
#include <string>
#include <vector>
#include <random>
#include "Episode.h"

//----------------------------------------------------------------------
// Period class
//----------------------------------------------------------------------
class Period
{
public:
	Period();
	~Period() {};
	
	std::map<size_t, size_t> getSubsectorMap() const { return m_subsectors; };
	size_t getSize() const { return m_size; };
	size_t getDays() const { return m_nDays; };
	size_t getID() const { return m_periodID; };
	double getColTime(int index) const { return m_episodes[index].colTime(); };
	size_t getSubsector(int index) const { return m_episodes[index].subsector(); };
	size_t getDischarge(int index) const { return m_episodes[index].discharge(); };
	std::vector<size_t> getAcquisitions() const { return m_acquisitions; };
	std::vector<size_t> getPositives() const { return m_positives; };
	std::vector<size_t> getNegatives() const { return m_negatives; };

	void newPeriod();
	void addEpisode(Episode episode, int period);
	void addCovidPrev(std::vector<std::vector<double>> prevalenceMatrix, int display=0);
	void addNBeds(std::vector<std::vector<double>> demographicsMatrix, int display=0);
	void initialColonizationTimes(std::mt19937_64& gen, int display=0);
	void compute_prevalences(int display = 0);
	void update_prevalences(int epi, int display = 0);
	void display_data();
	void display_prevalence_matrix();

	int newColTime(int index, std::mt19937_64& gen);
	void setColTime(int index, int colTime);

	double compute_log_lik(std::vector<double> parameters, 
		std::vector<int> selectedParam, 
		int display = 0);

private:
	size_t m_size;
	size_t m_nDays;
	int m_periodID;
	std::vector<Episode> m_episodes;
	std::map<size_t, size_t> m_subsectors;
	std::vector<std::vector<size_t>> m_idsPerSubsector;
	std::vector<std::vector<double>> m_colonizedPerSubsector;
	std::vector<std::vector<double>> m_patientPerSubsector;
	std::vector<double> m_colonizedPrev;
	//std::vector<double> m_patientsTot;
	std::vector<std::vector<double>> m_covidPrev;
	std::vector<std::vector<double>> m_nbeds;
	std::vector<size_t> m_positives;
	std::vector<size_t> m_negatives;
	std::vector<size_t> m_acquisitions;
	std::vector<size_t> m_subsectorsInInteraction;
};

#endif
