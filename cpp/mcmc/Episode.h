#ifndef DEF_EPISODE__H
#define DEF_EPISODE__H
#include <string>
#include <vector>
#include <map>
#include <random>


//----------------------------------------------------------------------
// Episode class
//----------------------------------------------------------------------
class Episode
{
public:
	Episode();
	~Episode() {};
	
	std::string type() const { return m_type; };
	size_t subsector() const { return m_subsector; };
	size_t id() const { return m_id; };
	size_t colTime() const { return m_colTime; };
	int firstPositiveTest() const { return m_firstPositiveTest; };
	int lastNegativeTest() const { return m_lastNegativeTest; };
	size_t admissionTime() const { return m_admissionTime; };
	size_t discharge() const { return m_discharge; };

	void newEpisode(size_t id, size_t subsector, int admission, int discharge, int detectionDay, int lastNegative, std::string type, std::string data_augmentation);
	int getTransmissionStatus(int day);
	int getStay(int day);

	int newColTime(std::mt19937_64& gen);
	void setColTime(int colTime);

	double compute_log_lik_neg(
		std::vector<double> parameters, 
		std::vector<int> selectedParam,
		int period,
		std::map<size_t, size_t> all_subsectors,
		std::vector<size_t> subsectorsInInteraction,
  		std::vector<double> colonizedPrev,
  		std::vector<std::vector<double>> covidPrev,
  		int display = 0
  		);

	double compute_log_lik_acq(
		std::vector<double> parameters, 
		std::vector<int> selectedParam,
		int period,
		std::map<size_t, size_t> all_subsectors,
		std::vector<size_t> subsectorsInInteraction,
  		std::vector<double> colonizedPrev,
  		std::vector<std::vector<double>> covidPrev,
  		int display = 0
  		);

private:
	size_t m_id;
	size_t m_subsector;
	int m_admissionTime;
	int m_discharge;
	int m_lastNegativeTest;
	int m_firstPositiveTest;
	int m_colTime;
	std::string m_type;

};

#endif