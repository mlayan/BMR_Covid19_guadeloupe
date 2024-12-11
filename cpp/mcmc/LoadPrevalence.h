#ifndef DEF_LOAD_PREV__H
#define DEF_LOAD_PREV__H
#include<vector>

std::vector<std::vector<std::vector<double>>> getPrevalenceMatrix(
	std::string fileName, 
	int numberOfPeriodsColonization, 
	int display=0
	);

#endif
