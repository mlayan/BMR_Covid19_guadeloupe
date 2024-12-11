#ifndef DEF_UTILS__H
#define DEF_UTILS__H
#include <random>

// Random generators
double rnorm(std::mt19937_64& gen, double mean, double sd);
double runif(std::mt19937_64& gen, double lower_value = 0.0, double upper_value = 1.0);
int rdiscrete(std::mt19937_64& gen, int range);

#endif
