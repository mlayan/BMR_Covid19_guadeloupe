#include <cmath>
#include <random>
#include <algorithm>
#include "utils.h"

using namespace std;


//----------------------------------------------------------------------
// Random generator
//----------------------------------------------------------------------
double rnorm(std::mt19937_64& gen, double mean, double sd) {
    std::normal_distribution<double> dist(mean, sd);
    return dist(gen);
}

double runif(std::mt19937_64& gen, double lower_value, double upper_value) {
    std::uniform_real_distribution<double> dist(lower_value, upper_value);
    return dist(gen);
}

int rdiscrete(std::mt19937_64& gen, int range) {
    std::vector<int> weights(range, 1);
    std::discrete_distribution<> dist(weights.begin(), weights.end());
    return dist(gen);
}
