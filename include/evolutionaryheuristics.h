#ifndef EVOLUTIONARYHEURISTICS_H
#define EVOLUTIONARYHEURISTICS_H

#include <vector>
#include <utility>

std::pair<Solution, int> evolutionaryAlgorithm(const std::vector<std::vector<double>> &distanceMatrix, int timeLimit);

#endif