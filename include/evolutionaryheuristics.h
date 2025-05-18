#ifndef EVOLUTIONARYHEURISTICS_H
#define EVOLUTIONARYHEURISTICS_H

#include <vector>
#include <utility>

bool insertIntoPopulation(std::vector<Solution> &population, const Solution &child);

Solution crossover(const Solution &parent1, const Solution &parent2, const std::vector<std::vector<double>> &distanceMatrix);

Solution proposed_crossover(const Solution &parent1, const Solution &parent2, const std::vector<std::vector<double>> &distanceMatrix);

std::pair<Solution, int> evolutionaryAlgorithm(const std::vector<std::vector<double>> &distanceMatrix, int timeLimit, bool localSearch);

#endif