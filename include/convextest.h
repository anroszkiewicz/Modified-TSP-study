#ifndef CONVEXTEST_H
#define CONVEXTEST_H

#include <vector>
#include <utility>

int common_vertices_metric(Solution &a, Solution &b);

int common_edges_metric(Solution &a, Solution &b);

void convex_test(Solution &goodSolution, std::vector<Solution> &solutions, const std::vector<std::vector<double>> &distanceMatrix, std::string metric);

#endif