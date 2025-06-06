#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <limits>
#include <random>
#include <chrono>
#include <queue>
#include <algorithm>

#include "point.h"
#include "solution.h"
#include "move.h"
#include "localsearch.h"
#include "optimization.h"
#include "extendedlocalsearch.h"
#include "utils.h"
#include "greedyheuristics.h"
#include "evolutionaryheuristics.h"

bool point_in_cycle(int point, std::vector <int> &cycle)
{
    return std::find(cycle.begin(), cycle.end(), point) != cycle.end();
}

int common_vertices_metric(Solution &a, Solution &b)
{
    int count = 0;
    for (size_t i=0; i<200; ++i)
    {
        for (size_t j=0; j<200; ++j)
        {
            if (i == j) continue;
            if (point_in_cycle(i, a.cycleIndices[0]) && point_in_cycle(j, a.cycleIndices[0]))
            {
                if (point_in_cycle(i, b.cycleIndices[0]) && point_in_cycle(j, b.cycleIndices[0])) count++;
                if (point_in_cycle(i, b.cycleIndices[1]) && point_in_cycle(j, b.cycleIndices[1])) count++;
            }
            if (point_in_cycle(i, a.cycleIndices[1]) && point_in_cycle(j, a.cycleIndices[1]))
            {
                if (point_in_cycle(i, b.cycleIndices[0]) && point_in_cycle(j, b.cycleIndices[0])) count++;
                if (point_in_cycle(i, b.cycleIndices[1]) && point_in_cycle(j, b.cycleIndices[1])) count++;
            }
        }
    }
    return count;
}

int common_edges_metric(Solution &a, Solution &b)
{
    int count = 0;
    for (int cycle= 0; cycle<2; ++cycle)
    {
        for (size_t i = 1; i < a.cycleIndices[cycle].size(); i++)
        {
            // get edge
            std::vector<int> sequence = {a.cycleIndices[cycle][i-1], a.cycleIndices[cycle][i]};

            // check first cycle
            auto it = std::search(b.cycleIndices[0].begin(), b.cycleIndices[0].end(), sequence.begin(), sequence.end());
            if (it != b.cycleIndices[0].end()) count++;

            // check second cycle
            it = std::search(b.cycleIndices[1].begin(), b.cycleIndices[1].end(), sequence.begin(), sequence.end());
            if (it != b.cycleIndices[1].end()) count++;

            // check reversed edge
            sequence = {a.cycleIndices[cycle][i], a.cycleIndices[cycle][i-1]};

            // check first cycle
            it = std::search(b.cycleIndices[0].begin(), b.cycleIndices[0].end(), sequence.begin(), sequence.end());
            if (it != b.cycleIndices[0].end()) count++;

            // check second cycle
            it = std::search(b.cycleIndices[1].begin(), b.cycleIndices[1].end(), sequence.begin(), sequence.end());
            if (it != b.cycleIndices[1].end()) count++;

            // check first and last element
            if (i-1 == b.cycleIndices[0][0] && i == b.cycleIndices[0][99]) count++;
            if (i == b.cycleIndices[0][0] && i-1 == b.cycleIndices[0][99]) count++;
        }
    }
    return count;
}

void convex_test(Solution &goodSolution, std::vector<Solution> &solutions, const std::vector<std::vector<double>> &distanceMatrix, std::string metric="vertices")
{
    // sort solutions
    sort(solutions.begin(), solutions.end());

    // get similarity between solutions and good solutions
    std::vector <double> similarities;
    double result;
    for (int i = 0; i < static_cast<int>(solutions.size()); ++i)
    {
        if (metric == "vertices") result = (double)common_vertices_metric(goodSolution, solutions[i]);
        else result = (double)common_edges_metric(goodSolution, solutions[i]);
            
        similarities.push_back(result);
    }

    // draw plot
    std::string title = "Podobieństwo do bardzo dobrego rozwiązania";
    if (metric == "vertices") title += ", miara wierzchołków";
    else title += ", miara krawędzi";
    plotSimilarity(solutions, similarities, distanceMatrix, title);

    // get average similarity betweens solutions
    std::vector <double> average_similarities;
    int similarity_sum;
    double average;
    std::cout << "Calculating similarity" << std::endl;
    for (int i = 0; i < static_cast<int>(solutions.size()); ++i)
    {
        std::cout << "Progress: " << i << "/1000" << std::endl;
        similarity_sum = 0;
        for (int j = 0; j < static_cast<int>(solutions.size()); ++j)
        {
            if (i == j) continue;

            if (metric == "vertices") similarity_sum += common_vertices_metric(solutions[i], solutions[j]);
            else similarity_sum = common_edges_metric(solutions[i], solutions[j]);
            
        }
        average = (double)similarity_sum / 1000;
        average_similarities.push_back(average);
    }

    // draw plot
    title = "Średnie podobieństwo";
    if (metric == "vertices") title += ", miara wierzchołków";
    else title += ", miara krawędzi";
    plotSimilarity(solutions, average_similarities, distanceMatrix, title);
}