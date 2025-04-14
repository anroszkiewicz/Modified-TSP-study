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
#include <utility>

#include "point.h"
#include "solution.h"
#include "move.h"
#include "localsearch.h"
#include "optimization.h"
#include "extendedlocalsearch.h"

Solution multipleStartLocalSearch(const std::vector<std::vector<double>> &distanceMatrix)
{
    Solution bestSolution;
    double bestScore = std::numeric_limits<double>::max();

    for(int i=0; i<200; i++)
    {
        std::srand(static_cast<unsigned int>(std::time(nullptr)));
        Solution initialRandom = randomCycle(distanceMatrix);
        Solution newSolution = localSearchMemory(initialRandom, distanceMatrix);

        newSolution.calculateScore(distanceMatrix);
        double newScore = newSolution.score;
        if (newScore < bestScore)
        {
            bestSolution = newSolution;
            bestScore = newScore;
        }
    }
    return bestSolution;
}

std::pair <Solution, int> iteratedLocalSearch(const std::vector<std::vector<double>> &distanceMatrix, int timeLimit)
{
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    Solution initialRandom = randomCycle(distanceMatrix);
    
    auto t1 = std::chrono::high_resolution_clock::now();
    long runtime = 0;
    int iterations = 0;

    // while(runtime < timeLimit)
    // {
    //     TODO
    //     auto t2 = std::chrono::high_resolution_clock::now();
    //     runtime = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
    //     iterations++;
    // }

    return std::make_pair(initialRandom, iterations);
}

std::pair <Solution, int> largeNeighborhoodSearch(const std::vector<std::vector<double>> &distanceMatrix, int timeLimit)
{
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    Solution initialRandom = randomCycle(distanceMatrix);

    auto t1 = std::chrono::high_resolution_clock::now();
    long runtime = 0;
    int iterations = 0;

    // while(runtime < timeLimit)
    // {
    //     TODO
    //     auto t2 = std::chrono::high_resolution_clock::now();
    //     runtime = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
    //     iterations++;
    // }

    return std::make_pair(initialRandom, iterations);
}