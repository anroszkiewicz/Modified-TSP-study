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

std::pair<Solution, int> evolutionaryAlgorithm(const std::vector<std::vector<double>> &distanceMatrix, int timeLimit)
{
    auto t1 = std::chrono::high_resolution_clock::now();

    Solution x = randomCycle(distanceMatrix);
    x = localSearchMemory(x, distanceMatrix);

    long runtime = 0;
    int iterations = 0;
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    while (runtime < timeLimit)
    {
        auto t2 = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

        Solution y = x;

        y = localSearchMemory(y, distanceMatrix);

        x.calculateScore(distanceMatrix);
        y.calculateScore(distanceMatrix);

        if (y.score < x.score)
            x = y;
        iterations++;
    }

    return std::make_pair(x, iterations);
}