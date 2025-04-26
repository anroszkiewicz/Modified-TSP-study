#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <limits>
#include <random>
#include <chrono>
#include <queue>
#include <utility>

#include "utils.h"
#include "solution.h"
#include "point.h"
#include "move.h"
#include "cut.h"
#include "vertexneighbourhood.h"
#include "greedyheuristics.h"
#include "localsearch.h"
#include "optimization.h"
#include "extendedlocalsearch.h"

void lab1(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix)
{
    // Initialize variables for tracking scores
    int iterations = 100;

    // GreedyNearestNeighbour Algorithm
    double minGreedyNNScore = std::numeric_limits<double>::max();
    double maxGreedyNNScore = std::numeric_limits<double>::lowest();
    double totalGreedyNNScore = 0.0;
    Solution bestGreedyNNSolution;

    // GreedyCycle Algorithm
    double minGreedyCycleScore = std::numeric_limits<double>::max();
    double maxGreedyCycleScore = std::numeric_limits<double>::lowest();
    double totalGreedyCycleScore = 0.0;
    Solution bestGreedyCycleSolution;

    // RegretCycle Algorithm
    double minRegretCycleScore = std::numeric_limits<double>::max();
    double maxRegretCycleScore = std::numeric_limits<double>::lowest();
    double totalRegretCycleScore = 0.0;
    Solution bestRegretCycleSolution;

    // RegretCycleWeighted Algorithm
    double minWeightedRegretCycleScore = std::numeric_limits<double>::max();
    double maxWeightedRegretCycleScore = std::numeric_limits<double>::lowest();
    double totalWeightedRegretCycleScore = 0.0;
    Solution bestWeightedRegretCycleSolution;

    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    for (int i = 0; i < iterations; ++i)
    {
        // std::cout << "Progress: " << i << "/" << iterations << std::endl;
        // Apply greedyNearestNeighbour algorithm
        Solution solution1 = greedyNearestNeighbour(distanceMatrix);
        solution1.calculateScore(distanceMatrix);
        double greedyNNScore = solution1.getScore();
        minGreedyNNScore = std::min(minGreedyNNScore, greedyNNScore);
        maxGreedyNNScore = std::max(maxGreedyNNScore, greedyNNScore);
        totalGreedyNNScore += greedyNNScore;

        // Check for best solution for GreedyNN
        if (greedyNNScore == minGreedyNNScore)
        {
            minGreedyNNScore = greedyNNScore;
            bestGreedyNNSolution = solution1;
        }

        // Apply greedyCycle algorithm
        Solution solution2 = greedyCycle(distanceMatrix);
        solution2.calculateScore(distanceMatrix);
        double greedyCycleScore = solution2.getScore();
        minGreedyCycleScore = std::min(minGreedyCycleScore, greedyCycleScore);
        maxGreedyCycleScore = std::max(maxGreedyCycleScore, greedyCycleScore);
        totalGreedyCycleScore += greedyCycleScore;

        // Check for best solution for GreedyCycle
        if (greedyCycleScore == minGreedyCycleScore)
        {
            minGreedyCycleScore = greedyCycleScore;
            bestGreedyCycleSolution = solution2;
        }

        // Apply regretCycle algorithm
        Solution empty;
        Solution solution3 = regretCycleWeighted(empty, distanceMatrix, 0.0, 1.0);
        solution3.calculateScore(distanceMatrix);
        double regretCycleScore = solution3.getScore();
        minRegretCycleScore = std::min(minRegretCycleScore, regretCycleScore);
        maxRegretCycleScore = std::max(maxRegretCycleScore, regretCycleScore);
        totalRegretCycleScore += regretCycleScore;

        // Check for best solution for RegretCycle
        if (regretCycleScore == minRegretCycleScore)
        {
            minRegretCycleScore = regretCycleScore;
            bestRegretCycleSolution = solution3;
        }

        // Apply regretCycleWeighted algorithm
        Solution solution4 = regretCycleWeighted(empty, distanceMatrix, 1.0, 1.0);
        solution4.calculateScore(distanceMatrix);
        double weightedRegretCycleScore = solution4.getScore();
        minWeightedRegretCycleScore = std::min(minWeightedRegretCycleScore, weightedRegretCycleScore);
        maxWeightedRegretCycleScore = std::max(maxWeightedRegretCycleScore, weightedRegretCycleScore);
        totalWeightedRegretCycleScore += weightedRegretCycleScore;

        // Check for best solution for RegretCycleWeighted
        if (weightedRegretCycleScore == minWeightedRegretCycleScore)
        {
            minWeightedRegretCycleScore = weightedRegretCycleScore;
            bestWeightedRegretCycleSolution = solution4;
        }
    }

    // Calculate average scores for each algorithm
    double avgGreedyNNScore = totalGreedyNNScore / iterations;
    double avgGreedyCycleScore = totalGreedyCycleScore / iterations;
    double avgRegretCycleScore = totalRegretCycleScore / iterations;
    double avgWeightedRegretCycleScore = totalWeightedRegretCycleScore / iterations;

    // Print results for each algorithm
    std::cout << "\nGreedyNearestNeighbour Algorithm Results:\n";
    std::cout << "Min: " << minGreedyNNScore << " Max: " << maxGreedyNNScore
              << " Avg: " << avgGreedyNNScore << "\n";
    plotSolution(bestGreedyNNSolution, points, distanceMatrix, "GreedyNearestNeighbour");

    std::cout << "\nGreedyCycle Algorithm Results:\n";
    std::cout << "Min: " << minGreedyCycleScore << " Max: " << maxGreedyCycleScore
              << " Avg: " << avgGreedyCycleScore << "\n";
    plotSolution(bestGreedyCycleSolution, points, distanceMatrix, "GreedyCycle");

    std::cout << "\nRegretCycle Algorithm Results:\n";
    std::cout << "Min: " << minRegretCycleScore << " Max: " << maxRegretCycleScore
              << " Avg: " << avgRegretCycleScore << "\n";
    plotSolution(bestRegretCycleSolution, points, distanceMatrix, "RegretCycle");

    std::cout << "\nRegretCycleWeighted Algorithm Results:\n";
    std::cout << "Min: " << minWeightedRegretCycleScore << " Max: " << maxWeightedRegretCycleScore
              << " Avg: " << avgWeightedRegretCycleScore << "\n";
    plotSolution(bestWeightedRegretCycleSolution, points, distanceMatrix, "RegretCycleWeighted");
}

void lab2(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix)
{
    int iterations = 100;

    /*
        In this experiment, we evaluate different approaches to solving the 2-cycle TSP.
        We have three independent choices, leading to 2^3 = 8 possible solution strategies:

        1. **Starting Solution**
           - `Random` (Completely randomized initial solution)
           - `Greedy` (Constructed using a heuristic like Nearest Neighbor or Greedy Cycle)

        2. **Local Search Strategy**
           - `Steepest Descent` (Evaluates all possible improvements and picks the best one)
           - `Greedy` (Applies the first improving move found)

        3. **Modification Type**
           - `Exchange Points` (Swaps individual points between cycles)
           - `Exchange Edges` (Swaps entire edges to optimize the path)
    */

    // 0. Random Walk (before all random-based approaches)
    double minRandomWalk1 = std::numeric_limits<double>::max();
    double maxRandomWalk1 = std::numeric_limits<double>::lowest();
    double totalRandomWalk1 = 0.0;
    Solution bestRandomWalk1;
    long minRandomWalk1time = std::numeric_limits<long>::max();
    long maxRandomWalk1time = std::numeric_limits<long>::lowest();
    long totalRandomWalk1time = 0.0;

    // 1. Random + Steepest Descent + Point Exchange
    double minRandomSteepestPointScore = std::numeric_limits<double>::max();
    double maxRandomSteepestPointScore = std::numeric_limits<double>::lowest();
    double totalRandomSteepestPointScore = 0.0;
    Solution bestRandomSteepestPointSolution;
    long minRandomSteepestPointTime = std::numeric_limits<long>::max();
    long maxRandomSteepestPointTime = std::numeric_limits<long>::lowest();
    long totalRandomSteepestPointTime = 0.0;

    // 2. Random + Steepest Descent + Edge Exchange
    double minRandomSteepestEdgeScore = std::numeric_limits<double>::max();
    double maxRandomSteepestEdgeScore = std::numeric_limits<double>::lowest();
    double totalRandomSteepestEdgeScore = 0.0;
    Solution bestRandomSteepestEdgeSolution;
    long minRandomSteepestEdgeTime = std::numeric_limits<long>::max();
    long maxRandomSteepestEdgeTime = std::numeric_limits<long>::lowest();
    long totalRandomSteepestEdgeTime = 0.0;

    // 3. Random + Greedy Local Search + Point Exchange
    double minRandomGreedyPointScore = std::numeric_limits<double>::max();
    double maxRandomGreedyPointScore = std::numeric_limits<double>::lowest();
    double totalRandomGreedyPointScore = 0.0;
    Solution bestRandomGreedyPointSolution;
    long minRandomGreedyPointTime = std::numeric_limits<long>::max();
    long maxRandomGreedyPointTime = std::numeric_limits<long>::lowest();
    long totalRandomGreedyPointTime = 0.0;

    // 4. Random + Greedy Local Search + Edge Exchange
    double minRandomGreedyEdgeScore = std::numeric_limits<double>::max();
    double maxRandomGreedyEdgeScore = std::numeric_limits<double>::lowest();
    double totalRandomGreedyEdgeScore = 0.0;
    Solution bestRandomGreedyEdgeSolution;
    long minRandomGreedyEdgeTime = std::numeric_limits<long>::max();
    long maxRandomGreedyEdgeTime = std::numeric_limits<long>::lowest();
    long totalRandomGreedyEdgeTime = 0.0;

    // 5. Random Walk (before all greedy-based approaches)
    double minRandomWalk2 = std::numeric_limits<double>::max();
    double maxRandomWalk2 = std::numeric_limits<double>::lowest();
    double totalRandomWalk2 = 0.0;
    Solution bestRandomWalk2;

    // 6. Greedy + Steepest Descent + Point Exchange
    double minGreedySteepestPointScore = std::numeric_limits<double>::max();
    double maxGreedySteepestPointScore = std::numeric_limits<double>::lowest();
    double totalGreedySteepestPointScore = 0.0;
    Solution bestGreedySteepestPointSolution;
    long minGreedySteepestPointTime = std::numeric_limits<long>::max();
    long maxGreedySteepestPointTime = std::numeric_limits<long>::lowest();
    long totalGreedySteepestPointTime = 0.0;

    // 7. Greedy + Steepest Descent + Edge Exchange
    double minGreedySteepestEdgeScore = std::numeric_limits<double>::max();
    double maxGreedySteepestEdgeScore = std::numeric_limits<double>::lowest();
    double totalGreedySteepestEdgeScore = 0.0;
    Solution bestGreedySteepestEdgeSolution;
    long minGreedySteepestEdgeTime = std::numeric_limits<long>::max();
    long maxGreedySteepestEdgeTime = std::numeric_limits<long>::lowest();
    long totalGreedySteepestEdgeTime = 0.0;

    // 8. Greedy + Greedy Local Search + Point Exchange
    double minGreedyGreedyPointScore = std::numeric_limits<double>::max();
    double maxGreedyGreedyPointScore = std::numeric_limits<double>::lowest();
    double totalGreedyGreedyPointScore = 0.0;
    Solution bestGreedyGreedyPointSolution;
    long minGreedyGreedyPointTime = std::numeric_limits<long>::max();
    long maxGreedyGreedyPointTime = std::numeric_limits<long>::lowest();
    long totalGreedyGreedyPointTime = 0.0;

    // 9. Greedy + Greedy Local Search + Edge Exchange
    double minGreedyGreedyEdgeScore = std::numeric_limits<double>::max();
    double maxGreedyGreedyEdgeScore = std::numeric_limits<double>::lowest();
    double totalGreedyGreedyEdgeScore = 0.0;
    Solution bestGreedyGreedyEdgeSolution;
    long minGreedyGreedyEdgeTime = std::numeric_limits<long>::max();
    long maxGreedyGreedyEdgeTime = std::numeric_limits<long>::lowest();
    long totalGreedyGreedyEdgeTime = 0.0;

    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    Solution empty;

    for (int i = 0; i < iterations; ++i)
    {
        std::cout << "Progress: " << i << "/" << iterations << std::endl;

        // Generate initial random and greedy solutions
        Solution initialRandom = randomCycle(distanceMatrix);
        Solution initialGreedy = regretCycleWeighted(empty, distanceMatrix, 1.0, 1.0);

        // 0. Random Walk 1
        auto t1 = std::chrono::high_resolution_clock::now();
        Solution solution0 = randomWalk(initialRandom, distanceMatrix, 236);
        auto t2 = std::chrono::high_resolution_clock::now();

        solution0.calculateScore(distanceMatrix);
        double score0 = solution0.getScore();
        minRandomWalk1 = std::min(minRandomWalk1, score0);
        maxRandomWalk1 = std::max(maxRandomWalk1, score0);
        totalRandomWalk1 += score0;
        if (score0 == minRandomWalk1)
            bestRandomWalk1 = solution0;

        auto time0 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minRandomWalk1time = std::min(minRandomWalk1time, static_cast<long>(time0));
        maxRandomWalk1time = std::max(maxRandomWalk1time, static_cast<long>(time0));
        totalRandomWalk1time += time0;

        // 1. Random + Steepest Descent + Point Exchange
        t1 = std::chrono::high_resolution_clock::now();
        Solution solution1 = localSearchVertex(initialRandom, distanceMatrix, "steepest");
        t2 = std::chrono::high_resolution_clock::now();

        solution1.calculateScore(distanceMatrix);
        double score1 = solution1.getScore();
        minRandomSteepestPointScore = std::min(minRandomSteepestPointScore, score1);
        maxRandomSteepestPointScore = std::max(maxRandomSteepestPointScore, score1);
        totalRandomSteepestPointScore += score1;
        if (score1 == minRandomSteepestPointScore)
            bestRandomSteepestPointSolution = solution1;

        auto time1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minRandomSteepestPointTime = std::min(minRandomSteepestPointTime, static_cast<long>(time1));
        maxRandomSteepestPointTime = std::max(maxRandomSteepestPointTime, static_cast<long>(time1));
        totalRandomSteepestPointTime += time1;

        // 2. Random + Steepest Descent + Edge Exchange
        t1 = std::chrono::high_resolution_clock::now();
        Solution solution2 = localSearchEdges(initialRandom, distanceMatrix, "steepest");
        t2 = std::chrono::high_resolution_clock::now();

        solution2.calculateScore(distanceMatrix);
        double score2 = solution2.getScore();
        minRandomSteepestEdgeScore = std::min(minRandomSteepestEdgeScore, score2);
        maxRandomSteepestEdgeScore = std::max(maxRandomSteepestEdgeScore, score2);
        totalRandomSteepestEdgeScore += score2;
        if (score2 == minRandomSteepestEdgeScore)
            bestRandomSteepestEdgeSolution = solution2;

        auto time2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minRandomSteepestEdgeTime = std::min(minRandomSteepestEdgeTime, static_cast<long>(time2));
        maxRandomSteepestEdgeTime = std::max(maxRandomSteepestEdgeTime, static_cast<long>(time2));
        totalRandomSteepestEdgeTime += time2;

        // 3. Random + Greedy Local Search + Point Exchange
        t1 = std::chrono::high_resolution_clock::now();
        Solution solution3 = localSearchVertex(initialRandom, distanceMatrix, "greedy");
        t2 = std::chrono::high_resolution_clock::now();

        solution3.calculateScore(distanceMatrix);
        double score3 = solution3.getScore();
        minRandomGreedyPointScore = std::min(minRandomGreedyPointScore, score3);
        maxRandomGreedyPointScore = std::max(maxRandomGreedyPointScore, score3);
        totalRandomGreedyPointScore += score3;
        if (score3 == minRandomGreedyPointScore)
            bestRandomGreedyPointSolution = solution3;

        auto time3 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minRandomGreedyPointTime = std::min(minRandomGreedyPointTime, static_cast<long>(time3));
        maxRandomGreedyPointTime = std::max(maxRandomGreedyPointTime, static_cast<long>(time3));
        totalRandomGreedyPointTime += time3;

        // 4. Random + Greedy Local Search + Edge Exchange
        t1 = std::chrono::high_resolution_clock::now();
        Solution solution4 = localSearchEdges(initialRandom, distanceMatrix, "greedy");
        t2 = std::chrono::high_resolution_clock::now();

        solution4.calculateScore(distanceMatrix);
        double score4 = solution4.getScore();
        minRandomGreedyEdgeScore = std::min(minRandomGreedyEdgeScore, score4);
        maxRandomGreedyEdgeScore = std::max(maxRandomGreedyEdgeScore, score4);
        totalRandomGreedyEdgeScore += score4;
        if (score4 == minRandomGreedyEdgeScore)
            bestRandomGreedyEdgeSolution = solution4;

        auto time4 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minRandomGreedyEdgeTime = std::min(minRandomGreedyEdgeTime, static_cast<long>(time4));
        maxRandomGreedyEdgeTime = std::max(maxRandomGreedyEdgeTime, static_cast<long>(time4));
        totalRandomGreedyEdgeTime += time4;

        // 5. Random Walk 2
        t1 = std::chrono::high_resolution_clock::now();
        Solution solution5 = randomWalk(initialGreedy, distanceMatrix, 236);
        t2 = std::chrono::high_resolution_clock::now();

        solution5.calculateScore(distanceMatrix);
        double score5 = solution5.getScore();
        minRandomWalk2 = std::min(minRandomWalk2, score5);
        maxRandomWalk2 = std::max(maxRandomWalk2, score5);
        totalRandomWalk2 += score5;
        if (score5 == minRandomWalk2)
            bestRandomWalk2 = solution5;

        auto time5 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

        // 6. Greedy + Steepest Descent + Point Exchange
        t1 = std::chrono::high_resolution_clock::now();
        Solution solution6 = localSearchVertex(initialGreedy, distanceMatrix, "steepest");
        t2 = std::chrono::high_resolution_clock::now();

        solution6.calculateScore(distanceMatrix);
        double score6 = solution6.getScore();
        minGreedySteepestPointScore = std::min(minGreedySteepestPointScore, score6);
        maxGreedySteepestPointScore = std::max(maxGreedySteepestPointScore, score6);
        totalGreedySteepestPointScore += score6;
        if (score6 == minGreedySteepestPointScore)
            bestGreedySteepestPointSolution = solution6;

        auto time6 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minGreedySteepestPointTime = std::min(minGreedySteepestPointTime, static_cast<long>(time6));
        maxGreedySteepestPointTime = std::max(maxGreedySteepestPointTime, static_cast<long>(time6));
        totalGreedySteepestPointTime += time6;

        // 7. Greedy + Steepest Descent + Edge Exchange
        t1 = std::chrono::high_resolution_clock::now();
        Solution solution7 = localSearchEdges(initialGreedy, distanceMatrix, "steepest");
        t2 = std::chrono::high_resolution_clock::now();

        solution7.calculateScore(distanceMatrix);
        double score7 = solution7.getScore();
        minGreedySteepestEdgeScore = std::min(minGreedySteepestEdgeScore, score7);
        maxGreedySteepestEdgeScore = std::max(maxGreedySteepestEdgeScore, score7);
        totalGreedySteepestEdgeScore += score7;
        if (score7 == minGreedySteepestEdgeScore)
            bestGreedySteepestEdgeSolution = solution7;

        auto time7 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minGreedySteepestEdgeTime = std::min(minGreedySteepestEdgeTime, static_cast<long>(time7));
        maxGreedySteepestEdgeTime = std::max(maxGreedySteepestEdgeTime, static_cast<long>(time7));
        totalGreedySteepestEdgeTime += time7;

        // 8. Greedy + Greedy Local Search + Point Exchange
        t1 = std::chrono::high_resolution_clock::now();
        Solution solution8 = localSearchVertex(initialGreedy, distanceMatrix, "greedy");
        t2 = std::chrono::high_resolution_clock::now();

        solution8.calculateScore(distanceMatrix);
        double score8 = solution8.getScore();
        minGreedyGreedyPointScore = std::min(minGreedyGreedyPointScore, score8);
        maxGreedyGreedyPointScore = std::max(maxGreedyGreedyPointScore, score8);
        totalGreedyGreedyPointScore += score8;
        if (score8 == minGreedyGreedyPointScore)
            bestGreedyGreedyPointSolution = solution8;

        auto time8 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minGreedyGreedyPointTime = std::min(minGreedyGreedyPointTime, static_cast<long>(time8));
        maxGreedyGreedyPointTime = std::max(maxGreedyGreedyPointTime, static_cast<long>(time8));
        totalGreedyGreedyPointTime += time8;

        // 9. Greedy + Greedy Local Search + Edge Exchange
        t1 = std::chrono::high_resolution_clock::now();
        Solution solution9 = localSearchEdges(initialGreedy, distanceMatrix, "greedy");
        t2 = std::chrono::high_resolution_clock::now();

        solution9.calculateScore(distanceMatrix);
        double score9 = solution9.getScore();
        minGreedyGreedyEdgeScore = std::min(minGreedyGreedyEdgeScore, score9);
        maxGreedyGreedyEdgeScore = std::max(maxGreedyGreedyEdgeScore, score9);
        totalGreedyGreedyEdgeScore += score9;
        if (score9 == minGreedyGreedyEdgeScore)
            bestGreedyGreedyEdgeSolution = solution9;

        auto time9 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minGreedyGreedyEdgeTime = std::min(minGreedyGreedyEdgeTime, static_cast<long>(time9));
        maxGreedyGreedyEdgeTime = std::max(maxGreedyGreedyEdgeTime, static_cast<long>(time9));
        totalGreedyGreedyEdgeTime += time9;
    }

    // Print results
    std::cout << "\nRandom Walk 1:\n";
    std::cout << "Min: " << minRandomWalk1 << " Max: " << maxRandomWalk1
              << " Avg: " << (totalRandomWalk1 / iterations) << "\n";
    std::cout << "Min: " << minRandomWalk1time << " Max: " << maxRandomWalk1time
              << " Avg: " << (totalRandomWalk1time / iterations) << "\n";
    plotSolution(bestRandomWalk1, points, distanceMatrix, "Random Walk 1");

    std::cout << "\nRandom + Steepest Descent + Point Exchange:\n";
    std::cout << "Min: " << minRandomSteepestPointScore << " Max: " << maxRandomSteepestPointScore
              << " Avg: " << (totalRandomSteepestPointScore / iterations) << "\n";
    std::cout << "Min: " << minRandomSteepestPointTime << " Max: " << maxRandomSteepestPointTime
              << " Avg: " << (totalRandomSteepestPointTime / iterations) << "\n";
    plotSolution(bestRandomSteepestPointSolution, points, distanceMatrix, "Random + Steepest Descent + Vertex Exchange");

    std::cout << "\nRandom + Steepest Descent + Edge Exchange:\n";
    std::cout << "Min: " << minRandomSteepestEdgeScore << " Max: " << maxRandomSteepestEdgeScore
              << " Avg: " << (totalRandomSteepestEdgeScore / iterations) << "\n";
    std::cout << "Min: " << minRandomSteepestEdgeTime << " Max: " << maxRandomSteepestEdgeTime
              << " Avg: " << (totalRandomSteepestEdgeTime / iterations) << "\n";
    plotSolution(bestRandomSteepestEdgeSolution, points, distanceMatrix, "Random + Steepest Descent + Edge Exchange");

    std::cout << "\nRandom + Greedy Local Search + Point Exchange:\n";
    std::cout << "Min: " << minRandomGreedyPointScore << " Max: " << maxRandomGreedyPointScore
              << " Avg: " << (totalRandomGreedyPointScore / iterations) << "\n";
    std::cout << "Min: " << minRandomGreedyPointTime << " Max: " << maxRandomGreedyPointTime
              << " Avg: " << (totalRandomGreedyPointTime / iterations) << "\n";
    plotSolution(bestRandomGreedyPointSolution, points, distanceMatrix, "Random + Greedy Local Search + Vertex Exchange");

    std::cout << "\nRandom + Greedy Local Search + Edge Exchange:\n";
    std::cout << "Min: " << minRandomGreedyEdgeScore << " Max: " << maxRandomGreedyEdgeScore
              << " Avg: " << (totalRandomGreedyEdgeScore / iterations) << "\n";
    std::cout << "Min: " << minRandomGreedyEdgeTime << " Max: " << maxRandomGreedyEdgeTime
              << " Avg: " << (totalRandomGreedyEdgeTime / iterations) << "\n";
    plotSolution(bestRandomGreedyEdgeSolution, points, distanceMatrix, "Random + Greedy Local Search + Edge Exchange");

    std::cout << "\nRandom Walk 2:\n";
    std::cout << "Min: " << minRandomWalk2 << " Max: " << maxRandomWalk2
              << " Avg: " << (totalRandomWalk2 / iterations) << "\n";
    plotSolution(bestRandomWalk2, points, distanceMatrix, "Random Walk 2");

    std::cout << "\nGreedy + Steepest Descent + Point Exchange:\n";
    std::cout << "Min: " << minGreedySteepestPointScore << " Max: " << maxGreedySteepestPointScore
              << " Avg: " << (totalGreedySteepestPointScore / iterations) << "\n";
    std::cout << "Min: " << minGreedySteepestPointTime << " Max: " << maxGreedySteepestPointTime
              << " Avg: " << (totalGreedySteepestPointTime / iterations) << "\n";
    plotSolution(bestGreedySteepestPointSolution, points, distanceMatrix, "Greedy + Steepest Descent + Vertex Exchange");

    std::cout << "\nGreedy + Steepest Descent + Edge Exchange:\n";
    std::cout << "Min: " << minGreedySteepestEdgeScore << " Max: " << maxGreedySteepestEdgeScore
              << " Avg: " << (totalGreedySteepestEdgeScore / iterations) << "\n";
    std::cout << "Min: " << minGreedySteepestEdgeTime << " Max: " << maxGreedySteepestEdgeTime
              << " Avg: " << (totalGreedySteepestEdgeTime / iterations) << "\n";
    plotSolution(bestGreedySteepestEdgeSolution, points, distanceMatrix, "Greedy + Steepest Descent + Edge Exchange");

    std::cout << "\nGreedy + Greedy Local Search + Point Exchange:\n";
    std::cout << "Min: " << minGreedyGreedyPointScore << " Max: " << maxGreedyGreedyPointScore
              << " Avg: " << (totalGreedyGreedyPointScore / iterations) << "\n";
    std::cout << "Min: " << minGreedyGreedyPointTime << " Max: " << maxGreedyGreedyPointTime
              << " Avg: " << (totalGreedyGreedyPointTime / iterations) << "\n";
    plotSolution(bestGreedyGreedyPointSolution, points, distanceMatrix, "Greedy + Greedy Local Search + Vertex Exchange");

    std::cout << "\nGreedy + Greedy Local Search + Edge Exchange:\n";
    std::cout << "Min: " << minGreedyGreedyEdgeScore << " Max: " << maxGreedyGreedyEdgeScore
              << " Avg: " << (totalGreedyGreedyEdgeScore / iterations) << "\n";
    std::cout << "Min: " << minGreedyGreedyEdgeTime << " Max: " << maxGreedyGreedyEdgeTime
              << " Avg: " << (totalGreedyGreedyEdgeTime / iterations) << "\n";
    plotSolution(bestGreedyGreedyEdgeSolution, points, distanceMatrix, "Greedy + Greedy Local Search + Edge Exchange");
}

void lab3(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix)
{
    int iterations = 100;

    // 0. Random + Steepest Descent + Edge Exchange
    double minRandomSteepestEdgeScore = std::numeric_limits<double>::max();
    double maxRandomSteepestEdgeScore = std::numeric_limits<double>::lowest();
    double totalRandomSteepestEdgeScore = 0.0;
    Solution bestRandomSteepestEdgeSolution;
    long minRandomSteepestEdgeTime = std::numeric_limits<long>::max();
    long maxRandomSteepestEdgeTime = std::numeric_limits<long>::lowest();
    long totalRandomSteepestEdgeTime = 0.0;

    // 1. Random + Priority Queue + Memory + Steepest Descent + Edge Exchange
    double minRandomMemorySteepestEdgeScore = std::numeric_limits<double>::max();
    double maxRandomMemorySteepestEdgeScore = std::numeric_limits<double>::lowest();
    double totalRandomMemorySteepestEdgeScore = 0.0;
    Solution bestRandomMemorySteepestEdgeSolution;
    long minRandomMemorySteepestEdgeTime = std::numeric_limits<long>::max();
    long maxRandomMemorySteepestEdgeTime = std::numeric_limits<long>::lowest();
    long totalRandomMemorySteepestEdgeTime = 0.0;

    // 2. Random + Steepest Descent + Candidate Moves
    double minRandomSteepestCandidatesScore = std::numeric_limits<double>::max();
    double maxRandomSteepestCandidatesScore = std::numeric_limits<double>::lowest();
    double totalRandomSteepestCandidatesScore = 0.0;
    Solution bestRandomSteepestCandidatesSolution;
    long minRandomSteepestCandidatesTime = std::numeric_limits<long>::max();
    long maxRandomSteepestCandidatesTime = std::numeric_limits<long>::lowest();
    long totalRandomSteepestCandidatesTime = 0.0;

    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    for (int i = 0; i < iterations; ++i)
    {
        std::cout << "Progress: " << i << "/" << iterations << std::endl;

        // Generate initial random and solution
        Solution initialRandom = randomCycle(distanceMatrix);

        // 0. Random + Steepest Descent + Edge Exchange
        auto t1 = std::chrono::high_resolution_clock::now();
        Solution solution0 = localSearchEdges(initialRandom, distanceMatrix, "steepest");
        auto t2 = std::chrono::high_resolution_clock::now();

        solution0.calculateScore(distanceMatrix);
        double score0 = solution0.getScore();
        minRandomSteepestEdgeScore = std::min(minRandomSteepestEdgeScore, score0);
        maxRandomSteepestEdgeScore = std::max(maxRandomSteepestEdgeScore, score0);
        totalRandomSteepestEdgeScore += score0;
        if (score0 == minRandomSteepestEdgeScore)
            bestRandomSteepestEdgeSolution = solution0;

        auto time0 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minRandomSteepestEdgeTime = std::min(minRandomSteepestEdgeTime, static_cast<long>(time0));
        maxRandomSteepestEdgeTime = std::max(maxRandomSteepestEdgeTime, static_cast<long>(time0));
        totalRandomSteepestEdgeTime += time0;

        // 1. Random + Priority Queue + Memory + Steepest Descent + Edge Exchange
        t1 = std::chrono::high_resolution_clock::now();
        Solution solution1 = localSearchMemory(initialRandom, distanceMatrix);
        t2 = std::chrono::high_resolution_clock::now();

        solution1.calculateScore(distanceMatrix);
        double score1 = solution1.getScore();
        minRandomMemorySteepestEdgeScore = std::min(minRandomMemorySteepestEdgeScore, score1);
        maxRandomMemorySteepestEdgeScore = std::max(maxRandomMemorySteepestEdgeScore, score1);
        totalRandomMemorySteepestEdgeScore += score1;
        if (score1 == minRandomMemorySteepestEdgeScore)
            bestRandomMemorySteepestEdgeSolution = solution1;

        auto time1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minRandomMemorySteepestEdgeTime = std::min(minRandomMemorySteepestEdgeTime, static_cast<long>(time1));
        maxRandomMemorySteepestEdgeTime = std::max(maxRandomMemorySteepestEdgeTime, static_cast<long>(time1));
        totalRandomMemorySteepestEdgeTime += time1;

        // 2. Random + Steepest Descent + Candidate Moves
        t1 = std::chrono::high_resolution_clock::now();
        Solution solution2 = localSearchCandidates(initialRandom, distanceMatrix);
        t2 = std::chrono::high_resolution_clock::now();

        solution2.calculateScore(distanceMatrix);
        double score2 = solution2.getScore();
        minRandomSteepestCandidatesScore = std::min(minRandomSteepestCandidatesScore, score2);
        maxRandomSteepestCandidatesScore = std::max(maxRandomSteepestCandidatesScore, score2);
        totalRandomSteepestCandidatesScore += score2;
        if (score2 == minRandomSteepestCandidatesScore)
            bestRandomSteepestCandidatesSolution = solution2;

        auto time2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minRandomSteepestCandidatesTime = std::min(minRandomSteepestCandidatesTime, static_cast<long>(time2));
        maxRandomSteepestCandidatesTime = std::max(maxRandomSteepestCandidatesTime, static_cast<long>(time2));
        totalRandomSteepestCandidatesTime += time2;
    }

    // Print results
    std::cout << "\nRandom + Steepest Descent + Edge Exchange:\n";
    std::cout << "Min: " << minRandomSteepestEdgeScore << " Max: " << maxRandomSteepestEdgeScore
              << " Avg: " << (totalRandomSteepestEdgeScore / iterations) << "\n";
    std::cout << "Min: " << minRandomSteepestEdgeTime << " Max: " << maxRandomSteepestEdgeTime
              << " Avg: " << (totalRandomSteepestEdgeTime / iterations) << "\n";
    plotSolution(bestRandomSteepestEdgeSolution, points, distanceMatrix, "Random + Steepest Descent + Edge Exchange");

    std::cout << "\nRandom + Memory + Steepest Descent + Edge Exchange:\n";
    std::cout << "Min: " << minRandomMemorySteepestEdgeScore << " Max: " << maxRandomMemorySteepestEdgeScore
              << " Avg: " << (totalRandomMemorySteepestEdgeScore / iterations) << "\n";
    std::cout << "Min: " << minRandomMemorySteepestEdgeTime << " Max: " << maxRandomMemorySteepestEdgeTime
              << " Avg: " << (totalRandomMemorySteepestEdgeTime / iterations) << "\n";
    plotSolution(bestRandomMemorySteepestEdgeSolution, points, distanceMatrix, "Random + Memory + Steepest Descent + Edge Exchange");

    std::cout << "\nRandom + Steepest Descent + Candidate Moves:\n";
    std::cout << "Min: " << minRandomSteepestCandidatesScore << " Max: " << maxRandomSteepestCandidatesScore
              << " Avg: " << (totalRandomSteepestCandidatesScore / iterations) << "\n";
    std::cout << "Min: " << minRandomSteepestCandidatesTime << " Max: " << maxRandomSteepestCandidatesTime
              << " Avg: " << (totalRandomSteepestCandidatesTime / iterations) << "\n";
    plotSolution(bestRandomSteepestCandidatesSolution, points, distanceMatrix, "Random + Steepest Descent + Candidate Moves");
}

void lab4(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix)
{
    int iterations = 10;

    // 1. Multiple start local search
    double minMultipleStartScore = std::numeric_limits<double>::max();
    double maxMultipleStartScore = std::numeric_limits<double>::lowest();
    double totalMultipleStartScore = 0.0;
    Solution bestMultipleStartSolution;
    long minMultipleStartTime = std::numeric_limits<long>::max();
    long maxMultipleStartTime = std::numeric_limits<long>::lowest();
    long totalMultipleStartTime = 0.0;

    // 2. Iterated local search
    double minIteratedScore = std::numeric_limits<double>::max();
    double maxIteratedScore = std::numeric_limits<double>::lowest();
    double totalIteratedScore = 0.0;
    Solution bestIteratedSolution;
    long minIteratedTime = std::numeric_limits<long>::max();
    long maxIteratedTime = std::numeric_limits<long>::lowest();
    long totalIteratedTime = 0.0;
    int minIteratedIter = std::numeric_limits<int>::max();
    int maxIteratedIter = std::numeric_limits<int>::lowest();
    int totalIteratedIter = 0.0;

    // 3. Large neighborhood search
    double minLargeScore = std::numeric_limits<double>::max();
    double maxLargeScore = std::numeric_limits<double>::lowest();
    double totalLargeScore = 0.0;
    Solution bestLargeSolution;
    long minLargeTime = std::numeric_limits<long>::max();
    long maxLargeTime = std::numeric_limits<long>::lowest();
    long totalLargeTime = 0.0;
    int minLargeIter = std::numeric_limits<int>::max();
    int maxLargeIter = std::numeric_limits<int>::lowest();
    int totalLargeIter = 0.0;

    for (int i = 0; i < iterations; ++i)
    {
        std::cout << "Progress: " << i << "/" << iterations << std::endl;

        // 1. Multiple start local search
        auto t1 = std::chrono::high_resolution_clock::now();
        Solution solution0 = multipleStartLocalSearch(distanceMatrix);
        auto t2 = std::chrono::high_resolution_clock::now();

        solution0.calculateScore(distanceMatrix);
        double score0 = solution0.getScore();
        minMultipleStartScore = std::min(minMultipleStartScore, score0);
        maxMultipleStartScore = std::max(maxMultipleStartScore, score0);
        totalMultipleStartScore += score0;
        if (score0 == minMultipleStartScore)
            bestMultipleStartSolution = solution0;

        auto time0 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minMultipleStartTime = std::min(minMultipleStartTime, static_cast<long>(time0));
        maxMultipleStartTime = std::max(maxMultipleStartTime, static_cast<long>(time0));
        totalMultipleStartTime += time0;

        // 2. Iterated local search
        t1 = std::chrono::high_resolution_clock::now();
        std::pair<Solution, int> result1 = iteratedLocalSearch(distanceMatrix, time0);
        t2 = std::chrono::high_resolution_clock::now();

        Solution solution1 = result1.first;
        int iter1 = result1.second;

        solution1.calculateScore(distanceMatrix);
        double score1 = solution1.getScore();
        minIteratedScore = std::min(minIteratedScore, score1);
        maxIteratedScore = std::max(maxIteratedScore, score1);
        totalIteratedScore += score1;
        if (score1 == minIteratedScore)
            bestIteratedSolution = solution1;

        auto time1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minIteratedTime = std::min(minIteratedTime, static_cast<long>(time1));
        maxIteratedTime = std::max(maxIteratedTime, static_cast<long>(time1));
        totalIteratedTime += time1;

        minIteratedIter = std::min(minIteratedIter, iter1);
        maxIteratedIter = std::max(maxIteratedIter, iter1);
        totalIteratedIter += iter1;

        // 3. Large neighborhood search
        t1 = std::chrono::high_resolution_clock::now();
        std::pair<Solution, int> result2 = largeNeighborhoodSearch(distanceMatrix, time0);
        t2 = std::chrono::high_resolution_clock::now();

        Solution solution2 = result2.first;
        int iter2 = result2.second;

        solution2.calculateScore(distanceMatrix);
        double score2 = solution2.getScore();
        minLargeScore = std::min(minLargeScore, score2);
        maxLargeScore = std::max(maxLargeScore, score2);
        totalLargeScore += score2;
        if (score2 == minLargeScore)
            bestLargeSolution = solution2;

        auto time2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        minLargeTime = std::min(minLargeTime, static_cast<long>(time2));
        maxLargeTime = std::max(maxLargeTime, static_cast<long>(time2));
        totalLargeTime += time2;

        minLargeIter = std::min(minLargeIter, iter2);
        maxLargeIter = std::max(maxLargeIter, iter2);
        totalLargeIter += iter2;
    }

    // Print results
    std::cout << "\nMultiple start local search:\n";
    std::cout << "Min: " << minMultipleStartScore << " Max: " << maxMultipleStartScore
              << " Avg: " << (totalMultipleStartScore / iterations) << "\n";
    std::cout << "Min: " << minMultipleStartTime << " Max: " << maxMultipleStartTime
              << " Avg: " << (totalMultipleStartTime / iterations) << "\n";
    plotSolution(bestMultipleStartSolution, points, distanceMatrix, "Multiple start local search");

    std::cout << "\nIterated local search:\n";
    std::cout << "Results\n";
    std::cout << "Min: " << minIteratedScore << " Max: " << maxIteratedScore
              << " Avg: " << (totalIteratedScore / iterations) << "\n";
    std::cout << "Time\n";
    std::cout << "Min: " << minIteratedTime << " Max: " << maxIteratedTime
              << " Avg: " << (totalIteratedTime / iterations) << "\n";
    std::cout << "Iterations\n";
    std::cout << "Min: " << minIteratedIter << " Max: " << maxIteratedIter
              << " Avg: " << (totalIteratedIter / iterations) << "\n";
    plotSolution(bestIteratedSolution, points, distanceMatrix, "Iterated local search");

    std::cout << "\nLarge neighborhood search:\n";
    std::cout << "Results\n";
    std::cout << "Min: " << minLargeScore << " Max: " << maxLargeScore
              << " Avg: " << (totalLargeScore / iterations) << "\n";
    std::cout << "Time\n";
    std::cout << "Min: " << minLargeTime << " Max: " << maxLargeTime
              << " Avg: " << (totalLargeTime / iterations) << "\n";
    std::cout << "Iterations\n";
    std::cout << "Min: " << minLargeIter << " Max: " << maxLargeIter
              << " Avg: " << (totalLargeIter / iterations) << "\n";
    plotSolution(bestLargeSolution, points, distanceMatrix, "Large neighborhood search");
}