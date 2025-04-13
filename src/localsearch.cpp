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
#include "localsearch.h"

Solution randomCycle(const std::vector<std::vector<double>> &distanceMatrix)
{
    Solution solution;
    int n = distanceMatrix.size();

    // Generate a shuffled list of indices
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::mt19937(std::random_device()()));

    // Create a vector with half 0s and half 1s
    std::vector<int> cycleAssignment(n, 0);
    std::fill(cycleAssignment.begin() + n / 2, cycleAssignment.end(), 1);
    std::shuffle(cycleAssignment.begin(), cycleAssignment.end(), std::mt19937(std::random_device()()));

    // Assign points to cycles based on the cycleAssignment vector
    for (int i = 0; i < n; ++i)
    {
        if (cycleAssignment[i] == 0)
            solution.cycleIndices[0].push_back(indices[i]);
        else
            solution.cycleIndices[1].push_back(indices[i]);
    }

    return solution;
}

bool stepPointExchange(Solution &solution, const std::vector<std::vector<double>> &distanceMatrix, const std::string &strategy)
{
    const double epsilon = 1e-6; // Avoid numerical errors
    size_t numberOfPoints = distanceMatrix.size();
    size_t pointsInCycle1 = solution.cycleIndices[0].size(); // Size of first cycle
    size_t pointsInCycle2 = solution.cycleIndices[1].size(); // Size of second cycle

    std::vector<size_t> cycleSizes = {pointsInCycle1, pointsInCycle2};

    double minDelta = 0.0;
    int besti = -1;
    int bestj = -1;
    bool foundSwap = false;

    int point1;
    int point2;
    int cycleOfPoint1;
    int cycleOfPoint2;

    // Loop over each pair of points in the solution
    for (size_t i = 0; i < numberOfPoints; ++i)
    {
        for (size_t j = i + 1; j < numberOfPoints; ++j)
        {
            // Calculate to which cycle these point belong to
            // Example 13/50 = 0 (integer division) so it is the 1st cycle
            // Example 63/50 = 1 (integer division) so it is the 2nd cycle
            cycleOfPoint1 = i / pointsInCycle1;
            cycleOfPoint2 = j / pointsInCycle1;
            // Calculate how much can we gain by swapping these 2 points
            double delta = 0;
            // Get all the neighbours of point 1
            int next1;
            int prev1;
            point1 = solution.cycleIndices[cycleOfPoint1][i % pointsInCycle1];
            next1 = solution.cycleIndices[cycleOfPoint1][(i + 1) % cycleSizes[cycleOfPoint1]];
            prev1 = solution.cycleIndices[cycleOfPoint1][(i - 1 + cycleSizes[cycleOfPoint1]) % cycleSizes[cycleOfPoint1]];
            // Get all the neighbours of point 2
            int next2;
            int prev2;
            point2 = solution.cycleIndices[cycleOfPoint2][j % pointsInCycle1];
            next2 = solution.cycleIndices[cycleOfPoint2][(j + 1) % cycleSizes[cycleOfPoint2]];
            prev2 = solution.cycleIndices[cycleOfPoint2][(j - 1 + cycleSizes[cycleOfPoint2]) % cycleSizes[cycleOfPoint2]];
            // Calculate delta
            delta -= distanceMatrix[point1][next1];
            delta -= distanceMatrix[point1][prev1];
            delta -= distanceMatrix[point2][next2];
            delta -= distanceMatrix[point2][prev2];
            // Handle edge case when we have A----->point1----->point2----->B
            // Then distance[point1][prev2] is distance[point1][point1] and it should be distance[point1][point2]
            if (point2 == next1)
            {
                prev2 = point2;
                next1 = point1;
            }
            // Handle edge case when we have A----->point2----->point1----->B
            // Then distance[point1][next2] is distance[point1][point1] and it should be distance[point1][point2]
            if (point2 == prev1)
            {
                prev1 = point1;
                next2 = point2;
            }
            delta += distanceMatrix[point1][next2];
            delta += distanceMatrix[point1][prev2];
            delta += distanceMatrix[point2][next1];
            delta += distanceMatrix[point2][prev1];

            // Check if this is the best swap so far
            if (delta < minDelta - epsilon)
            {
                besti = i;
                bestj = j;
                minDelta = delta;
                foundSwap = true;
            }
            // But greedy LS means take 1st found
            if (strategy == "greedy" && delta < -epsilon)
            {
                // A swap was found
                solution.cycleIndices[cycleOfPoint1][i % pointsInCycle1] = point2;
                solution.cycleIndices[cycleOfPoint2][j % pointsInCycle1] = point1;
                return true;
            }
        }
    }
    if (!foundSwap)
    {
        return false;
    }
    else
    {
        // A swap was found
        cycleOfPoint1 = besti / pointsInCycle1;
        cycleOfPoint2 = bestj / pointsInCycle1;
        point1 = solution.cycleIndices[cycleOfPoint1][besti % pointsInCycle1];
        point2 = solution.cycleIndices[cycleOfPoint2][bestj % pointsInCycle1];
        solution.cycleIndices[cycleOfPoint1][besti % pointsInCycle1] = point2;
        solution.cycleIndices[cycleOfPoint2][bestj % pointsInCycle1] = point1;
        return true;
    }

    return foundSwap;
}

bool stepEdgeExchange(Solution &solution, const std::vector<std::vector<double>> &distanceMatrix, const std::string &strategy)
{
    const double epsilon = 1e-6; // Avoid numerical error
    size_t numberOfPoints = distanceMatrix.size();
    size_t pointsInCycle1 = solution.cycleIndices[0].size();
    size_t pointsInCycle2 = solution.cycleIndices[1].size();
    std::vector<size_t> cycleSizes = {pointsInCycle1, pointsInCycle2};

    double minDelta = 0.0;
    int besti = -1;
    int bestj = -1;
    int bestc = -1;
    bool foundSwap = false;

    int point1, point2, cycleOfPoint1, cycleOfPoint2;

    // Loop over each cycle
    for (size_t c = 0; c < 2; ++c)
    {
        // Loop over each pair of points in the solution
        for (size_t i = 0; i < cycleSizes[c]; ++i)
        {
            for (size_t j = 1; j < cycleSizes[c] - 1; ++j) // Avoid reversing whole cycle
            {
                double delta = 0;
                int prev1, next2;
                point1 = solution.cycleIndices[c][i];
                prev1 = solution.cycleIndices[c][(i - 1 + cycleSizes[c]) % cycleSizes[c]];
                point2 = solution.cycleIndices[c][(i + j) % cycleSizes[c]];
                next2 = solution.cycleIndices[c][(i + j + 1) % cycleSizes[c]];

                delta -= distanceMatrix[point1][prev1];
                delta -= distanceMatrix[point2][next2];
                delta += distanceMatrix[point1][next2];
                delta += distanceMatrix[point2][prev1];

                if (delta < minDelta - epsilon)
                {
                    bestc = c;
                    besti = i;
                    bestj = j;
                    minDelta = delta;
                    foundSwap = true;
                }

                if (strategy == "greedy" && delta < -epsilon)
                {
                    // Perform the swap immediately
                    size_t endIdx = (i + j) % cycleSizes[c];

                    if (endIdx > i)
                    {
                        // Normal case: reverse directly
                        std::reverse(solution.cycleIndices[c].begin() + i, solution.cycleIndices[c].begin() + endIdx + 1);
                    }
                    else
                    {
                        // Wrap-around case
                        std::vector<int> temp;

                        // Move front elements to the back
                        temp.insert(temp.end(), solution.cycleIndices[c].begin() + endIdx + 1, solution.cycleIndices[c].end());
                        temp.insert(temp.end(), solution.cycleIndices[c].begin(), solution.cycleIndices[c].begin() + endIdx + 1);

                        // Reverse the required section
                        std::reverse(temp.begin(), temp.begin() + (cycleSizes[c] - i));

                        // Restore back to the original vector
                        for (size_t k = 0; k < temp.size(); ++k)
                        {
                            solution.cycleIndices[c][(i + k) % cycleSizes[c]] = temp[k];
                        }
                    }
                    return true;
                }
            }
        }
    }

    if (!foundSwap)
    {
        return false;
    }
    else
    {
        // Perform the best swap found
        size_t endIdx = (besti + bestj) % cycleSizes[bestc];

        if (endIdx > besti)
        {
            // Normal case: reverse directly
            std::reverse(solution.cycleIndices[bestc].begin() + besti, solution.cycleIndices[bestc].begin() + endIdx + 1);
        }
        else
        {
            // Wrap-around case
            std::vector<int> temp;

            // Move front elements to the back
            temp.insert(temp.end(), solution.cycleIndices[bestc].begin() + endIdx + 1, solution.cycleIndices[bestc].end());
            temp.insert(temp.end(), solution.cycleIndices[bestc].begin(), solution.cycleIndices[bestc].begin() + endIdx + 1);

            // Reverse the required section
            std::reverse(temp.begin(), temp.begin() + (cycleSizes[bestc] - besti));

            // Restore back to the original vector
            for (size_t k = 0; k < temp.size(); ++k)
            {
                solution.cycleIndices[bestc][(besti + k) % cycleSizes[bestc]] = temp[k];
            }
            return true;
        }
    }
    return foundSwap;
}

Solution localSearchVertex(Solution solution, const std::vector<std::vector<double>> &distanceMatrix, const std::string &strategy)
{
    if (strategy == "greedy")
    {
        // Randomize the process slightly
        std::random_device rd;
        std::mt19937 rng(rd());

        // Shuffle starting points for both cycles
        for (int i = 0; i < 2; ++i)
        {
            if (!solution.cycleIndices[i].empty())
            {
                std::uniform_int_distribution<int> dist(0, solution.cycleIndices[i].size() - 1);
                int randomIndex = dist(rng); // Pick a random index
                // Rotate cycle to start from the chosen random index
                std::rotate(solution.cycleIndices[i].begin(),
                            solution.cycleIndices[i].begin() + randomIndex,
                            solution.cycleIndices[i].end());
            }
        }
    }

    bool improvement = true;
    while (improvement)
    {
        improvement = stepPointExchange(solution, distanceMatrix, strategy);
    }
    return solution;
}

Solution localSearchEdges(Solution solution, const std::vector<std::vector<double>> &distanceMatrix, const std::string &strategy)
{
    if (strategy == "greedy")
    {
        // Randomize the process slightly
        std::random_device rd;
        std::mt19937 rng(rd());

        // Shuffle starting points for both cycles
        for (int i = 0; i < 2; ++i)
        {
            if (!solution.cycleIndices[i].empty())
            {
                std::uniform_int_distribution<int> dist(0, solution.cycleIndices[i].size() - 1);
                int randomIndex = dist(rng); // Pick a random index
                // Rotate cycle to start from the chosen random index
                std::rotate(solution.cycleIndices[i].begin(),
                            solution.cycleIndices[i].begin() + randomIndex,
                            solution.cycleIndices[i].end());
            }
        }
    }

    bool improvement = true;
    while (improvement)
    {
        improvement = stepEdgeExchange(solution, distanceMatrix, strategy);
    }
    return solution;
}

Solution randomWalk(Solution &solution, const std::vector<std::vector<double>> &distanceMatrix, int timeLimit = 1000)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    Solution bestSolution = solution;
    solution.calculateScore(distanceMatrix);
    double bestScore = solution.score;

    Solution currentSolution = solution;

    size_t pointsInCycle1 = solution.cycleIndices[0].size();
    size_t pointsInCycle2 = solution.cycleIndices[1].size();
    std::vector<size_t> cycleSizes = {pointsInCycle1, pointsInCycle2};

    auto t1 = std::chrono::high_resolution_clock::now();
    long runtime = 0;
    while(runtime < timeLimit)
    {
        Solution tempSolution = currentSolution;
        std::uniform_int_distribution<int> choiceDist(0, 1);
        bool swapVertices = choiceDist(gen); // Randomly choose to swap vertices or edges

        if (swapVertices)
        {
            // Swap two random vertices
            std::uniform_int_distribution<int> pointDist(0, pointsInCycle1 + pointsInCycle2 - 1);
            int index1 = pointDist(gen);
            int index2 = pointDist(gen);

            while (index1 == index2) // Ensure different indices
                index2 = pointDist(gen);

            int cycle1 = (index1 < pointsInCycle1) ? 0 : 1;
            int cycle2 = (index2 < pointsInCycle1) ? 0 : 1;
            index1 %= cycleSizes[cycle1];
            index2 %= cycleSizes[cycle2];

            std::swap(tempSolution.cycleIndices[cycle1][index1], tempSolution.cycleIndices[cycle2][index2]);
        }
        else
        {
            // Swap edges
            std::uniform_int_distribution<int> cycleDist(0, 1);
            int c = cycleDist(gen);
            std::uniform_int_distribution<int> indexDist(0, cycleSizes[c] - 1);

            int i = indexDist(gen);
            int j = indexDist(gen);
            while (j == 0 || j == cycleSizes[c] - 1) // Ensure valid range
                j = indexDist(gen);

            size_t endIdx = (i + j) % cycleSizes[c];

            if (endIdx > i)
            {
                // Normal case: reverse directly
                std::reverse(tempSolution.cycleIndices[c].begin() + i, tempSolution.cycleIndices[c].begin() + endIdx + 1);
            }
            else
            {
                // Wrap-around case
                std::vector<int> temp;

                // Move front elements to the back
                temp.insert(temp.end(), tempSolution.cycleIndices[c].begin() + endIdx + 1, tempSolution.cycleIndices[c].end());
                temp.insert(temp.end(), tempSolution.cycleIndices[c].begin(), tempSolution.cycleIndices[c].begin() + endIdx + 1);

                // Reverse the required section
                std::reverse(temp.begin(), temp.begin() + (cycleSizes[c] - i));

                // Restore back to the original vector
                for (size_t k = 0; k < temp.size(); ++k)
                {
                    tempSolution.cycleIndices[c][(i + k) % cycleSizes[c]] = temp[k];
                }
            }
        }

        tempSolution.calculateScore(distanceMatrix);
        double newScore = tempSolution.score;
        if (newScore < bestScore)
        {
            bestSolution = tempSolution;
            bestScore = newScore;
        }

        currentSolution = tempSolution; // Update current solution for the next iteration
        auto t2 = std::chrono::high_resolution_clock::now();
        runtime = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
    }

    return bestSolution;
}