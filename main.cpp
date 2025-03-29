#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <limits>
#include <random>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

struct Point
{
    double x, y;
};

struct Solution
{
    std::vector<int> cycleIndices[2];
    double cycle1Score;
    double cycle2Score;
    double score;

    Solution() : cycle1Score(0.0), cycle2Score(0.0), score(0.0) {}

    void calculateScore(const std::vector<std::vector<double>> &distanceMatrix)
    {
        cycle1Score = calculateCycleDistance(cycleIndices[0], distanceMatrix);
        cycle2Score = calculateCycleDistance(cycleIndices[1], distanceMatrix);
        score = cycle1Score + cycle2Score;
    }

    void modifyScore(double change)
    {
        score += change;
    }

    double getScore() const
    {
        return score;
    }

    double getCycle1Score(const std::vector<std::vector<double>> &distanceMatrix)
    {
        cycle1Score = calculateCycleDistance(cycleIndices[0], distanceMatrix);
        return cycle1Score;
    }

    double getCycle2Score(const std::vector<std::vector<double>> &distanceMatrix)
    {
        cycle2Score = calculateCycleDistance(cycleIndices[1], distanceMatrix);
        return cycle2Score;
    }

    double calculateCycleDistance(const std::vector<int> &cycle, const std::vector<std::vector<double>> &distanceMatrix) const
    {
        double totalDistance = 0.0;
        int cycleSize = cycle.size();
        for (int i = 0; i < cycleSize; ++i)
        {
            int current = cycle[i];
            int next = cycle[(i + 1) % cycleSize]; // Wrap around to form a cycle
            totalDistance += distanceMatrix[current][next];
        }
        return totalDistance;
    }
};

std::vector<Point> loadPointsFromFile(const std::string &filepath)
{
    std::vector<Point> points;
    std::ifstream file(filepath);
    if (!file.is_open())
    {
        std::cerr << "Error opening file!" << std::endl;
        return points;
    }
    std::string line;
    bool nodeSection = false;
    // Read the file line by line
    while (std::getline(file, line))
    {
        // Skip lines until we reach NODE_COORD_SECTION
        if (line.find("NODE_COORD_SECTION") != std::string::npos)
        {
            nodeSection = true;
            continue;
        }

        // Once we're in the NODE_COORD_SECTION, parse the points
        if (nodeSection)
        {
            // If we reach EOF or "EOF" line, stop parsing
            if (line.find("EOF") != std::string::npos)
            {
                break;
            }

            // Read the x and y values from the line
            int node_id;
            double x, y;
            std::istringstream ss(line);
            ss >> node_id >> x >> y;

            // Add point to the list
            points.push_back({x, y});
        }
    }

    file.close();
    return points;
}

double euclideanDistance(Point a, Point b)
{
    return std::sqrt(std::pow(b.x - a.x, 2) + std::pow(b.y - a.y, 2));
}

std::vector<std::vector<double>> createdistanceMatrix(const std::vector<Point> &points)
{
    int n = points.size();
    std::vector<std::vector<double>> distanceMatrix(n, std::vector<double>(n, 0.0));
    // Calculate the Euclidean distance for each pair of points
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            double dist = euclideanDistance(points[i], points[j]);
            distanceMatrix[i][j] = dist;
            distanceMatrix[j][i] = dist; // Symmetrical distance
        }
    }
    return distanceMatrix;
}

void plotPoints(const std::vector<Point> &points)
{
    std::vector<double> x_coords, y_coords;
    for (const auto &point : points)
    {
        x_coords.push_back(point.x);
        y_coords.push_back(point.y);
    }
    plt::scatter(x_coords, y_coords, 10.0);
    plt::show();
}

int findFurthestPointIndex(const std::vector<std::vector<double>> &distanceMatrix, int start)
{
    int n = distanceMatrix.size();
    int furthestPointIndex = -1;
    double maxDistance = -1;

    // Iterate over all points to find the one furthest from the start
    for (int i = 0; i < n; ++i)
    {
        if (i != start)
        {
            double distance = distanceMatrix[start][i];
            if (distance > maxDistance)
            {
                maxDistance = distance;
                furthestPointIndex = i;
            }
        }
    }

    return furthestPointIndex;
}

Solution greedyNearestNeighbour(const std::vector<std::vector<double>> &distanceMatrix)
{
    Solution solution;
    int n = distanceMatrix.size();

    // Randomly select the starting points for the 1st cycle
    int start1 = std::rand() % n;
    int start2 = findFurthestPointIndex(distanceMatrix, start1);

    // Greedily build the first cycle starting from start1 and start2
    solution.cycleIndices[0].push_back(start1);
    solution.cycleIndices[1].push_back(start2);
    std::vector<bool> visited(n, false);
    visited[start1] = true;
    visited[start2] = true;
    int current1 = start1;
    int current2 = start2;
    // If there is an odd number of points then in 1 cycle there must be 1 vertex more
    // Example: 3 vertices = 2 + 1
    // And (3+1)/2 = 2
    // That's why +1
    int maxCycleSize = (n + 1) / 2;
    // 2 vertices are already visited so we check n-2 times
    for (int i = 2; i < n; ++i)
    {
        double minDistance = std::numeric_limits<double>::max();
        int nextPoint = -1;
        int addToCycle = -1;

        // Find the nearest unvisited neighbor using the distanceMatrix
        for (int j = 0; j < n; ++j)
        {
            if (!visited[j] && distanceMatrix[current1][j] < minDistance && solution.cycleIndices[0].size() < maxCycleSize)
            {
                minDistance = distanceMatrix[current1][j];
                nextPoint = j;
                addToCycle = 1;
            }
            if (!visited[j] && distanceMatrix[current2][j] < minDistance && solution.cycleIndices[1].size() < maxCycleSize)
            {
                minDistance = distanceMatrix[current2][j];
                nextPoint = j;
                addToCycle = 2;
            }
        }
        // Add best vertex to the appropriate cycle
        visited[nextPoint] = true;
        if (addToCycle == 1)
        {
            solution.cycleIndices[0].push_back(nextPoint);
            current1 = nextPoint;
        }
        else
        {
            solution.cycleIndices[1].push_back(nextPoint);
            current2 = nextPoint;
        }
    }

    return solution;
}

Solution greedyCycle(const std::vector<std::vector<double>> &distanceMatrix)
{
    Solution solution;
    int n = distanceMatrix.size();

    // Randomly select the starting point for the 1st cycle
    int start1 = std::rand() % n;
    int start2 = findFurthestPointIndex(distanceMatrix, start1);

    // Initialize the cycles with the starting points
    solution.cycleIndices[0].push_back(start1);
    solution.cycleIndices[1].push_back(start2);
    std::vector<bool> visited(n, false);
    visited[start1] = true;
    visited[start2] = true;

    int current1 = start1;
    int current2 = start2;
    int maxCycleSize = (n + 1) / 2;

    // Greedily build the cycles, inserting the best vertex at the best position
    for (int i = 2; i < n; ++i)
    {
        double minIncrease = std::numeric_limits<double>::max();
        int nextPoint = -1;
        int addToCycle = -1;
        int insertPosition = -1;

        // Try to insert the new vertex in all possible positions in both cycles
        for (int j = 0; j < n; ++j)
        {
            if (!visited[j] && solution.cycleIndices[0].size() < maxCycleSize)
            {
                // Try inserting in cycle 1
                for (int k = 0; k < solution.cycleIndices[0].size(); ++k)
                {
                    int prev = solution.cycleIndices[0][k];
                    int next = solution.cycleIndices[0][(k + 1) % solution.cycleIndices[0].size()];
                    double increase = distanceMatrix[prev][j] + distanceMatrix[j][next] - distanceMatrix[prev][next];

                    if (increase < minIncrease)
                    {
                        minIncrease = increase;
                        nextPoint = j;
                        addToCycle = 1;
                        insertPosition = k + 1;
                    }
                }
            }

            if (!visited[j] && solution.cycleIndices[1].size() < maxCycleSize)
            {
                // Try inserting in cycle 2
                for (int k = 0; k < solution.cycleIndices[1].size(); ++k)
                {
                    int prev = solution.cycleIndices[1][k];
                    int next = solution.cycleIndices[1][(k + 1) % solution.cycleIndices[1].size()];
                    double increase = distanceMatrix[prev][j] + distanceMatrix[j][next] - distanceMatrix[prev][next];

                    if (increase < minIncrease)
                    {
                        minIncrease = increase;
                        nextPoint = j;
                        addToCycle = 2;
                        insertPosition = k + 1;
                    }
                }
            }
        }

        // Insert the best vertex into the chosen cycle at the determined position
        visited[nextPoint] = true;
        if (addToCycle == 1)
        {
            solution.cycleIndices[0].insert(solution.cycleIndices[0].begin() + insertPosition, nextPoint);
            current1 = nextPoint;
        }
        else
        {
            solution.cycleIndices[1].insert(solution.cycleIndices[1].begin() + insertPosition, nextPoint);
            current2 = nextPoint;
        }
    }

    return solution;
}

Solution regretCycleWeighted(const std::vector<std::vector<double>> &distanceMatrix, double weightCost = 0.0, double weightRegret = 1.0)
{
    Solution solution;
    int n = distanceMatrix.size();

    // Randomly select the starting points for the two cycles
    int start1 = std::rand() % n;
    int start2 = findFurthestPointIndex(distanceMatrix, start1);

    // Initialize cycles with starting points
    solution.cycleIndices[0].push_back(start1);
    solution.cycleIndices[1].push_back(start2);
    std::vector<bool> visited(n, false);
    visited[start1] = true;
    visited[start2] = true;

    int maxCycleSize = (n + 1) / 2;

    // Greedily build the cycles using regret heuristic
    for (int i = 2; i < n; ++i)
    {
        int nextPoint = -1;
        int bestCycle = -1;
        int bestInsertPos = -1;
        double minScore = std::numeric_limits<double>::max();

        // Try inserting each unvisited point
        for (int j = 0; j < n; ++j)
        {
            if (visited[j])
                continue;

            double bestCost1 = std::numeric_limits<double>::max();
            double bestCost2 = std::numeric_limits<double>::max();
            int bestPos1 = -1, bestPos2 = -1;

            // Find best insertion cost for cycle 1
            if (solution.cycleIndices[0].size() < maxCycleSize)
            {
                for (size_t idx = 0; idx < solution.cycleIndices[0].size(); ++idx)
                {
                    int prev = solution.cycleIndices[0][idx];
                    int next = solution.cycleIndices[0][(idx + 1) % solution.cycleIndices[0].size()];
                    double cost = distanceMatrix[prev][j] + distanceMatrix[j][next] - distanceMatrix[prev][next];

                    if (cost < bestCost1)
                    {
                        bestCost1 = cost;
                        bestPos1 = idx + 1;
                    }
                }
            }

            // Find best insertion cost for cycle 2
            if (solution.cycleIndices[1].size() < maxCycleSize)
            {
                for (size_t idx = 0; idx < solution.cycleIndices[1].size(); ++idx)
                {
                    int prev = solution.cycleIndices[1][idx];
                    int next = solution.cycleIndices[1][(idx + 1) % solution.cycleIndices[1].size()];
                    double cost = distanceMatrix[prev][j] + distanceMatrix[j][next] - distanceMatrix[prev][next];

                    if (cost < bestCost2)
                    {
                        bestCost2 = cost;
                        bestPos2 = idx + 1;
                    }
                }
            }

            // Compute regret (difference between best insertions in cycle1 and cycle2)
            double score = std::min(bestCost1, bestCost2) * weightCost - std::abs(bestCost1 - bestCost2) * weightRegret;

            // Choose the city with the highest regret
            if (score < minScore)
            {
                minScore = score;
                nextPoint = j;

                if (bestCost1 < bestCost2)
                {
                    bestCycle = 1;
                    bestInsertPos = bestPos1;
                }
                else
                {
                    bestCycle = 2;
                    bestInsertPos = bestPos2;
                }
            }
        }

        // Insert the selected city at the best position in the chosen cycle
        visited[nextPoint] = true;
        if (bestCycle == 1)
        {
            solution.cycleIndices[0].insert(solution.cycleIndices[0].begin() + bestInsertPos, nextPoint);
        }
        else
        {
            solution.cycleIndices[1].insert(solution.cycleIndices[1].begin() + bestInsertPos, nextPoint);
        }
    }

    return solution;
}

void displayDistanceMatrix(const std::vector<std::vector<double>> &distanceMatrix)
{
    std::cout << "Distance Matrix (distances between points):\n";
    for (int i = 0; i < distanceMatrix.size(); ++i)
    {
        for (int j = 0; j < distanceMatrix[i].size(); ++j)
        {
            std::cout << distanceMatrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void plotSolution(Solution &solution, const std::vector<Point> &points, const std::vector<std::vector<double>> &distanceMatrix)
{
    // Calculate distances
    double cycle1Distance = solution.getCycle1Score(distanceMatrix);
    double cycle2Distance = solution.getCycle2Score(distanceMatrix);
    double totalDistance = cycle1Distance + cycle2Distance;

    // Print cycles and distances to console
    std::cout << "Cycle 1 (size " << solution.cycleIndices[0].size() << "): ";
    for (int index : solution.cycleIndices[0])
    {
        std::cout << index << " ";
    }
    std::cout << "\nTotal Distance: " << cycle1Distance << std::endl;

    std::cout << "Cycle 2 (size " << solution.cycleIndices[1].size() << "): ";
    for (int index : solution.cycleIndices[1])
    {
        std::cout << index << " ";
    }
    std::cout << "\nTotal Distance: " << cycle2Distance << std::endl;

    std::cout << "Total Distance for both cycles: " << totalDistance << std::endl;

    // Prepare data for plotting
    std::vector<double> x_coords, y_coords;
    for (const auto &point : points)
    {
        x_coords.push_back(point.x);
        y_coords.push_back(point.y);
    }

    // Scatter plot of points
    plt::scatter(x_coords, y_coords, 10.0);

    // Plot Cycle 1 in Red
    for (size_t i = 0; i < solution.cycleIndices[0].size(); ++i)
    {
        int a = solution.cycleIndices[0][i];
        int b = solution.cycleIndices[0][(i + 1) % solution.cycleIndices[0].size()];
        plt::plot({points[a].x, points[b].x}, {points[a].y, points[b].y}, "r");
    }

    // Plot Cycle 2 in Blue
    for (size_t i = 0; i < solution.cycleIndices[1].size(); ++i)
    {
        int a = solution.cycleIndices[1][i];
        int b = solution.cycleIndices[1][(i + 1) % solution.cycleIndices[1].size()];
        plt::plot({points[a].x, points[b].x}, {points[a].y, points[b].y}, "b");
    }

    // Display score under the plot
    plt::text(0.05, -0.1, "Total Distance: " + std::to_string(totalDistance));

    // Show plot
    plt::show();
}

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
        Solution solution3 = regretCycleWeighted(distanceMatrix, 0.0, 1.0);
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
        Solution solution4 = regretCycleWeighted(distanceMatrix, 1.0, 1.0);
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
    std::cout << "\nGreedyNearestNeighbour Algorithm Results:" << std::endl;
    std::cout << "Min Score: " << minGreedyNNScore << std::endl;
    std::cout << "Max Score: " << maxGreedyNNScore << std::endl;
    std::cout << "Average Score: " << avgGreedyNNScore << std::endl;

    std::cout << "\nGreedyCycle Algorithm Results:" << std::endl;
    std::cout << "Min Score: " << minGreedyCycleScore << std::endl;
    std::cout << "Max Score: " << maxGreedyCycleScore << std::endl;
    std::cout << "Average Score: " << avgGreedyCycleScore << std::endl;

    std::cout << "\nRegretCycle Algorithm Results:" << std::endl;
    std::cout << "Min Score: " << minRegretCycleScore << std::endl;
    std::cout << "Max Score: " << maxRegretCycleScore << std::endl;
    std::cout << "Average Score: " << avgRegretCycleScore << std::endl;

    std::cout << "\nRegretCycleWeighted Algorithm Results:" << std::endl;
    std::cout << "Min Score: " << minWeightedRegretCycleScore << std::endl;
    std::cout << "Max Score: " << maxWeightedRegretCycleScore << std::endl;
    std::cout << "Average Score: " << avgWeightedRegretCycleScore << std::endl;

    // Plot best solutions for each algorithm
    std::cout << "\nPlotting Best Solutions:" << std::endl;

    std::cout << "Best GreedyNearestNeighbour Solution:" << std::endl;
    plotSolution(bestGreedyNNSolution, points, distanceMatrix);

    std::cout << "Best GreedyCycle Solution:" << std::endl;
    plotSolution(bestGreedyCycleSolution, points, distanceMatrix);

    std::cout << "Best RegretCycle Solution:" << std::endl;
    plotSolution(bestRegretCycleSolution, points, distanceMatrix);

    std::cout << "Best RegretCycleWeighted Solution:" << std::endl;
    plotSolution(bestWeightedRegretCycleSolution, points, distanceMatrix);
}

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
            if (delta < minDelta)
            {
                besti = i;
                bestj = j;
                minDelta = delta;
                foundSwap = true;
            }
            // But greedy LS means take 1st found
            if (strategy == "greedy")
            {
                if (delta < 0)
                {
                    // A swap was found
                    solution.cycleIndices[cycleOfPoint1][i % pointsInCycle1] = point2;
                    solution.cycleIndices[cycleOfPoint2][j % pointsInCycle1] = point1;
                    solution.calculateScore(distanceMatrix);
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

Solution localSearchVertex(Solution solution, const std::vector<std::vector<double>> &distanceMatrix, const std::string &strategy)
{
    bool improvement = true;

    while (improvement)
    {
        improvement = stepPointExchange(solution, distanceMatrix, strategy);
    }
    return solution;
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

    // 1. Random + Steepest Descent + Point Exchange
    double minRandomSteepestPointScore = std::numeric_limits<double>::max();
    double maxRandomSteepestPointScore = std::numeric_limits<double>::lowest();
    double totalRandomSteepestPointScore = 0.0;
    Solution bestRandomSteepestPointSolution;

    // 2. Random + Steepest Descent + Edge Exchange
    double minRandomSteepestEdgeScore = std::numeric_limits<double>::max();
    double maxRandomSteepestEdgeScore = std::numeric_limits<double>::lowest();
    double totalRandomSteepestEdgeScore = 0.0;
    Solution bestRandomSteepestEdgeSolution;

    // 3. Random + Greedy Local Search + Point Exchange
    double minRandomGreedyPointScore = std::numeric_limits<double>::max();
    double maxRandomGreedyPointScore = std::numeric_limits<double>::lowest();
    double totalRandomGreedyPointScore = 0.0;
    Solution bestRandomGreedyPointSolution;

    // 4. Random + Greedy Local Search + Edge Exchange
    double minRandomGreedyEdgeScore = std::numeric_limits<double>::max();
    double maxRandomGreedyEdgeScore = std::numeric_limits<double>::lowest();
    double totalRandomGreedyEdgeScore = 0.0;
    Solution bestRandomGreedyEdgeSolution;

    // 5. Greedy + Steepest Descent + Point Exchange
    double minGreedySteepestPointScore = std::numeric_limits<double>::max();
    double maxGreedySteepestPointScore = std::numeric_limits<double>::lowest();
    double totalGreedySteepestPointScore = 0.0;
    Solution bestGreedySteepestPointSolution;

    // 6. Greedy + Steepest Descent + Edge Exchange
    double minGreedySteepestEdgeScore = std::numeric_limits<double>::max();
    double maxGreedySteepestEdgeScore = std::numeric_limits<double>::lowest();
    double totalGreedySteepestEdgeScore = 0.0;
    Solution bestGreedySteepestEdgeSolution;

    // 7. Greedy + Greedy Local Search + Point Exchange
    double minGreedyGreedyPointScore = std::numeric_limits<double>::max();
    double maxGreedyGreedyPointScore = std::numeric_limits<double>::lowest();
    double totalGreedyGreedyPointScore = 0.0;
    Solution bestGreedyGreedyPointSolution;

    // 8. Greedy + Greedy Local Search + Edge Exchange
    double minGreedyGreedyEdgeScore = std::numeric_limits<double>::max();
    double maxGreedyGreedyEdgeScore = std::numeric_limits<double>::lowest();
    double totalGreedyGreedyEdgeScore = 0.0;
    Solution bestGreedyGreedyEdgeSolution;

    for (int i = 0; i < iterations; ++i)
    {
        std::cout << i << std::endl;
        // Generate initial random and greedy solutions
        Solution initialRandom = randomCycle(distanceMatrix);
        Solution initialGreedy = regretCycleWeighted(distanceMatrix, 1.0, 1.0);

        // 1. Random + Steepest Descent + Point Exchange
        Solution solution1 = localSearchVertex(initialRandom, distanceMatrix, "steepest");
        solution1.calculateScore(distanceMatrix);
        double score1 = solution1.getScore();
        minRandomSteepestPointScore = std::min(minRandomSteepestPointScore, score1);
        maxRandomSteepestPointScore = std::max(maxRandomSteepestPointScore, score1);
        totalRandomSteepestPointScore += score1;
        if (score1 == minRandomSteepestPointScore)
            bestRandomSteepestPointSolution = solution1;

        // 2. Random + Steepest Descent + Edge Exchange (Placeholder)
        Solution solution2; // Implement Edge Exchange here
        // solution2.calculateScore(distanceMatrix);
        // double score2 = solution2.getScore();
        // minRandomSteepestEdgeScore = std::min(minRandomSteepestEdgeScore, score2);
        // maxRandomSteepestEdgeScore = std::max(maxRandomSteepestEdgeScore, score2);
        // totalRandomSteepestEdgeScore += score2;
        // if (score2 == minRandomSteepestEdgeScore) bestRandomSteepestEdgeSolution = solution2;

        // 3. Random + Greedy Local Search + Point Exchange
        Solution solution3 = localSearchVertex(initialRandom, distanceMatrix, "greedy");
        solution3.calculateScore(distanceMatrix);
        double score3 = solution3.getScore();
        minRandomGreedyPointScore = std::min(minRandomGreedyPointScore, score3);
        maxRandomGreedyPointScore = std::max(maxRandomGreedyPointScore, score3);
        totalRandomGreedyPointScore += score3;
        if (score3 == minRandomGreedyPointScore)
            bestRandomGreedyPointSolution = solution3;

        // 4. Random + Greedy Local Search + Edge Exchange (Placeholder)
        Solution solution4; // Implement Edge Exchange here

        // 5. Greedy + Steepest Descent + Point Exchange
        Solution solution5 = localSearchVertex(initialGreedy, distanceMatrix, "steepest");
        solution5.calculateScore(distanceMatrix);
        double score5 = solution5.getScore();
        minGreedySteepestPointScore = std::min(minGreedySteepestPointScore, score5);
        maxGreedySteepestPointScore = std::max(maxGreedySteepestPointScore, score5);
        totalGreedySteepestPointScore += score5;
        if (score5 == minGreedySteepestPointScore)
            bestGreedySteepestPointSolution = solution5;

        // 6. Greedy + Steepest Descent + Edge Exchange (Placeholder)
        Solution solution6; // Implement Edge Exchange here

        // 7. Greedy + Greedy Local Search + Point Exchange
        Solution solution7 = localSearchVertex(initialGreedy, distanceMatrix, "greedy");
        solution7.calculateScore(distanceMatrix);
        double score7 = solution7.getScore();
        minGreedyGreedyPointScore = std::min(minGreedyGreedyPointScore, score7);
        maxGreedyGreedyPointScore = std::max(maxGreedyGreedyPointScore, score7);
        totalGreedyGreedyPointScore += score7;
        if (score7 == minGreedyGreedyPointScore)
            bestGreedyGreedyPointSolution = solution7;

        // 8. Greedy + Greedy Local Search + Edge Exchange (Placeholder)
        Solution solution8; // Implement Edge Exchange here
    }

    // Print results
    std::cout << "\nRandom + Steepest Descent + Point Exchange:\n";
    std::cout << "Min: " << minRandomSteepestPointScore << " Max: " << maxRandomSteepestPointScore
              << " Avg: " << (totalRandomSteepestPointScore / iterations) << "\n";
    plotSolution(bestRandomSteepestPointSolution, points, distanceMatrix);

    std::cout << "\nRandom + Steepest Descent + Edge Exchange:\n";
    std::cout << "Min: " << minRandomSteepestEdgeScore << " Max: " << maxRandomSteepestEdgeScore
              << " Avg: " << (totalRandomSteepestEdgeScore / iterations) << "\n";
    plotSolution(bestRandomSteepestEdgeSolution, points, distanceMatrix);

    std::cout << "\nRandom + Greedy Local Search + Point Exchange:\n";
    std::cout << "Min: " << minRandomGreedyPointScore << " Max: " << maxRandomGreedyPointScore
              << " Avg: " << (totalRandomGreedyPointScore / iterations) << "\n";
    plotSolution(bestRandomGreedyPointSolution, points, distanceMatrix);

    std::cout << "\nRandom + Greedy Local Search + Edge Exchange:\n";
    std::cout << "Min: " << minRandomGreedyEdgeScore << " Max: " << maxRandomGreedyEdgeScore
              << " Avg: " << (totalRandomGreedyEdgeScore / iterations) << "\n";
    plotSolution(bestRandomGreedyEdgeSolution, points, distanceMatrix);

    std::cout << "\nGreedy + Steepest Descent + Point Exchange:\n";
    std::cout << "Min: " << minGreedySteepestPointScore << " Max: " << maxGreedySteepestPointScore
              << " Avg: " << (totalGreedySteepestPointScore / iterations) << "\n";
    plotSolution(bestGreedySteepestPointSolution, points, distanceMatrix);

    std::cout << "\nGreedy + Steepest Descent + Edge Exchange:\n";
    std::cout << "Min: " << minGreedySteepestEdgeScore << " Max: " << maxGreedySteepestEdgeScore
              << " Avg: " << (totalGreedySteepestEdgeScore / iterations) << "\n";
    plotSolution(bestGreedySteepestEdgeSolution, points, distanceMatrix);

    std::cout << "\nGreedy + Greedy Local Search + Point Exchange:\n";
    std::cout << "Min: " << minGreedyGreedyPointScore << " Max: " << maxGreedyGreedyPointScore
              << " Avg: " << (totalGreedyGreedyPointScore / iterations) << "\n";
    plotSolution(bestGreedyGreedyPointSolution, points, distanceMatrix);

    std::cout << "\nGreedy + Greedy Local Search + Edge Exchange:\n";
    std::cout << "Min: " << minGreedyGreedyEdgeScore << " Max: " << maxGreedyGreedyEdgeScore
              << " Avg: " << (totalGreedyGreedyEdgeScore / iterations) << "\n";
    plotSolution(bestGreedyGreedyEdgeSolution, points, distanceMatrix);
}

int main(int argc, char *argv[])
{
    std::string filepath = argv[1];
    std::vector<Point> points = loadPointsFromFile(filepath);
    if (points.empty())
    {
        std::cerr << "No points loaded from the file!" << std::endl;
        return -1;
    }

    // Create the adjacency list
    std::vector<std::vector<double>> distanceMatrix = createdistanceMatrix(points);

    lab2(points, distanceMatrix);

    return 0;
}