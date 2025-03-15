#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <limits>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

struct Point
{
    double x, y;
};

struct Solution
{
    std::vector<int> cycle1Indexes;
    std::vector<int> cycle2Indexes;
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
    solution.cycle1Indexes.push_back(start1);
    solution.cycle2Indexes.push_back(start2);
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
            if (!visited[j] && distanceMatrix[current1][j] < minDistance && solution.cycle1Indexes.size() < maxCycleSize)
            {
                minDistance = distanceMatrix[current1][j];
                nextPoint = j;
                addToCycle = 1;
            }
            if (!visited[j] && distanceMatrix[current2][j] < minDistance && solution.cycle2Indexes.size() < maxCycleSize)
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
            solution.cycle1Indexes.push_back(nextPoint);
            current1 = nextPoint;
        }
        else
        {
            solution.cycle2Indexes.push_back(nextPoint);
            current2 = nextPoint;
        }
    }

    return solution;
}

// Greedy cycle algorithm for solving TSP with insertion at the best position
Solution greedyCycle(const std::vector<std::vector<double>> &distanceMatrix)
{
    Solution solution;
    int n = distanceMatrix.size();

    // Randomly select the starting point for the 1st cycle
    int start1 = std::rand() % n;
    int start2 = findFurthestPointIndex(distanceMatrix, start1);

    // Initialize the cycles with the starting points
    solution.cycle1Indexes.push_back(start1);
    solution.cycle2Indexes.push_back(start2);
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
            if (!visited[j] && solution.cycle1Indexes.size() < maxCycleSize)
            {
                // Try inserting in cycle 1
                for (int k = 0; k < solution.cycle1Indexes.size(); ++k)
                {
                    int prev = solution.cycle1Indexes[k];
                    int next = solution.cycle1Indexes[(k + 1) % solution.cycle1Indexes.size()];
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

            if (!visited[j] && solution.cycle2Indexes.size() < maxCycleSize)
            {
                // Try inserting in cycle 2
                for (int k = 0; k < solution.cycle2Indexes.size(); ++k)
                {
                    int prev = solution.cycle2Indexes[k];
                    int next = solution.cycle2Indexes[(k + 1) % solution.cycle2Indexes.size()];
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
            solution.cycle1Indexes.insert(solution.cycle1Indexes.begin() + insertPosition, nextPoint);
            current1 = nextPoint;
        }
        else
        {
            solution.cycle2Indexes.insert(solution.cycle2Indexes.begin() + insertPosition, nextPoint);
            current2 = nextPoint;
        }
    }

    return solution;
}

// Regret-based greedy cycle algorithm for TSP
Solution regretCycleWeighted(const std::vector<std::vector<double>> &distanceMatrix, double weightCost = 0.0, double weightRegret = 1.0)
{
    Solution solution;
    int n = distanceMatrix.size();

    // Randomly select the starting points for the two cycles
    int start1 = std::rand() % n;
    int start2 = findFurthestPointIndex(distanceMatrix, start1);

    // Initialize cycles with starting points
    solution.cycle1Indexes.push_back(start1);
    solution.cycle2Indexes.push_back(start2);
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
            if (solution.cycle1Indexes.size() < maxCycleSize)
            {
                for (size_t idx = 0; idx < solution.cycle1Indexes.size(); ++idx)
                {
                    int prev = solution.cycle1Indexes[idx];
                    int next = solution.cycle1Indexes[(idx + 1) % solution.cycle1Indexes.size()];
                    double cost = distanceMatrix[prev][j] + distanceMatrix[j][next] - distanceMatrix[prev][next];

                    if (cost < bestCost1)
                    {
                        bestCost1 = cost;
                        bestPos1 = idx + 1;
                    }
                }
            }

            // Find best insertion cost for cycle 2
            if (solution.cycle2Indexes.size() < maxCycleSize)
            {
                for (size_t idx = 0; idx < solution.cycle2Indexes.size(); ++idx)
                {
                    int prev = solution.cycle2Indexes[idx];
                    int next = solution.cycle2Indexes[(idx + 1) % solution.cycle2Indexes.size()];
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
            solution.cycle1Indexes.insert(solution.cycle1Indexes.begin() + bestInsertPos, nextPoint);
        }
        else
        {
            solution.cycle2Indexes.insert(solution.cycle2Indexes.begin() + bestInsertPos, nextPoint);
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

double calculateCycleDistance(const std::vector<int> &cycle, const std::vector<std::vector<double>> &distanceMatrix)
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

void plotSolution(const Solution &solution, const std::vector<Point> &points, const std::vector<std::vector<double>> &distanceMatrix)
{
    // Calculate distances
    double cycle1Distance = calculateCycleDistance(solution.cycle1Indexes, distanceMatrix);
    double cycle2Distance = calculateCycleDistance(solution.cycle2Indexes, distanceMatrix);
    double totalDistance = cycle1Distance + cycle2Distance;

    // Print cycles and distances to console
    std::cout << "Cycle 1 (size " << solution.cycle1Indexes.size() << "): ";
    for (int index : solution.cycle1Indexes)
    {
        std::cout << index << " ";
    }
    std::cout << "\nTotal Distance: " << cycle1Distance << std::endl;

    std::cout << "Cycle 2 (size " << solution.cycle2Indexes.size() << "): ";
    for (int index : solution.cycle2Indexes)
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
    for (size_t i = 0; i < solution.cycle1Indexes.size(); ++i)
    {
        int a = solution.cycle1Indexes[i];
        int b = solution.cycle1Indexes[(i + 1) % solution.cycle1Indexes.size()];
        plt::plot({points[a].x, points[b].x}, {points[a].y, points[b].y}, "r");
    }

    // Plot Cycle 2 in Blue
    for (size_t i = 0; i < solution.cycle2Indexes.size(); ++i)
    {
        int a = solution.cycle2Indexes[i];
        int b = solution.cycle2Indexes[(i + 1) % solution.cycle2Indexes.size()];
        plt::plot({points[a].x, points[b].x}, {points[a].y, points[b].y}, "b");
    }

    // Display score under the plot
    plt::text(0.05, -0.1, "Total Distance: " + std::to_string(totalDistance));

    // Show plot
    plt::show();
}

double calculateSolutionScore(const Solution &solution, const std::vector<std::vector<double>> &distanceMatrix)
{
    double cycle1Distance = calculateCycleDistance(solution.cycle1Indexes, distanceMatrix);
    double cycle2Distance = calculateCycleDistance(solution.cycle2Indexes, distanceMatrix);
    return cycle1Distance + cycle2Distance;
}

int main()
{
    std::string filepath = "TSPlib95/kroB200.tsp";
    std::vector<Point> points = loadPointsFromFile(filepath);
    if (points.empty())
    {
        std::cerr << "No points loaded from the file!" << std::endl;
        return -1;
    }

    // Create the adjacency list
    std::vector<std::vector<double>> distanceMatrix = createdistanceMatrix(points);

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

    // RegretCycleWeighted Algorithm
    double minRegretCycleScore = std::numeric_limits<double>::max();
    double maxRegretCycleScore = std::numeric_limits<double>::lowest();
    double totalRegretCycleScore = 0.0;
    Solution bestRegretCycleSolution;

    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    for (int i = 0; i < iterations; ++i)
    {
        // Apply greedyNearestNeighbour algorithm
        Solution solution1 = greedyNearestNeighbour(distanceMatrix);
        double greedyNNScore = calculateSolutionScore(solution1, distanceMatrix);
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
        double greedyCycleScore = calculateSolutionScore(solution2, distanceMatrix);
        minGreedyCycleScore = std::min(minGreedyCycleScore, greedyCycleScore);
        maxGreedyCycleScore = std::max(maxGreedyCycleScore, greedyCycleScore);
        totalGreedyCycleScore += greedyCycleScore;

        // Check for best solution for GreedyCycle
        if (greedyCycleScore == minGreedyCycleScore)
        {
            minGreedyCycleScore = greedyCycleScore;
            bestGreedyCycleSolution = solution2;
        }

        // Apply regretCycleWeighted algorithm
        Solution solution3 = regretCycleWeighted(distanceMatrix, 0.0, 1.0);
        double regretCycleScore = calculateSolutionScore(solution3, distanceMatrix);
        minRegretCycleScore = std::min(minRegretCycleScore, regretCycleScore);
        maxRegretCycleScore = std::max(maxRegretCycleScore, regretCycleScore);
        totalRegretCycleScore += regretCycleScore;

        // Check for best solution for RegretCycleWeighted
        if (regretCycleScore == minRegretCycleScore)
        {
            minRegretCycleScore = regretCycleScore;
            bestRegretCycleSolution = solution3;
        }
    }

    // Calculate average scores for each algorithm
    double avgGreedyNNScore = totalGreedyNNScore / iterations;
    double avgGreedyCycleScore = totalGreedyCycleScore / iterations;
    double avgRegretCycleScore = totalRegretCycleScore / iterations;

    // Print results for each algorithm
    std::cout << "\nGreedyNearestNeighbour Algorithm Results:" << std::endl;
    std::cout << "Min Score: " << minGreedyNNScore << std::endl;
    std::cout << "Max Score: " << maxGreedyNNScore << std::endl;
    std::cout << "Average Score: " << avgGreedyNNScore << std::endl;

    std::cout << "\nGreedyCycle Algorithm Results:" << std::endl;
    std::cout << "Min Score: " << minGreedyCycleScore << std::endl;
    std::cout << "Max Score: " << maxGreedyCycleScore << std::endl;
    std::cout << "Average Score: " << avgGreedyCycleScore << std::endl;

    std::cout << "\nRegretCycleWeighted Algorithm Results:" << std::endl;
    std::cout << "Min Score: " << minRegretCycleScore << std::endl;
    std::cout << "Max Score: " << maxRegretCycleScore << std::endl;
    std::cout << "Average Score: " << avgRegretCycleScore << std::endl;

    // Plot best solutions for each algorithm
    std::cout << "\nPlotting Best Solutions:" << std::endl;

    std::cout << "Best GreedyNearestNeighbour Solution:" << std::endl;
    plotSolution(bestGreedyNNSolution, points, distanceMatrix);

    std::cout << "Best GreedyCycle Solution:" << std::endl;
    plotSolution(bestGreedyCycleSolution, points, distanceMatrix);

    std::cout << "Best RegretCycleWeighted Solution:" << std::endl;
    plotSolution(bestRegretCycleSolution, points, distanceMatrix);

    return 0;
}