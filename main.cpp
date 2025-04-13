#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <limits>
#include <random>
#include <chrono>
#include <queue>
#include <algorithm>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

struct Point
{
    double x, y;
};

std::vector<Point> points;

struct Solution
{
    std::vector<int> cycleIndices[2];
    double cycle1Score;
    double cycle2Score;
    double score;

    // pointPositions[i] = {cycleIndex, positionInCycle}
    std::vector<std::pair<int, int>> pointPositions;

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

    void updatePointPositions()
    {
        int totalPoints = cycleIndices[0].size() + cycleIndices[1].size();
        pointPositions.clear();
        pointPositions.resize(totalPoints, {-1, -1});

        for (int cycle = 0; cycle < 2; ++cycle)
        {
            for (int pos = 0; pos < cycleIndices[cycle].size(); ++pos)
            {
                int point = cycleIndices[cycle][pos];
                pointPositions[point] = {cycle, pos};
            }
        }
    }

    std::pair<int, int> getPointPosition(int pointIndex) const
    {
        return pointPositions[pointIndex];
    }
};

enum moveType
{
    VERTEX = 0,
    EDGE = 1
};

struct VertexNeighbourhood
{
    int prev;
    int mid;
    int next;

    VertexNeighbourhood(int p, int m, int n)
        : prev(p), mid(m), next(n) {}

    bool operator==(const VertexNeighbourhood &other) const
    {
        return mid == other.mid &&
               ((prev == other.prev && next == other.next) ||
                (prev == other.next && next == other.prev));
    }

    bool operator!=(const VertexNeighbourhood &other) const
    {
        return !(*this == other);
    }
};

struct Cut
{
    int prev;
    int next;

    Cut() : prev(-1), next(-1) {}
    Cut(int p, int n)
        : prev(p), next(n) {}

    bool operator==(const Cut &other) const
    {
        return (prev == other.prev && next == other.next) ||
               (prev == other.next && next == other.prev);
    }

    bool operator!=(const Cut &other) const
    {
        return !(*this == other);
    }

    bool isReversedOf(const Cut &other) const
    {
        return (prev == other.next && next == other.prev);
    }

    bool isSameDirectionAs(const Cut &other) const
    {
        return (prev == other.prev && next == other.next);
    }
};

struct Move
{
    moveType type;
    double delta;

    // Optional fields depending on type
    VertexNeighbourhood neighbourhood1{0, 0, 0}, neighbourhood2{0, 0, 0};
    Cut cut1{0, 0}, cut2{0, 0};

    Move() : type(VERTEX), delta(0.0) {}
    Move(moveType t, double d) : type(t), delta(d) {}

    // Sorting by delta (min-heap priority queue)
    bool operator<(const Move &other) const
    {
        return delta > other.delta;
    }
};

struct VertexMove : public Move
{
    VertexMove(const VertexNeighbourhood &v1,
               const VertexNeighbourhood &v2,
               double d)
        : Move(VERTEX, d)
    {
        neighbourhood1 = v1;
        neighbourhood2 = v2;
    }
};

struct EdgeMove : public Move
{
    EdgeMove(const Cut &c1,
             const Cut &c2,
             double d)
        : Move(EDGE, d)
    {
        cut1 = c1;
        cut2 = c2;
    }

    bool operator==(const EdgeMove &other) const
    {
        // Case 1: same order
        if (cut1 == other.cut1 && cut2 == other.cut2)
        {
            return ((cut1.isSameDirectionAs(other.cut1) && cut2.isSameDirectionAs(other.cut2)) || (cut1.isReversedOf(other.cut1) && cut2.isReversedOf(other.cut2)));
        }

        // Case 2: reversed order
        if (cut1 == other.cut2 && cut2 == other.cut1)
        {
            return ((cut1.isSameDirectionAs(other.cut2) && cut2.isSameDirectionAs(other.cut1)) || (cut1.isReversedOf(other.cut2) && cut2.isReversedOf(other.cut1)));
        }

        return false;
    }

    bool operator!=(const EdgeMove &other) const
    {
        return !(*this == other);
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

        // Find the nearest unvisited neighbour using the distanceMatrix
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
    solution.updatePointPositions();
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
    solution.updatePointPositions();
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
    solution.updatePointPositions();
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

void plotSolution(Solution &solution, const std::vector<Point> &points, const std::vector<std::vector<double>> &distanceMatrix, const std::string &title)
{
    // Calculate distances
    double cycle1Distance = solution.getCycle1Score(distanceMatrix);
    double cycle2Distance = solution.getCycle2Score(distanceMatrix);
    double totalDistance = cycle1Distance + cycle2Distance;

    // // Print cycles and distances to console
    // std::cout << "Cycle 1 (size " << solution.cycleIndices[0].size() << "): ";
    // for (int index : solution.cycleIndices[0])
    // {
    //     std::cout << index << " ";
    // }
    // std::cout << "\nTotal Distance: " << cycle1Distance << std::endl;

    // std::cout << "Cycle 2 (size " << solution.cycleIndices[1].size() << "): ";
    // for (int index : solution.cycleIndices[1])
    // {
    //     std::cout << index << " ";
    // }
    // std::cout << "\nTotal Distance: " << cycle2Distance << std::endl;

    // std::cout << "Total Distance for both cycles: " << totalDistance << std::endl;

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

    // Add title and display score under the plot
    plt::title(title);
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

Solution randomCycle(const std::vector<std::vector<double>> &distanceMatrix)
{
    Solution solution;
    int n = distanceMatrix.size();

    // int seed = 41; // dowolna liczba całkowita jako seed
    // std::mt19937 rng(seed);

    // Generate a shuffled list of indices
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::mt19937(std::random_device()()));
    // std::shuffle(indices.begin(), indices.end(), rng);

    // Create a vector with half 0s and half 1s
    std::vector<int> cycleAssignment(n, 0);
    std::fill(cycleAssignment.begin() + n / 2, cycleAssignment.end(), 1);
    std::shuffle(cycleAssignment.begin(), cycleAssignment.end(), std::mt19937(std::random_device()()));
    // std::shuffle(cycleAssignment.begin(), cycleAssignment.end(), rng);

    // Assign points to cycles based on the cycleAssignment vector
    for (int i = 0; i < n; ++i)
    {
        if (cycleAssignment[i] == 0)
            solution.cycleIndices[0].push_back(indices[i]);
        else
            solution.cycleIndices[1].push_back(indices[i]);
    }
    solution.updatePointPositions();
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
    int bestPos1 = -1;
    int bestPos2 = -1;
    int bestCycle = -1;
    bool foundSwap = false;
    bool swappingVertices = true;

    int point1, point2;
    int cycleOfPoint1 = 0;
    int cycleOfPoint2 = 1;

    Move bestMove;
    // Add vertex swap between 2 cycles
    // For some reason it should always be a legal move
    for (size_t i = 0; i < pointsInCycle1; ++i)
    {
        for (size_t j = 0; j < pointsInCycle2; ++j)
        {
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
            delta += distanceMatrix[point1][next2];
            delta += distanceMatrix[point1][prev2];
            delta += distanceMatrix[point2][next1];
            delta += distanceMatrix[point2][prev1];
            // Only add moves that can improve the solution
            if (strategy == "greedy" && delta < -epsilon)
            {
                int pos1 = i;
                int pos2 = j;
                // Ensure the points are in separate cycles
                if (cycleOfPoint1 != cycleOfPoint2)
                {
                    std::swap(solution.cycleIndices[cycleOfPoint1][pos1], solution.cycleIndices[cycleOfPoint2][pos2]);
                    std::swap(solution.pointPositions[point1], solution.pointPositions[point2]);
                }
                return true;
            }
            if (delta < minDelta - epsilon)
            {
                bestPos1 = i;
                bestPos2 = j;
                minDelta = delta;
                foundSwap = true;
            }
        }
    }
    for (size_t c = 0; c < 2; ++c)
    {
        // Loop over each pair of points in the solution
        for (size_t i = 0; i < cycleSizes[c]; ++i)
        {
            for (size_t j = 0; j < cycleSizes[c]; ++j) // Avoid reversing whole cycle
            {
                // Skip reversing the entire cycle
                if (i == j || (j + 1) % cycleSizes[c] == i || (j + 2) % cycleSizes[c] == i)
                {
                    continue;
                }
                double delta = 0;
                int prev1, next2;
                point1 = solution.cycleIndices[c][i];
                prev1 = solution.cycleIndices[c][(i - 1 + cycleSizes[c]) % cycleSizes[c]];
                point2 = solution.cycleIndices[c][j];
                next2 = solution.cycleIndices[c][(j + 1) % cycleSizes[c]];

                delta -= distanceMatrix[point1][prev1];
                delta -= distanceMatrix[point2][next2];
                delta += distanceMatrix[point1][next2];
                delta += distanceMatrix[point2][prev1];

                if (delta < minDelta - epsilon)
                {
                    bestMove = EdgeMove({prev1, point1}, {point2, next2}, delta);
                    bestCycle = c;
                    bestPos1 = i;
                    bestPos2 = j;
                    minDelta = delta;
                    foundSwap = true;
                    swappingVertices = false; // Change strategy
                }

                if (strategy == "greedy" && delta < -epsilon)
                {
                    bestPos1 = i;
                    bestPos2 = j;
                    if (bestPos1 < bestPos2)
                    {
                        // Normal case: reverse directly
                        std::reverse(solution.cycleIndices[bestCycle].begin() + bestPos1, solution.cycleIndices[bestCycle].begin() + bestPos2 + 1);
                        solution.updatePointPositions();
                    }
                    else
                    {
                        // Wrap-around case
                        std::reverse(solution.cycleIndices[bestCycle].begin() + bestPos2 + 1, solution.cycleIndices[bestCycle].begin() + bestPos1);
                        solution.updatePointPositions();
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
        // Check which type of a move to do
        if (swappingVertices)
        {
            //  Swap the points between the two cycles
            int pointToSwap1 = solution.cycleIndices[0][bestPos1];
            int pointToSwap2 = solution.cycleIndices[1][bestPos2];
            // std::cout << "Move vertices " << pointToSwap1 << " " << pointToSwap2 << " score " << minDelta << std::endl;
            std::swap(solution.cycleIndices[0][bestPos1], solution.cycleIndices[1][bestPos2]);
            std::swap(solution.pointPositions[pointToSwap1], solution.pointPositions[pointToSwap2]);
        }
        else
        {
            // Get current positions of the points in the cycles
            auto [cycle1, pos1] = solution.getPointPosition(bestMove.cut1.next);
            auto [cycle2, pos2] = solution.getPointPosition(bestMove.cut2.prev);
            int c = cycle1; // Both cycles are the same
            // Ensure position 1 is before position 2
            if (pos1 > pos2)
            {
                std::swap(pos1, pos2);
                std::swap(bestMove.cut1, bestMove.cut2);
            }
            std::reverse(solution.cycleIndices[c].begin() + pos1, solution.cycleIndices[c].begin() + pos2 + 1);
            solution.updatePointPositions();
            // std::cout << ": (" << bestMove.cut1.prev << " " << bestMove.cut1.next << ") and (" << bestMove.cut2.prev << " " << bestMove.cut2.next << ")" << " score " << bestMove.delta << std::endl;
        }
        return true;
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
    solution.updatePointPositions();
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
    solution.updatePointPositions();
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
    bestSolution.updatePointPositions();
    return bestSolution;
}

void lab2(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix)
{
    int iterations = 100;
    /*
        In this experiment, we evaluate different approaches to solving the 2-cycle TSP.
        We have three independent choices, leading to 2^3 = 8 possible solution strategies:

        1. **Starting Solution**
           - `Random` (Completely randomized initial solution)
           - `Greedy` (Constructed using a heuristic like Nearest neighbour or Greedy Cycle)

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

    for (int i = 0; i < iterations; ++i)
    {
        std::cout << "Progress: " << i << "/" << iterations << std::endl;

        // Generate initial random and greedy solutions
        Solution initialRandom = randomCycle(distanceMatrix);
        Solution initialGreedy = regretCycleWeighted(distanceMatrix, 1.0, 1.0);

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

        auto time0 = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
        minRandomWalk1time = std::min(minRandomWalk1time, time0);
        maxRandomWalk1time = std::max(maxRandomWalk1time, time0);
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

        auto time1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
        minRandomSteepestPointTime = std::min(minRandomSteepestPointTime, time1);
        maxRandomSteepestPointTime = std::max(maxRandomSteepestPointTime, time1);
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

        auto time2 = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
        minRandomSteepestEdgeTime = std::min(minRandomSteepestEdgeTime, time2);
        maxRandomSteepestEdgeTime = std::max(maxRandomSteepestEdgeTime, time2);
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

        auto time3 = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
        minRandomGreedyPointTime = std::min(minRandomGreedyPointTime, time3);
        maxRandomGreedyPointTime = std::max(maxRandomGreedyPointTime, time3);
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

        auto time4 = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
        minRandomGreedyEdgeTime = std::min(minRandomGreedyEdgeTime, time4);
        maxRandomGreedyEdgeTime = std::max(maxRandomGreedyEdgeTime, time4);
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

        auto time5 = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();

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

        auto time6 = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
        minGreedySteepestPointTime = std::min(minGreedySteepestPointTime, time6);
        maxGreedySteepestPointTime = std::max(maxGreedySteepestPointTime, time6);
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

        auto time7 = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
        minGreedySteepestEdgeTime = std::min(minGreedySteepestEdgeTime, time7);
        maxGreedySteepestEdgeTime = std::max(maxGreedySteepestEdgeTime, time7);
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

        auto time8 = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
        minGreedyGreedyPointTime = std::min(minGreedyGreedyPointTime, time8);
        maxGreedyGreedyPointTime = std::max(maxGreedyGreedyPointTime, time8);
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

        auto time9 = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
        minGreedyGreedyEdgeTime = std::min(minGreedyGreedyEdgeTime, time9);
        maxGreedyGreedyEdgeTime = std::max(maxGreedyGreedyEdgeTime, time9);
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

void enqueueImprovingVertexMoves(
    const std::vector<int> &pointsToUpdate,
    const Solution &solution,
    const std::vector<std::vector<double>> &distanceMatrix,
    const std::vector<size_t> &cycleSizes,
    std::priority_queue<Move> &pq,
    double epsilon)
{
    for (int updatedPoint : pointsToUpdate)
    {
        auto [cycle, position] = solution.getPointPosition(updatedPoint);
        int otherCycle = 1 - cycle;

        for (int j = 0; j < cycleSizes[otherCycle]; j++)
        {
            double delta = 0.0;

            int otherPoint = solution.cycleIndices[otherCycle][j];

            int newNext1 = solution.cycleIndices[cycle][(position + 1) % cycleSizes[cycle]];
            int newPrev1 = solution.cycleIndices[cycle][(position - 1 + cycleSizes[cycle]) % cycleSizes[cycle]];

            int newNext2 = solution.cycleIndices[otherCycle][(j + 1) % cycleSizes[otherCycle]];
            int newPrev2 = solution.cycleIndices[otherCycle][(j - 1 + cycleSizes[otherCycle]) % cycleSizes[otherCycle]];

            delta -= distanceMatrix[updatedPoint][newNext1];
            delta -= distanceMatrix[updatedPoint][newPrev1];
            delta -= distanceMatrix[otherPoint][newNext2];
            delta -= distanceMatrix[otherPoint][newPrev2];

            delta += distanceMatrix[updatedPoint][newNext2];
            delta += distanceMatrix[updatedPoint][newPrev2];
            delta += distanceMatrix[otherPoint][newNext1];
            delta += distanceMatrix[otherPoint][newPrev1];

            if (delta < -epsilon)
            {
                pq.push(VertexMove(
                    {newPrev1, updatedPoint, newNext1},
                    {newPrev2, otherPoint, newNext2},
                    delta));
            }
        }
    }
}

void enqueueImprovingEdgeMoves(
    const std::vector<Cut> &cutsToUpdate,
    const Solution &solution,
    const std::vector<std::vector<double>> &distanceMatrix,
    const std::vector<size_t> &cycleSizes,
    std::priority_queue<Move> &pq,
    double epsilon)
{
    for (size_t i = 0; i < cycleSizes[0]; i++)
    {
        for (const Cut &cut : cutsToUpdate)
        {
            auto [cycleOfCut, positionOfCut] = solution.getPointPosition(cut.prev);

            // Skip reversing the entire cycle
            if (i == positionOfCut ||
                (i + 1) % cycleSizes[cycleOfCut] == positionOfCut ||
                (positionOfCut + 1) % cycleSizes[cycleOfCut] == i)
            {
                continue;
            }

            int point1 = solution.cycleIndices[cycleOfCut][i];
            int next1 = solution.cycleIndices[cycleOfCut][(i + 1) % cycleSizes[cycleOfCut]];
            int point2 = solution.cycleIndices[cycleOfCut][positionOfCut];
            int next2 = solution.cycleIndices[cycleOfCut][(positionOfCut + 1) % cycleSizes[cycleOfCut]];

            double delta = 0.0;
            delta -= distanceMatrix[point1][next1];
            delta -= distanceMatrix[point2][next2];
            delta += distanceMatrix[next1][next2];
            delta += distanceMatrix[point1][point2];

            if (delta < -epsilon)
            {
                pq.push(EdgeMove({point1, next1}, {point2, next2}, delta));
            }
        }
    }
}

Solution localSearchMemory(Solution solution, const std::vector<std::vector<double>> &distanceMatrix)
{
    std::priority_queue<Move> pq;
    // Push all the legal moves to the queue
    const double epsilon = 1e-6; // Avoid numerical errors
    size_t numberOfPoints = distanceMatrix.size();
    size_t pointsInCycle1 = solution.cycleIndices[0].size(); // Size of first cycle
    size_t pointsInCycle2 = solution.cycleIndices[1].size(); // Size of second cycle
    std::vector<size_t> cycleSizes = {pointsInCycle1, pointsInCycle2};

    int point1;
    int point2;
    int cycleOfPoint1 = 0;
    int cycleOfPoint2 = 1;
    double delta;

    // Add vertex swap between 2 cycles
    // For some reason it should always be a legal move
    for (size_t i = 0; i < pointsInCycle1; ++i)
    {
        for (size_t j = 0; j < pointsInCycle2; ++j)
        {
            delta = 0;
            // Get all the neighbours of point 1
            int next1;
            int prev1;
            point1 = solution.cycleIndices[cycleOfPoint1][i % pointsInCycle1];
            next1 = solution.cycleIndices[cycleOfPoint1][(i + 1) % cycleSizes[cycleOfPoint1]];
            prev1 = solution.cycleIndices[cycleOfPoint1][(i - 1 + cycleSizes[cycleOfPoint1]) % cycleSizes[cycleOfPoint1]];
            // Get all the neighbours of point 2
            int next2;
            int prev2;
            point2 = solution.cycleIndices[cycleOfPoint2][j % pointsInCycle2];
            next2 = solution.cycleIndices[cycleOfPoint2][(j + 1) % cycleSizes[cycleOfPoint2]];
            prev2 = solution.cycleIndices[cycleOfPoint2][(j - 1 + cycleSizes[cycleOfPoint2]) % cycleSizes[cycleOfPoint2]];
            // Calculate delta
            delta -= distanceMatrix[point1][next1];
            delta -= distanceMatrix[point1][prev1];
            delta -= distanceMatrix[point2][next2];
            delta -= distanceMatrix[point2][prev2];
            delta += distanceMatrix[point1][next2];
            delta += distanceMatrix[point1][prev2];
            delta += distanceMatrix[point2][next1];
            delta += distanceMatrix[point2][prev1];
            // Only add moves that can improve the solution
            if (delta < -epsilon)
            {
                pq.push(VertexMove({prev1, point1, next1}, {prev2, point2, next2}, delta));
            }
        }
    }
    // Add edge exchange
    for (size_t c = 0; c < 2; ++c)
    {
        // Loop over each pair of points in the solution
        for (size_t i = 0; i < cycleSizes[c]; ++i)
        {
            for (size_t j = 0; j < cycleSizes[c]; ++j) // Avoid reversing whole cycle
            {
                // Skip reversing the entire cycle
                if (i == j || (j + 1) % cycleSizes[c] == i || (j + 2) % cycleSizes[c] == i)
                {
                    continue;
                }
                double delta = 0;
                int prev1, prev2;
                int next1, next2;
                point1 = solution.cycleIndices[c][i];
                prev1 = solution.cycleIndices[c][(i - 1 + cycleSizes[c]) % cycleSizes[c]];
                next1 = solution.cycleIndices[c][(i + 1) % cycleSizes[c]];
                point2 = solution.cycleIndices[c][j];
                prev2 = solution.cycleIndices[c][(j - 1 + cycleSizes[c]) % cycleSizes[c]];
                next2 = solution.cycleIndices[c][(j + 1) % cycleSizes[c]];

                delta -= distanceMatrix[point1][prev1];
                delta -= distanceMatrix[point2][next2];
                delta += distanceMatrix[point1][next2];
                delta += distanceMatrix[point2][prev1];

                if (delta < -epsilon)
                {
                    // This time we only do cuts (prev1 X point1) and (point2 X next2)
                    pq.push(EdgeMove({prev1, point1}, {point2, next2}, delta));
                }
            }
        }
    }
    int no = 0;
    while (!pq.empty())
    {
        Move best = pq.top();
        pq.pop();

        if (best.type == moveType::VERTEX)
        {
            // Get current positions of the points in the cycles
            auto [cycle1, pos1] = solution.getPointPosition(best.neighbourhood1.mid);
            auto [cycle2, pos2] = solution.getPointPosition(best.neighbourhood2.mid);
            point1 = best.neighbourhood1.mid;
            point2 = best.neighbourhood2.mid;

            // Ensure the points are in separate cycles
            if (cycle1 != cycle2)
            {
                // Get the prev and next indices for both points
                int prev1 = solution.cycleIndices[cycle1][(pos1 - 1 + cycleSizes[cycle1]) % cycleSizes[cycle1]];
                int next1 = solution.cycleIndices[cycle1][(pos1 + 1) % cycleSizes[cycle1]];
                int prev2 = solution.cycleIndices[cycle2][(pos2 - 1 + cycleSizes[cycle2]) % cycleSizes[cycle2]];
                int next2 = solution.cycleIndices[cycle2][(pos2 + 1) % cycleSizes[cycle2]];

                VertexNeighbourhood neighbourhood1(prev1, point1, next1);
                VertexNeighbourhood neighbourhood2(prev2, point2, next2);

                // Check legality by comparing prevs and nexts, including reversed cases
                if ((neighbourhood1 == best.neighbourhood1 && neighbourhood2 == best.neighbourhood2) || (neighbourhood2 == best.neighbourhood1 && neighbourhood1 == best.neighbourhood2))
                {
                    no++;
                    // std::cout << "Move vertices " << no << ": " << best.neighbourhood1.mid << " " << best.neighbourhood2.mid << " score " << best.delta << "   ---------   ";
                    // std::cout << "neighbours1 " << best.neighbourhood1.prev << " " << best.neighbourhood1.next << " neighbours2 " << best.neighbourhood1.prev << " " << best.neighbourhood2.next << std::endl;

                    //  Swap the points between the two cycles
                    std::swap(solution.cycleIndices[cycle1][pos1], solution.cycleIndices[cycle2][pos2]);
                    std::swap(solution.pointPositions[point1], solution.pointPositions[point2]);
                    // And new potential moves
                    // That is both swapped points and their neighbours
                    std::vector<int> pointsToUpdate = {point1, prev1, next1, point2, prev2, next2};
                    Cut cut1 = {prev1, point1};
                    Cut cut2 = {point1, next1};
                    Cut cut3 = {prev2, point2};
                    Cut cut4 = {point2, next2};
                    std::vector<Cut> cutsToUpdate = {cut1, cut2, cut3, cut4};
                    enqueueImprovingEdgeMoves(cutsToUpdate, solution, distanceMatrix, cycleSizes, pq, epsilon);
                    enqueueImprovingVertexMoves(pointsToUpdate, solution, distanceMatrix, cycleSizes, pq, epsilon);
                }
            }
        }
        else
        {
            // Get current positions of the points in the cycles
            auto [cycle1, pos1] = solution.getPointPosition(best.cut1.prev);
            auto [cycle2, pos2] = solution.getPointPosition(best.cut2.prev);
            int c = cycle1; // Both cycles are the same
            // Ensure position 1 is before position 2
            if (pos1 > pos2)
            {
                std::swap(pos1, pos2);
                std::swap(best.cut1, best.cut2);
            }
            point1 = solution.cycleIndices[c][pos1];
            point2 = solution.cycleIndices[c][pos2];
            //  Get the prev and next indices for both points
            int prev1 = solution.cycleIndices[c][(pos1 - 1 + cycleSizes[c]) % cycleSizes[c]];
            int next1 = solution.cycleIndices[c][(pos1 + 1) % cycleSizes[c]];
            int prev2 = solution.cycleIndices[c][(pos2 - 1 + cycleSizes[c]) % cycleSizes[c]];
            int next2 = solution.cycleIndices[c][(pos2 + 1) % cycleSizes[c]];

            Cut cut1;
            Cut cut2;
            int innerStart, innerEnd;
            std::vector<int> pointsToUpdate = {point1, point2};
            if (prev1 == best.cut1.next)
            {
                cut1 = Cut(prev1, point1);
                innerStart = pos1;
                pointsToUpdate.push_back(best.cut1.next);
            }
            else if (next1 == best.cut1.next)
            {
                cut1 = Cut(point1, next1);
                innerStart = (pos1 + 1) % cycleSizes[c];
                pointsToUpdate.push_back(best.cut1.next);
            }
            if (prev2 == best.cut2.next)
            {
                cut2 = Cut(prev2, point2);
                innerEnd = (pos2 - 1 + cycleSizes[c]) % cycleSizes[c];
                pointsToUpdate.push_back(best.cut2.next);
            }
            else if (next2 == best.cut2.next)
            {
                cut2 = Cut(point2, next2);
                innerEnd = pos2;
                pointsToUpdate.push_back(best.cut2.next);
            }
            if ((cut1.isSameDirectionAs(best.cut1) && cut2.isSameDirectionAs(best.cut2)) ||
                (cut1.isReversedOf(best.cut1) && cut2.isReversedOf(best.cut2)))
            {
                no++;
                std::reverse(solution.cycleIndices[c].begin() + innerStart, solution.cycleIndices[c].begin() + innerEnd + 1);
                solution.updatePointPositions();
                // Add new legal moves
                int outerIterator = (innerEnd + 1) % cycleSizes[c];
                while (outerIterator != innerStart)
                {
                    int innerIterator = innerStart % cycleSizes[c]; // skip the first cut (it has to update whole array)
                    while (innerIterator != innerEnd)               // skip the last cut for the same reason
                    {
                        // Skip reversing the entire cycle
                        if (outerIterator == innerIterator || (innerIterator + 1) % cycleSizes[c] == outerIterator || (outerIterator + 1) % cycleSizes[c] == innerIterator)
                        {
                            innerIterator++;
                            innerIterator = innerIterator % cycleSizes[c];
                            continue;
                        }

                        double delta = 0.0;
                        int next1, next2;
                        point1 = solution.cycleIndices[c][outerIterator];
                        next1 = solution.cycleIndices[c][(outerIterator + 1) % cycleSizes[c]];
                        point2 = solution.cycleIndices[c][innerIterator];
                        next2 = solution.cycleIndices[c][(innerIterator + 1) % cycleSizes[c]];

                        delta -= distanceMatrix[point1][next1];
                        delta -= distanceMatrix[point2][next2];
                        delta += distanceMatrix[next1][next2];
                        delta += distanceMatrix[point1][point2];
                        if (delta < -epsilon)
                        {
                            pq.push(EdgeMove({point1, next1}, {point2, next2}, delta));
                        }

                        innerIterator++;
                        innerIterator = innerIterator % cycleSizes[c];
                    }
                    outerIterator++;
                    outerIterator = outerIterator % cycleSizes[c];
                }
                Cut afterExchange1 = {best.cut1.prev, best.cut2.prev};
                Cut afterExchange2 = {best.cut1.next, best.cut2.next};
                std::vector<Cut> cutsToUpdate = {afterExchange1, afterExchange2};
                enqueueImprovingEdgeMoves(cutsToUpdate, solution, distanceMatrix, cycleSizes, pq, epsilon);
                enqueueImprovingVertexMoves(pointsToUpdate, solution, distanceMatrix, cycleSizes, pq, epsilon);
            }
        }
    }

    return solution;
}

std::vector<std::vector<int>> computeNearestneighbours(
    const std::vector<std::vector<double>> &distanceMatrix,
    int n = 10)
{
    size_t numberOfPoints = distanceMatrix.size();
    std::vector<std::vector<int>> nearestneighbours(numberOfPoints);

    for (size_t i = 0; i < numberOfPoints; ++i)
    {
        std::vector<std::pair<double, int>> distances;

        for (size_t j = 0; j < numberOfPoints; ++j)
        {
            if (i == j)
                continue;
            distances.emplace_back(distanceMatrix[i][j], j);
        }

        std::partial_sort(
            distances.begin(),
            distances.begin() + std::min(n, (int)distances.size()),
            distances.end());

        for (int k = 0; k < std::min(n, (int)distances.size()); ++k)
        {
            nearestneighbours[i].push_back(distances[k].second);
        }
    }

    return nearestneighbours;
}

Solution localSearchCandidates(Solution solution, const std::vector<std::vector<double>> &distanceMatrix)
{
    const double epsilon = 1e-6; // Avoid numerical errors
    size_t numberOfPoints = distanceMatrix.size();
    size_t pointsInCycle1 = solution.cycleIndices[0].size(); // Size of first cycle
    size_t pointsInCycle2 = solution.cycleIndices[1].size(); // Size of second cycle
    std::vector<size_t> cycleSizes = {pointsInCycle1, pointsInCycle2};
    int numberOfNearestPoints = 10;

    auto nearestneighbours = computeNearestneighbours(distanceMatrix, numberOfNearestPoints);

    bool improvement = true;
    while (improvement)
    {
        improvement = false;
        double minDelta = 0.0;
        int whichPoint;
        int insertedNeighbour;
        int whichSide;
        int removedPoint;
        int placedAfter;
        bool sameCycle = true;
        double oldCost, newCost;
        int prev, next;

        // For each cycle
        for (int c = 0; c < 2; c++)
        {
            // Check each point in cycle
            for (int i = 0; i < cycleSizes[c]; i++)
            {
                int currentPoint = solution.cycleIndices[c][i];

                // And its nearest neighbours
                for (int j = 0; j < numberOfNearestPoints; j++)
                {
                    int neighbour = nearestneighbours[currentPoint][j];
                    auto [neighbourCycle, neighbourPosition] = solution.getPointPosition(neighbour);
                    double delta = 0.0;

                    if (neighbourCycle == c)
                    {
                        if ((i + 1) % cycleSizes[c] == neighbourPosition || (neighbourPosition + 1) % cycleSizes[c] == i)
                        {
                            continue; // Skip if the point is already there
                        }
                        // Check both +1 and -1 position
                        std::vector<int> sides = {-1, 1};
                        for (const int side : sides)
                        {
                            int insertingPosition = (i + side + cycleSizes[c]) % cycleSizes[c];
                            int otherCurrentlyThere = solution.cycleIndices[c][insertingPosition];
                            newCost = distanceMatrix[currentPoint][neighbour];
                            oldCost = distanceMatrix[currentPoint][otherCurrentlyThere];
                            double insertNeighbourDelta = newCost - oldCost;
                            int indexConnectedToNeighbour = (neighbourPosition + side + cycleSizes[c]) % cycleSizes[c];
                            int connectedToNeighbour = solution.cycleIndices[c][indexConnectedToNeighbour];
                            newCost = distanceMatrix[otherCurrentlyThere][connectedToNeighbour];
                            oldCost = distanceMatrix[neighbour][connectedToNeighbour];
                            double secondEdgeDelta = newCost - oldCost;
                            delta = insertNeighbourDelta + secondEdgeDelta;
                            if (delta < minDelta - epsilon)
                            {
                                minDelta = delta;
                                improvement = true;
                                whichPoint = currentPoint;
                                insertedNeighbour = neighbour;
                                whichSide = side;
                                sameCycle = true;
                            }
                        }
                    }
                    else
                    {
                        // Calculate removal gain
                        prev = solution.cycleIndices[neighbourCycle][(neighbourPosition - 1 + cycleSizes[neighbourCycle]) % cycleSizes[neighbourCycle]];
                        next = solution.cycleIndices[neighbourCycle][(neighbourPosition + 1) % cycleSizes[neighbourCycle]];
                        oldCost = distanceMatrix[prev][neighbour] + distanceMatrix[neighbour][next];
                        newCost = distanceMatrix[prev][next];
                        double neighbourRemovalDelta = newCost - oldCost;

                        // Check both +1 and -1 position
                        std::vector<int> sides = {-1, 1};
                        for (const int side : sides)
                        {
                            int insertingPosition = (i + side + cycleSizes[c]) % cycleSizes[c];
                            int other = solution.cycleIndices[c][insertingPosition];
                            newCost = distanceMatrix[other][neighbour] + distanceMatrix[neighbour][currentPoint];
                            oldCost = distanceMatrix[currentPoint][other];
                            double insertNeighbourDelta = newCost - oldCost;
                            // Calculate the cost of removing any other point
                            for (int k = 0; k < cycleSizes[c]; k++)
                            {
                                int pointToRemove = solution.cycleIndices[c][k];
                                if (pointToRemove == currentPoint)
                                {
                                    continue; // It makes no sense to remove the point whose nearest neighbor we are trying to insert
                                }
                                prev = solution.cycleIndices[c][(k - 1 + cycleSizes[c]) % cycleSizes[c]];
                                if (prev == currentPoint && side == 1)
                                {
                                    prev = neighbour;
                                }
                                next = solution.cycleIndices[c][(k + 1) % cycleSizes[c]];
                                if (next == currentPoint && side == -1)
                                {
                                    next = neighbour;
                                }
                                oldCost = distanceMatrix[prev][pointToRemove] + distanceMatrix[pointToRemove][next];
                                newCost = distanceMatrix[prev][next];
                                double newRemovalDelta = newCost - oldCost;
                                // Calculate where we can insert this removed point
                                for (int m = 0; m < cycleSizes[neighbourCycle]; m++)
                                {
                                    prev = solution.cycleIndices[neighbourCycle][m];
                                    if (prev == neighbour)
                                    {
                                        prev = solution.cycleIndices[neighbourCycle][(m - 1 + cycleSizes[neighbourCycle]) % cycleSizes[neighbourCycle]];
                                    }
                                    next = solution.cycleIndices[neighbourCycle][(m + 1) % cycleSizes[neighbourCycle]];
                                    if (next == neighbour)
                                    {
                                        next = solution.cycleIndices[neighbourCycle][(m + 2) % cycleSizes[neighbourCycle]];
                                    }
                                    newCost = distanceMatrix[prev][pointToRemove] + distanceMatrix[pointToRemove][next];
                                    oldCost = distanceMatrix[prev][next];
                                    double insertionDelta = newCost - oldCost;
                                    delta = neighbourRemovalDelta + insertNeighbourDelta + newRemovalDelta + insertionDelta;
                                    if (delta < minDelta - epsilon)
                                    {
                                        minDelta = delta;
                                        improvement = true;
                                        whichPoint = currentPoint;
                                        insertedNeighbour = neighbour;
                                        whichSide = side;
                                        removedPoint = pointToRemove;
                                        placedAfter = prev;
                                        sameCycle = false;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if (improvement)
        {
            if (sameCycle == true)
            {
                // It's basically the same as edge exchange (rotation)
                auto [cycle1, smallerPos] = solution.getPointPosition(whichPoint);
                auto [cycle2, greaterPos] = solution.getPointPosition(insertedNeighbour);
                int c = cycle1; // Both cycles are the same
                // Ensure position 1 is before position 2
                if (smallerPos > greaterPos)
                {
                    std::swap(smallerPos, greaterPos);
                }
                int innerStart = smallerPos;
                int innerEnd = greaterPos;
                if (whichSide == 1)
                {
                    innerStart = smallerPos + 1;
                    innerEnd = greaterPos;
                }
                else // -1
                {
                    innerStart = smallerPos;
                    innerEnd = greaterPos - 1;
                }

                std::reverse(solution.cycleIndices[c].begin() + innerStart, solution.cycleIndices[c].begin() + innerEnd + 1);
                solution.updatePointPositions();
            }
            else
            {
                // 0. Remove insertedNeighbour from its current cycle
                auto [insertedCycle, insertedPos] = solution.getPointPosition(insertedNeighbour);
                solution.cycleIndices[insertedCycle].erase(solution.cycleIndices[insertedCycle].begin() + insertedPos);

                // 1. Find position of whichPoint
                auto [targetCycle, targetPos] = solution.getPointPosition(whichPoint);

                // 2. Insert insertedNeighbour before or after whichPoint depending on whichSide
                if (whichSide == -1)
                {
                    solution.cycleIndices[targetCycle].insert(solution.cycleIndices[targetCycle].begin() + targetPos, insertedNeighbour);
                }
                else // whichSide == 1
                {
                    solution.cycleIndices[targetCycle].insert(solution.cycleIndices[targetCycle].begin() + targetPos + 1, insertedNeighbour);
                }

                solution.updatePointPositions();

                // 4. Remove removedPoint from its current cycle
                auto [removedCycle, removedPos] = solution.getPointPosition(removedPoint);
                solution.cycleIndices[removedCycle].erase(solution.cycleIndices[removedCycle].begin() + removedPos);

                // 5. Insert removedPoint after placedAfter
                auto [placedAfterCycle, placedAfterPos] = solution.getPointPosition(placedAfter);
                solution.cycleIndices[placedAfterCycle].insert(solution.cycleIndices[placedAfterCycle].begin() + placedAfterPos + 1, removedPoint);

                solution.calculateScore(distanceMatrix);

                // 6. Update point positions again
                solution.updatePointPositions();
            }
        }
    }

    return solution;
}

void lab3(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix)
{
    int iterations = 1;

    // 0. Regret Cycle Weighted Heuristic
    double minHeuristicScore = std::numeric_limits<double>::max();
    double maxHeuristicScore = std::numeric_limits<double>::lowest();
    double totalHeuristicScore = 0.0;
    Solution bestHeuristicSolution;

    // 1. Random + Priority Queue + Memory + Steepest Descent + Edge Exchange
    double minRandomMemorySteepestEdgeScore = std::numeric_limits<double>::max();
    double maxRandomMemorySteepestEdgeScore = std::numeric_limits<double>::lowest();
    double totalRandomMemorySteepestEdgeScore = 0.0;
    Solution bestRandomMemorySteepestEdgeSolution;

    // 2. Random + Steepest Descent + Candidate Moves
    double minRandomSteepestCandidatesScore = std::numeric_limits<double>::max();
    double maxRandomSteepestCandidatesScore = std::numeric_limits<double>::lowest();
    double totalRandomSteepestCandidatesScore = 0.0;
    Solution bestRandomSteepestCandidatesSolution;

    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    for (int i = 0; i < iterations; ++i)
    {
        std::cout << "Progress: " << i << "/" << iterations << std::endl;
        // 0. Regret Cycle Weighted Heuristic
        Solution heuristic = regretCycleWeighted(distanceMatrix, 1.0, 1.0);
        heuristic.calculateScore(distanceMatrix);
        double score0 = heuristic.getScore();
        minHeuristicScore = std::min(minHeuristicScore, score0);
        maxHeuristicScore = std::max(maxHeuristicScore, score0);
        totalHeuristicScore += score0;
        if (score0 == minHeuristicScore)
            bestHeuristicSolution = heuristic;

        // Generate initial random and solution
        Solution initialRandom = randomCycle(distanceMatrix);

        // 1. Random + Priority Queue + Memory + Steepest Descent + Edge Exchange
        Solution solution1 = localSearchMemory(initialRandom, distanceMatrix);

        solution1.calculateScore(distanceMatrix);
        double score1 = solution1.getScore();
        minRandomMemorySteepestEdgeScore = std::min(minRandomMemorySteepestEdgeScore, score1);
        maxRandomMemorySteepestEdgeScore = std::max(maxRandomMemorySteepestEdgeScore, score1);
        totalRandomMemorySteepestEdgeScore += score1;
        if (score1 == minRandomMemorySteepestEdgeScore)
            bestRandomMemorySteepestEdgeSolution = solution1;
        // 2. Random + Steepest Descent + Candidate Moves
        Solution solution2 = localSearchCandidates(initialRandom, distanceMatrix);
        solution2.calculateScore(distanceMatrix);
        double score2 = solution2.getScore();
        minRandomSteepestCandidatesScore = std::min(minRandomSteepestCandidatesScore, score2);
        maxRandomSteepestCandidatesScore = std::max(maxRandomSteepestCandidatesScore, score2);
        totalRandomSteepestCandidatesScore += score2;
        if (score2 == minRandomSteepestCandidatesScore)
            bestRandomSteepestCandidatesSolution = solution2;
    }

    // Print results
    std::cout << "\nHeuristic (Regret Cycle Weighted):\n";
    std::cout << "Min: " << minHeuristicScore << " Max: " << maxHeuristicScore
              << " Avg: " << (totalHeuristicScore / iterations) << "\n";
    plotSolution(bestHeuristicSolution, points, distanceMatrix, "Heuristic (Regret Cycle Weighted)");

    std::cout << "\nRandom + Memory + Steepest Descent + Edge Exchange:\n";
    std::cout << "Min: " << minRandomMemorySteepestEdgeScore << " Max: " << maxRandomMemorySteepestEdgeScore
              << " Avg: " << (totalRandomMemorySteepestEdgeScore / iterations) << "\n";
    plotSolution(bestRandomMemorySteepestEdgeSolution, points, distanceMatrix, "Random + Memory + Steepest Descent + Edge Exchange");

    std::cout << "\nRandom + Steepest Descent + Candidate Moves:\n";
    std::cout << "Min: " << minRandomSteepestCandidatesScore << " Max: " << maxRandomSteepestCandidatesScore
              << " Avg: " << (totalRandomSteepestCandidatesScore / iterations) << "\n";
    plotSolution(bestRandomSteepestCandidatesSolution, points, distanceMatrix, "Random + Steepest Descent + Candidate Moves");
}

int main(int argc, char *argv[])
{
    std::string filepath = argv[1];
    points = loadPointsFromFile(filepath);
    if (points.empty())
    {
        std::cerr << "No points loaded from the file!" << std::endl;
        return -1;
    }

    // Create the adjacency list
    std::vector<std::vector<double>> distanceMatrix = createdistanceMatrix(points);

    lab3(points, distanceMatrix);

    return 0;
}