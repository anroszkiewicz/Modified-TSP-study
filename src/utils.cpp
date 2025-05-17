#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "matplotlibcpp.h"
#include "utils.h"
#include "solution.h"
#include "point.h"

namespace plt = matplotlibcpp;

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

void displayDistanceMatrix(const std::vector<std::vector<double>> &distanceMatrix)
{
    std::cout << "Distance Matrix (distances between points):\n";
    for (size_t i = 0; i < distanceMatrix.size(); ++i)
    {
        for (size_t j = 0; j < distanceMatrix[i].size(); ++j)
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

    // Print cycles and distances to console
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