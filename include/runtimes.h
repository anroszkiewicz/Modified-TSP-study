#ifndef RUNTIMES_H
#define RUNTIMES_H

#include <vector>
#include "point.h"

void lab1(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix);

void lab2(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix);

void lab3(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix);

void lab4(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix);

void lab5(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix);

void lab5_heuristic(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix);

void lab5_evolutionary(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix);

void lab6(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix);

void convex_test_runtime(std::vector<Point> &points, std::vector<std::vector<double>> &distanceMatrix);

#endif