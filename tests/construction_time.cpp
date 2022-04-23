#include "../KDTree.hpp"
#include <algorithm>
#include <chrono>
#include <ctime>
#include <iostream>
#include <numeric>
#include <random>
#include <ratio>
#include <vector>

#define DIM 3

double getNum()
{
    return ((double)rand() / (RAND_MAX));
}

std::vector<double> generateVector()
{
    std::vector<double> temp(DIM);
    for (size_t idx = 0; idx < DIM; idx++) {
        temp[idx] = getNum();
    }
    return temp;
}

std::vector<std::vector<double>> getListofGeneratedVectors(size_t length)
{
    std::vector<std::vector<double>> temp(length);
    for (size_t idx = 0; idx < length; idx++) {
        temp[idx] = generateVector();
    }
    return temp;
}

int main()
{
    // seed
    srand(5);

    size_t npoints = 400000;
    std::cout << "constructing KDTree with " << npoints << " points."
              << std::endl;

    const auto points = getListofGeneratedVectors(npoints);
    const auto points2 = getListofGeneratedVectors(npoints);

    auto start = std::chrono::high_resolution_clock::now();
    KDTree tree(points);
    auto stop = std::chrono::high_resolution_clock::now();
    auto timespan = std::chrono::duration<double>(stop - start);
    std::cout << "it took " << timespan.count() << " seconds." << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (const auto& point : points2) {
        tree.insertPoint(point);
    }

    stop = std::chrono::high_resolution_clock::now();
    timespan = std::chrono::duration<double>(stop - start);
    std::cout << "Inserting took " << timespan.count() << " seconds."
              << std::endl;

    const auto points3 = getListofGeneratedVectors(npoints);

    const auto repeat = npoints;

    {
        KDTree tree(points);
        start = std::chrono::high_resolution_clock::now();

        for (const auto& point : points3) {
            auto nearest = tree.nearestIndex(point);
        }

        stop = std::chrono::high_resolution_clock::now();
        timespan = std::chrono::duration<double>(stop - start);
        std::cout << "Searching with nearest within balanced took "
                  << timespan.count() << " seconds." << std::endl;
    }

    {
        KDTree tree;
        for (const auto& point : points2) {
            tree.insertPoint(point);
        }

        start = std::chrono::high_resolution_clock::now();
        for (const auto& point : points3) {
            auto nearest = tree.nearestIndex(point);
        }

        stop = std::chrono::high_resolution_clock::now();
        timespan = std::chrono::duration<double>(stop - start);
        std::cout << "Searching with nearest within unbalanced took "
                  << timespan.count() << " seconds." << std::endl;
    }

    {
        KDTree tree(points);

        start = std::chrono::high_resolution_clock::now();
        const auto epsilon = 0.01;
        for (const auto& point : points3) {
            auto nearest = tree.neighborhood(point, epsilon);
        }
        stop = std::chrono::high_resolution_clock::now();
        timespan = std::chrono::duration<double>(stop - start);
        std::cout << "Searching with neighborhood within balanced 0.1 took "
                  << timespan.count() << " seconds." << std::endl;

        start = std::chrono::high_resolution_clock::now();
        for (const auto& point : points3) {
            auto nearest = tree.firstNeighbor(point, epsilon);
        }
        stop = std::chrono::high_resolution_clock::now();
        timespan = std::chrono::duration<double>(stop - start);
        std::cout << "First Neighbor within balanced " << epsilon << " took "
                  << timespan.count() << " seconds." << std::endl;
    }

    {
        KDTree tree;
        for (const auto& point : points2) {
            tree.insertPoint(point);
        }

        start = std::chrono::high_resolution_clock::now();
        const auto epsilon = 0.01;
        for (const auto& point : points3) {
            auto nearest = tree.neighborhood(point, epsilon);
        }
        stop = std::chrono::high_resolution_clock::now();
        timespan = std::chrono::duration<double>(stop - start);
        std::cout << "Searching with neighborhood within unbalanced " << epsilon
                  << " took " << timespan.count() << " seconds." << std::endl;

        start = std::chrono::high_resolution_clock::now();
        for (const auto& point : points3) {
            auto nearest = tree.firstNeighbor(point, epsilon);
        }
        stop = std::chrono::high_resolution_clock::now();
        timespan = std::chrono::duration<double>(stop - start);
        std::cout << "First Neighbor within unbalanced " << epsilon << " took "
                  << timespan.count() << " seconds." << std::endl;
    }

    {
        KDTree tree;
        for (const auto& point : points2) {
            tree.insertPoint(point);
        }

        start = std::chrono::high_resolution_clock::now();
        const auto epsilon = 0.001;
        for (const auto& point : points3) {
            auto nearest = tree.neighborhood(point, epsilon);
        }
        stop = std::chrono::high_resolution_clock::now();
        timespan = std::chrono::duration<double>(stop - start);
        std::cout << "Searching with neighborhood within unbalanced " << epsilon
                  << " took " << timespan.count() << " seconds." << std::endl;

        start = std::chrono::high_resolution_clock::now();
        for (const auto& point : points3) {
            auto nearest = tree.firstNeighbor(point, epsilon);
        }
        stop = std::chrono::high_resolution_clock::now();
        timespan = std::chrono::duration<double>(stop - start);
        std::cout << "First Neighbor within unbalanced " << epsilon << " took "
                  << timespan.count() << " seconds." << std::endl;
    }

    return 0;
}
