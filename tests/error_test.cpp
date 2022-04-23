#include "../KDTree.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>
#include <random>
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

std::vector<std::vector<double>> getListofGeneratedVectors(int length)
{
    std::vector<std::vector<double>> temp(length);
    for (size_t idx = 0; idx < length; idx++) {
        temp[idx] = generateVector();
    }
    return temp;
}

double sumSqrdErr(std::vector<double>& p1, std::vector<double>& p2)
{
    std::vector<double> diff(DIM);
    std::vector<double> square(DIM);
    std::transform(p1.begin(), p1.end(), p2.begin(), diff.begin(),
                   std::minus<double>());
    std::transform(diff.begin(), diff.end(), diff.begin(), square.begin(),
                   std::multiplies<double>());
    return std::accumulate(square.begin(), square.end(), 0.0);
}

int main()
{
    // seed
    srand(5);

    const std::vector<int> dataPointSizes = {50,  100, 200, 300, 400,  500,
                                             600, 700, 800, 900, 1000, 2000};
    constexpr int nIter = 50;

    std::cout << "Total number of iterations ran: " << nIter << std::endl;

    for (auto& sizes : dataPointSizes) {
        std::vector<std::vector<double>> points(sizes,
                                                std::vector<double>(DIM)),
            pointToRetrieve(sizes, std::vector<double>(DIM));
        long long kdTreeRetTotalTime = 0.0;
        long long bruteForceRetTotalTime = 0.0;

        int correct = 0;
        for (int i = 0; i < nIter; i++) {
            // generate test points to build a tree
            points = getListofGeneratedVectors(sizes);

            // genereate KD Tree
            KDTree tree;

            for (const auto& point : points) {
                tree.insert_point(point);
            }

            // generate retrieve test data points
            pointToRetrieve = getListofGeneratedVectors(sizes);

            for (auto& vals : pointToRetrieve) {
                double minSumSqdErr = std::numeric_limits<double>::max();
                std::vector<double> groundTruthVec(DIM);
                for (auto& gtvals :
                     points) // loop through all the points that built KDTRee
                {
                    double sumSqdErr = sumSqrdErr(gtvals, vals);
                    if (sumSqdErr < minSumSqdErr) {
                        minSumSqdErr = sumSqdErr;
                        groundTruthVec = gtvals;
                    }
                }
                std::vector<double> checkVec(DIM);
                checkVec = tree.nearest_point(vals);
                if (std::equal(groundTruthVec.begin(), groundTruthVec.end(),
                               checkVec.begin())) {
                    correct += 1;
                }
            }
        }
        std::cout << "\n\nAccuracy (tested with " << sizes
                  << " datasets per iter) = "
                  << ((correct * 100.0) / (sizes * nIter))
                  << " %.  Avg. Total Number Correct: "
                  << (int)(correct / nIter) << " / " << sizes;
    }
    return 0;
}
