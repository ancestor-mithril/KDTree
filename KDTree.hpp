#pragma once

/*
 * file: KDTree.hpp
 * author: G. Stoica and S. Amihaesei
 *
 * Based on the implemetation from https://github.com/crvs/KDTree by J.
 * Frederico Carvalho with elements from rosetta code from
 * https://github.com/crvs/KDTree. It is a reimplementation of the above
 * metioned versions, refactored and with methods added to support the
 * following:
 *  * insertion of single points
 *  * getting the first neighbor of a point within range
 */

#include <algorithm>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

using point_t = std::vector<double>;
using indexArr = std::vector<std::size_t>;
using pointIndex = typename std::pair<std::vector<double>, std::size_t>;

class KDNode
{
    using KDNodePtr = std::unique_ptr<KDNode>;

  public:
    const point_t x;
    const std::size_t index;
    KDNodePtr left;
    KDNodePtr right;

    // initializer
    KDNode() = delete;
    KDNode(const point_t&, std::size_t, KDNodePtr&&, KDNodePtr&&);
    KDNode(const pointIndex&, KDNodePtr&&, KDNodePtr&&);

    // getter
    double coord(std::size_t) const;
    const point_t& getPoint() const;

    // conversions
    explicit operator bool() const;
    explicit operator std::size_t() const;
};

using KDNodePtr = std::unique_ptr<KDNode>;

KDNodePtr nullKDNodePtr();

using pointIndexArr = typename std::vector<pointIndex>;

using pointVec = std::vector<point_t>;

class KDTree
{

  public:
    KDTree() = default;
    explicit KDTree(const pointVec& point_array);

    void insertPoint(const point_t& pt);
    void unsafeInsertPoint(const point_t& pt);
    const point_t& nearestPoint(const point_t& pt) const;
    std::size_t nearestIndex(const point_t& pt) const;
    std::pair<std::size_t, double>
    nearestIndexAndValue(const point_t& pt) const;
    pointIndex nearestPointIndex(const point_t& pt) const;

    indexArr neighborhood(const point_t& pt, double rad) const;

    // TODO: Move vector of points here
    // add neighbor points

    indexArr neighborhoodIndices(const point_t& pt, double rad) const;

    std::optional<std::size_t>
    firstNeighbor(const point_t& pt, double rad) const;

  private:
    KDNodePtr
    makeTree(pointIndexArr::iterator begin, pointIndexArr::iterator end,
             std::size_t length, std::size_t level);

    std::pair<const KDNode*, double>
    nearest_(const KDNode* branch, const point_t& pt, std::size_t level,
             const KDNode* best, double best_dist) const;

    // default caller
    std::pair<const KDNode*, double> nearest_(const point_t& pt) const;

    indexArr neighborhood_(const KDNode* branch, const point_t& pt, double rad,
                           std::size_t level) const;

    std::optional<std::size_t>
    firstNeighbor_(const KDNode* branch, const point_t& pt, double rad,
                   std::size_t level) const;

    std::size_t array_size = 0;
    KDNodePtr root;
};
