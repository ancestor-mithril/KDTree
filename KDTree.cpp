/*
 * file: KDTree.cpp
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

#include "KDTree.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
constexpr auto level0 = 0;

double dist2(const point_t& a, const point_t& b)
{
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0,
                              std::plus<double>(),
                              [](auto x, auto y) { return (x - y) * (x - y); });
}

class comparer
{
  public:
    const std::size_t idx;
    comparer() = delete;
    comparer(std::size_t index) : idx{index} {};

    bool
    operator()(const std::vector<double>& a, const std::vector<double>& b) const
    {
        return a[idx] < b[idx];
    }
};

void sort_on_idx(pointVec::iterator begin, pointVec::iterator end,
                 std::size_t idx)
{
    std::nth_element(begin, begin + std::distance(begin, end) / 2, end,
                     comparer{idx});
}

KDNodePtr newKDNodePtr()
{
    return std::unique_ptr<KDNode>(nullptr);
}

KDNodePtr newKDNodePtr(const point_t& x, std::size_t idx)
{
    return std::make_unique<KDNode>(x, idx, newKDNodePtr(), newKDNodePtr());
}

} // namespace

KDNode::KDNode(const point_t& pt, std::size_t idx_, KDNodePtr&& left_,
               KDNodePtr&& right_)
    : x{pt}, index{idx_}, left{std::move(left_)}, right{std::move(right_)}
{
}

KDNode::KDNode(const pointIndex& pi, KDNodePtr&& left_, KDNodePtr&& right_)
    : x{pi.first}, index{pi.second}, left{std::move(left_)}, right{std::move(
                                                                 right_)}
{
}

double KDNode::coord(std::size_t idx) const
{
    return x[idx];
}

KDNode::operator std::size_t() const
{
    return index;
}

const point_t& KDNode::getPoint() const
{
    return x;
}

KDNodePtr
KDTree::makeTree(pointVec::iterator begin, pointVec::iterator end,
                 std::size_t beginIndex, std::size_t length, std::size_t level)
{
    if (begin == end) {
        return newKDNodePtr();
    }

    if (length > 1) {
        sort_on_idx(begin, end, level);
    }

    const auto middle = begin + (length / 2);
    const auto middleIndex = beginIndex + (length / 2);
    const auto l_len = length / 2;

    auto left = [&]() {
        const auto l_begin = begin;
        const auto l_end = middle;

        if (l_len > 0 and dim > 0) {
            return makeTree(l_begin, l_end, beginIndex, l_len,
                            (level + 1) % dim);
        }
        return newKDNodePtr();
    }();

    auto right = [&]() {
        const auto r_begin = middle + 1;
        const auto r_end = end;
        const auto r_len = length - l_len - 1;

        if (r_len > 0 && dim > 0) {
            return makeTree(r_begin, r_end, middleIndex + 1, r_len,
                            (level + 1) % dim);
        }
        return newKDNodePtr();
    }();

    return std::make_unique<KDNode>(*middle, middleIndex, std::move(left),
                                    std::move(right));
}

KDTree::KDTree(std::size_t dimensions)
    : dim{dimensions}, array_size{0}, root{newKDNodePtr()}
{
}

KDTree::KDTree(const pointVec& point_array, std::size_t dimensions)
    : dim{dimensions}, array_size{point_array.size()}
{
    rebuild(point_array);
}

void KDTree::rebuild(const pointVec& point_array)
{
    auto pointArray = point_array; // copy because we need to sort it
    root = makeTree(pointArray.begin(), pointArray.end(), 0, pointArray.size(),
                    level0);
    // begin, end, length, starting level
}

std::pair<const KDNode*, double>
KDTree::nearest_(const KDNode* branch, const point_t& pt, std::size_t level,
                 const KDNode* best, double best_dist, double rad) const
{
    if (not branch) {
        return {nullptr, best_dist};
    }

    const auto r2 = std::min(rad * rad, best_dist);
    const auto& branch_pt = branch->getPoint();
    const auto d = dist2(branch_pt, pt);

    if (d < r2) {
        best_dist = d;
        best = branch;
    }

    const auto dx = branch_pt[level] - pt[level];

    // select which branch makes sense to check
    const auto [section, other] = [=]() {
        if (dx > 0) {
            return std::make_pair(branch->left.get(), branch->right.get());
        }
        return std::make_pair(branch->right.get(), branch->left.get());
    }();

    const auto next_lv = (level + 1) % dim;
    // keep nearest neighbor from further down the tree
    const auto [further, furtherBest] =
        nearest_(section, pt, next_lv, best, best_dist, rad);
    if (further and furtherBest < r2) {
        best_dist = furtherBest;
        best = further;
    }

    const auto dx2 = dx * dx;
    // only check the other branch if it makes sense to do so
    if (dx2 < r2) {
        const auto [further, furtherBest] =
            nearest_(other, pt, next_lv, best, best_dist, rad);
        if (further and furtherBest < r2) {
            best_dist = furtherBest;
            best = further;
        }
    }

    if (best_dist < r2) {
        return {best, best_dist};
    }
    return {nullptr, best_dist};
};

// default caller
std::pair<const KDNode*, double> KDTree::nearest_(const point_t& pt) const
{
    if (not root) {
        throw std::logic_error("tree is empty");
    }

    const auto dist = dist2(root.get()->x, pt);

    return nearest_(root.get(),       // beginning of tree
                    pt,               // point we are querying
                    level0,           // start from level 0
                    root.get(),       // best is the root
                    dist,             // branch_dist
                    std::sqrt(dist)); // radius
    // radius must be sqrt(worst possible distance)
};

void KDTree::unsafeInsertPoint(const point_t& pt)
{
    auto current = root.get();
    auto level = level0;

    const auto current_size = array_size++;

    while (true) {
        if (pt[level] < current->x[level]) {
            if (not current->left) {
                current->left = newKDNodePtr(pt, current_size);
                return;
            } else {
                current = current->left.get();
            }
        } else {
            if (not current->right) {
                current->right = newKDNodePtr(pt, current_size);
                return;
            } else {
                current = current->right.get();
            }
        }

        level = (level + 1) % dim;
    }
}

void KDTree::insertPoint(const point_t& pt)
{
    if (not root) {
        // increasing array size
        root = newKDNodePtr(pt, array_size++);
        return;
    }
    unsafeInsertPoint(pt);
}

std::pair<std::size_t, double>
KDTree::nearestIndexAndValue(const point_t& pt) const
{
    const auto [node, dist] = nearest_(pt);
    // if root exists, then node is not null
    // if root does not exist, then error is thrown
    // we don't have to check for nullptr
    return {node->index, dist};
}

point_t KDTree::nearestPoint(const point_t& pt) const
{
    return nearest_(pt).first->getPoint();
}

std::size_t KDTree::nearestIndex(const point_t& pt) const
{
    return std::size_t(*nearest_(pt).first);
}

std::optional<std::size_t>
KDTree::nearestIndexWithinRange(const point_t& pt, double range) const
{
    if (not root) {
        throw std::logic_error("tree is empty");
    }

    const auto [node, dist] = nearest_(root.get(), // beginning of tree
                                       pt,         // point we are querying
                                       level0,     // start from level 0
                                       root.get(), // best is the root
                                       dist2(root.get()->x, pt), // branch_dist
                                       range);                   // radius
    if (node) {
        return node->index;
    }
    return std::nullopt;
}

pointIndex KDTree::nearestPointIndex(const point_t& pt) const
{
    const auto [Nearest, _] = nearest_(pt);
    return pointIndex(Nearest->getPoint(), std::size_t(*Nearest));
}

indexArr KDTree::neighborhood_(const KDNode* branch, const point_t& pt,
                               double rad, std::size_t level) const
{
    if (not branch) {
        return indexArr{};
    }

    const auto r2 = rad * rad;
    const auto d = dist2(branch->x, pt);

    indexArr nbh;

    if (d <= r2) {
        nbh.push_back(std::size_t(*branch));
    }

    const auto dx = branch->coord(level) - pt[level];
    const auto [section, other] = [=]() {
        if (dx > 0) {
            return std::make_pair(branch->left.get(), branch->right.get());
        }
        return std::make_pair(branch->right.get(), branch->left.get());
    }();

    const auto next_lv = (level + 1) % dim;
    const auto nbh_s = neighborhood_(section, pt, rad, next_lv);
    nbh.insert(nbh.end(), nbh_s.begin(), nbh_s.end());

    const auto dx2 = dx * dx;
    if (dx2 < r2) {
        const auto nbh_o = neighborhood_(other, pt, rad, next_lv);
        nbh.insert(nbh.end(), nbh_o.begin(), nbh_o.end());
    }

    return nbh;
}

indexArr KDTree::neighborhood(const point_t& pt, double rad) const
{
    return neighborhood_(root.get(), pt, rad, level0);
}

pointVec KDTree::neighborhoodPoints_(const KDNode* branch, const point_t& pt,
                                     double rad, std::size_t level) const
{
    if (not branch) {
        return pointVec{};
    }

    const auto r2 = rad * rad;
    const auto d = dist2(branch->x, pt);

    pointVec nbh;

    if (d <= r2) {
        nbh.push_back(branch->getPoint());
    }

    const auto dx = branch->coord(level) - pt[level];
    const auto [section, other] = [=]() {
        if (dx > 0) {
            return std::make_pair(branch->left.get(), branch->right.get());
        }
        return std::make_pair(branch->right.get(), branch->left.get());
    }();

    const auto dim = pt.size();
    const auto nbh_s = neighborhoodPoints_(section, pt, rad, (level + 1) % dim);
    nbh.insert(nbh.end(), nbh_s.begin(), nbh_s.end());

    const auto dx2 = dx * dx;
    if (dx2 < r2) {
        const auto nbh_o =
            neighborhoodPoints_(other, pt, rad, (level + 1) % dim);
        nbh.insert(nbh.end(), nbh_o.begin(), nbh_o.end());
    }

    return nbh;
}

// TODO: better implementation:
// have a list of points in the tree and use neighborhoodIndices to get the
// points
pointVec KDTree::neighborhoodPoints(const point_t& pt, double rad) const
{
    return neighborhoodPoints_(root.get(), pt, rad, level0);
}

indexArr KDTree::neighborhoodIndices(const point_t& pt, double rad) const
{
    return neighborhood_(root.get(), pt, rad, level0);
}

std::optional<std::size_t>
KDTree::firstNeighbor(const point_t& pt, double rad) const
{
    return firstNeighbor_(root.get(), pt, rad, level0);
}

std::optional<std::size_t>
KDTree::firstNeighbor_(const KDNode* branch, const point_t& pt, double rad,
                       std::size_t level) const
{
    if (not branch) {
        return std::nullopt;
    }

    const auto r2 = rad * rad;
    const auto d = dist2(branch->getPoint(), pt);

    if (d <= r2) {
        return branch->index;
    }

    const auto dx = branch->coord(level) - pt[level];

    const auto [section, other] = [=]() {
        if (dx > 0) {
            return std::make_pair(branch->left.get(), branch->right.get());
        }
        return std::make_pair(branch->right.get(), branch->left.get());
    }();

    const auto dim = pt.size();
    const auto index = firstNeighbor_(section, pt, rad, (level + 1) % dim);

    if (index) {
        return index;
    }

    if (dx * dx < r2) {
        return firstNeighbor_(other, pt, rad, (level + 1) % dim);
    }

    return std::nullopt;
}
