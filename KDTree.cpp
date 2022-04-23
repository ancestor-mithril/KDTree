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
    const auto min = std::min(a.size(), b.size());
    const auto end = std::next(a.begin(), min);

    return std::inner_product(a.begin(), end, b.begin(), 0.0,
                              std::plus<double>(),
                              [](auto x, auto y) { return (x - y) * (x - y); });
}

class comparer
{
  public:
    const std::size_t idx;
    comparer() = delete;
    comparer(std::size_t index) : idx{index} {};

    bool operator()(const std::pair<std::vector<double>, std::size_t>& a,
                    const std::pair<std::vector<double>, std::size_t>& b)
    {
        return a.first[idx] < b.first[idx];
    }
};

void sort_on_idx(pointIndexArr::iterator begin, pointIndexArr::iterator end,
                 std::size_t idx)
{
    std::nth_element(begin, begin + std::distance(begin, end) / 2, end,
                     comparer{idx});
}

KDNodePtr newKDNodePtr()
{
    return std::unique_ptr<KDNode>(nullptr);
}

KDNodePtr newKDNodePtr(std::size_t idx)
{
    return std::make_unique<KDNode>(idx, newKDNodePtr(), newKDNodePtr());
}

} // namespace

KDNode::KDNode(std::size_t idx_, KDNodePtr&& left_,
               KDNodePtr&& right_)
    : index{idx_}, left{std::move(left_)}, right{std::move(right_)}
{
}

KDNode::operator std::size_t() const
{
    return index;
}

KDNodePtr
KDTree::makeTree(pointIndexArr::iterator begin, pointIndexArr::iterator end,
                 std::size_t length, std::size_t level)
{
    if (begin == end) {
        return newKDNodePtr(); // empty tree
    }

    const auto dim = begin->first.size();

    if (length > 1) {
        sort_on_idx(begin, end, level);
    }

    const auto middle = begin + (length / 2);
    const auto l_len = length / 2;

    auto left = [=, this]() {
        const auto l_begin = begin;
        const auto l_end = middle;

        if (l_len > 0 and dim > 0) {
            return makeTree(l_begin, l_end, l_len, (level + 1) % dim);
        }
        return newKDNodePtr();
    }();

    auto right = [=, this]() {
        const auto r_begin = middle + 1;
        const auto r_end = end;
        const auto r_len = length - l_len - 1;

        if (r_len > 0 && dim > 0) {
            return makeTree(r_begin, r_end, r_len, (level + 1) % dim);
        }
        return newKDNodePtr();
    }();

    return std::make_unique<KDNode>(middle->second, std::move(left), std::move(right));
}

KDTree::KDTree(const pointVec& point_array) : array_size(point_array.size())
{
    pointIndexArr arr;
    arr.reserve(point_array.size()); // allocating memory once
    std::transform(
        point_array.begin(), point_array.end(), std::back_inserter(arr),
        [index = std::size_t{}](
            const auto& point) mutable { // index is a member of the lambda
            // mutable to be able to change it from inside the lambda
            return std::make_pair(point, index++);
        });

    points = point_array;
    // TODO: This is not good, because points will actually be sorted. Or is it good? 
    root = KDTree::makeTree(arr.begin(), arr.end(), arr.size(), 0);
    // begin, end, length, starting level
}

void KDTree::reserve(std::size_t size)
{
    points.reserve(size);
}

std::pair<const KDNode*, double>
KDTree::nearest_(const KDNode* branch, const point_t& pt, std::size_t level,
                 const KDNode* best, double best_dist) const
{
    if (not branch) {
        return {nullptr, best_dist}; // basically, null
    }

    const auto& branch_pt = points[branch->index];

    const auto d = dist2(branch_pt, pt);
    // TODO: Measure if it makes sense to use this if to speed up same point retrieval
    // if (d == 0.0) {
    //     return {branch, d};
    // }

    auto best_l = best;
    auto best_dist_l = best_dist;

    if (d < best_dist) {
        best_dist_l = d;
        best_l = branch;
    }

    const auto dx = branch_pt[level] - pt[level];

    // select which branch makes sense to check
    const auto [section, other] = [=]() {
        if (dx > 0) {
            return std::make_pair(branch->left.get(), branch->right.get());
        }
        return std::make_pair(branch->right.get(), branch->left.get());
    }();

    const auto dim = branch_pt.size();
    const auto next_lv = (level + 1) % dim;

    // keep nearest neighbor from further down the tree
    const auto [further, furtherBest] =
        nearest_(section, pt, next_lv, best_l, best_dist_l);
    if (further) {
        // TODO: Measure if it makes sense to use this if to speed up same point retrieval
        // if (furtherBest == 0.0) {
        //     return {further, furtherBest};
        // }
        if (furtherBest < best_dist_l) {
            best_dist_l = furtherBest;
            best_l = further;
        }
    }

    // only check the other branch if it makes sense to do so
    if (dx * dx < best_dist_l) {
        const auto [further, furtherBest] =
            nearest_(other, pt, next_lv, best_l, best_dist_l);
        if (further) {
            if (furtherBest < best_dist_l) {
                best_dist_l = furtherBest;
                best_l = further;
            }
        }
    }

    return {best_l, best_dist_l};
};

// default caller
std::pair<const KDNode*, double> KDTree::nearest_(const point_t& pt) const
{
    if (not root) {
        throw std::logic_error("tree is empty");
    }

    return nearest_(root.get(),                // beginning of tree
                    pt,                        // point we are querying
                    level0,                    // start from level 0
                    root.get(),                // best is the root
                    dist2(points[root.get()->index], pt)); // branch_dist
};

double KDTree::coord(std::size_t index, std::size_t dimension) const
{
    return points[index][dimension];
}
double KDTree::coord(const KDNode* node, std::size_t dimension) const
{
    return points[node->index][dimension];
}

void KDTree::unsafeInsertPoint(const point_t& pt)
{
    points.push_back(pt);
    auto current = root.get();
    auto level = level0;

    const auto current_size = points.size() - 1;
    const auto dim = pt.size();

    while (true) {
        if (pt[level] < coord(current, level)) {
            if (not current->left) {
                current->left = newKDNodePtr(current_size);
                return;
            } else {
                current = current->left.get();
            }
        } else {
            if (not current->right) {
                current->right = newKDNodePtr(current_size);
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
        points.push_back(pt);
        root = newKDNodePtr(points.size() - 1);
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

const point_t& KDTree::point(const KDNode* node) const
{
    return points[node->index];
}

const point_t& KDTree::nearestPoint(const point_t& pt) const
{
    return point(nearest_(pt).first);
}

std::size_t KDTree::nearestIndex(const point_t& pt) const
{
    return std::size_t(*nearest_(pt).first);
}

pointIndex KDTree::nearestPointIndex(const point_t& pt) const
{
    const auto [nearest, _] = nearest_(pt);
    return pointIndex(point(nearest), std::size_t(*nearest));
}

indexArr KDTree::neighborhood_(const KDNode* branch, const point_t& pt,
                               double rad, std::size_t level) const
{
    if (not branch) {
        // check against empty branch ( nullptr or default constructed unique
        // pointer )
        return indexArr{};
    }

    const auto r2 = rad * rad;
    const auto d = dist2(point(branch), pt);

    indexArr nbh;

    if (d <= r2) {
        nbh.push_back(std::size_t(*branch));
    }

    const auto dx = coord(branch, level) - pt[level];
    const auto [section, other] = [=]() {
        if (dx > 0) {
            return std::make_pair(branch->left.get(), branch->right.get());
        }
        return std::make_pair(branch->right.get(), branch->left.get());
    }();

    const auto dim = pt.size();
    const auto nbh_s = neighborhood_(section, pt, rad, (level + 1) % dim);
    nbh.insert(nbh.end(), nbh_s.begin(), nbh_s.end());

    const auto dx2 = dx * dx;
    if (dx2 < r2) {
        const auto nbh_o = neighborhood_(other, pt, rad, (level + 1) % dim);
        nbh.insert(nbh.end(), nbh_o.begin(), nbh_o.end());
    }

    return nbh;
}

indexArr KDTree::neighborhood(const point_t& pt, double rad) const
{
    return neighborhood_(root.get(), pt, rad, level0);
}

indexArr KDTree::neighborhoodIndices(const point_t& pt, double rad) const
{
    return neighborhood_(root.get(), pt, rad, level0);
    ;
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
        // check against empty branch ( nullptr or default constructed unique
        // pointer )
        return std::nullopt;
    }

    const auto r2 = rad * rad;
    const auto d = dist2(point(branch), pt);

    if (d <= r2) {
        return branch->index;
    }

    const auto dx = coord(branch, level) - pt[level];

    const auto [section, other] = [=]() {
        if (dx > 0) {
            return std::make_pair(branch->left.get(), branch->right.get());
        }
        return std::make_pair(branch->right.get(), branch->left.get());
    }();

    const auto dim = pt.size();
    const auto index = firstNeighbor_(section, pt, rad, (level + 1) % dim);

    if (not index and dx * dx < r2) {
        return firstNeighbor_(other, pt, rad, (level + 1) % dim);
    }

    return std::nullopt;
}
