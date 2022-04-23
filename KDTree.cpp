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

double KDNode::coord(size_t idx)
{
    return x[idx];
}
KDNode::operator bool() const
{
    return (!x.empty());
}
KDNode::operator point_t() const
{
    return x;
}
KDNode::operator size_t() const
{
    return index;
}
KDNode::operator pointIndex() const
{
    return pointIndex(x, index);
}

KDNodePtr NewKDNodePtr()
{
    return std::unique_ptr<KDNode>(nullptr);
}

double dist2(const point_t& a, const point_t& b)
{
    const auto min = std::min(a.size(), b.size());
    const auto end = std::next(a.begin(), min);

    return std::inner_product(a.begin(), end, b.begin(), 0.0,
                              std::plus<double>(),
                              [](auto x, auto y) { return (x - y) * (x - y); });
}

double dist2(const KDNode* a, const KDNode* b)
{
    return dist2(a->x, b->x);
}

double dist(const point_t& a, const point_t& b)
{
    return std::sqrt(dist2(a, b));
}

double dist(const KDNode* a, const KDNode* b)
{
    return std::sqrt(dist2(a, b));
}

comparer::comparer(size_t idx_) : idx{idx_} {};

bool comparer::compare_idx(const pointIndex& a, const pointIndex& b)
{
    return (a.first[idx] < b.first[idx]);
}

void sort_on_idx(pointIndexArr::iterator begin, pointIndexArr::iterator end,
                 size_t idx)
{
    const auto comp = comparer{idx};

    using std::placeholders::_1;
    using std::placeholders::_2;

    std::nth_element(begin, begin + std::distance(begin, end) / 2, end,
                     std::bind(&comparer::compare_idx, comp, _1, _2));
}

KDNodePtr
KDTree::make_tree(pointIndexArr::iterator begin, pointIndexArr::iterator end,
                  size_t length, size_t level)
{
    if (begin == end) {
        return NewKDNodePtr(); // empty tree
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
            return make_tree(l_begin, l_end, l_len, (level + 1) % dim);
        }
        return NewKDNodePtr();
    }();

    auto right = [=, this]() {
        const auto r_begin = middle + 1;
        const auto r_end = end;
        const auto r_len = length - l_len - 1;

        if (r_len > 0 && dim > 0) {
            return make_tree(r_begin, r_end, r_len, (level + 1) % dim);
        }
        return NewKDNodePtr();
    }();

    return std::make_unique<KDNode>(*middle, std::move(left), std::move(right));
}

KDTree::KDTree(const pointVec& point_array)
    : leaf{NewKDNodePtr()}, array_size(point_array.size())
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

    root = KDTree::make_tree(arr.begin(), arr.end(), arr.size(), 0);
    // begin, end, length, starting level
}

const KDNode*
KDTree::nearest_(const KDNode* branch, const point_t& pt, size_t level,
                 const KDNode* best, double best_dist)
{
    if (not branch) {
        return nullptr; // basically, null
    }

    const auto& branch_pt = branch->x;
    const auto dim = branch_pt.size();

    const auto d = dist2(branch_pt, pt);
    const auto dx = branch_pt[level] - pt[level];
    const auto dx2 = dx * dx;

    auto best_l = best;
    auto best_dist_l = best_dist;

    if (d < best_dist) {
        best_dist_l = d;
        best_l = branch;
    }

    const auto next_lv = (level + 1) % dim;
    // select which branch makes sense to check
    const auto [section, other] = [&]() -> std::pair<KDNode*, KDNode*> {
        if (dx > 0) {
            return {branch->left.get(), branch->right.get()};
        }
        return {branch->right.get(), branch->left.get()};
    }();

    // keep nearest neighbor from further down the tree
    const auto further = nearest_(section, pt, next_lv, best_l, best_dist_l);
    if (further) {
        const auto dl = dist2(further->x, pt);
        if (dl < best_dist_l) {
            best_dist_l = dl;
            best_l = further;
        }
    }

    // only check the other branch if it makes sense to do so
    if (dx2 < best_dist_l) {
        const auto further = nearest_(other, pt, next_lv, best_l, best_dist_l);
        if (further) {
            const auto dl = dist2(further->x, pt);
            if (dl < best_dist_l) {
                best_dist_l = dl;
                best_l = further;
            }
        }
    }

    return best_l;
};

// default caller
const KDNode* KDTree::nearest_(const point_t& pt)
{
    if (not root) {
        throw std::logic_error("tree is empty");
    }
    const auto level = std::size_t{};
    const auto branch_dist = dist2(root.get()->x, pt);

    return nearest_(root.get(),   // beginning of tree
                    pt,           // point we are querying
                    level,        // start from level 0
                    root.get(),   // best is the root
                    branch_dist); // best_dist = branch_dist
};

void KDTree::insert_point(const point_t& pt)
{
    auto current = root.get();
    auto level = level0;
    auto dim = pt.size();

    const auto current_size = array_size++;

    if (!root) {
        root = std::make_unique<KDNode>(pt, current_size, NewKDNodePtr(),
                                        NewKDNodePtr());
        return;
    }

    while (true) {
        if (pt[level] < current->x[level]) {
            if (not current->left) {
                current->left = std::make_unique<KDNode>(
                    pt, current_size, NewKDNodePtr(), NewKDNodePtr());
                return;
            } else {
                current = current->left.get();
            }
        } else {
            if (not current->right) {
                current->right = std::make_unique<KDNode>(
                    pt, current_size, NewKDNodePtr(), NewKDNodePtr());
                return;
            } else {
                current = current->right.get();
            }
        }

        level = (level + 1) % dim;
    }
}

point_t KDTree::nearest_point(const point_t& pt)
{
    const auto nearest = nearest_(pt);
    return point_t(*nearest);
};
size_t KDTree::nearest_index(const point_t& pt)
{
    return size_t(*nearest_(pt));
};

pointIndex KDTree::nearest_pointIndex(const point_t& pt)
{
    const auto Nearest = nearest_(pt);
    return pointIndex(point_t(*Nearest), size_t(*Nearest));
}

pointIndexArr KDTree::neighborhood_(const KDNode* branch, const point_t& pt,
                                    double rad, size_t level)
{
    if (not branch) {
        // check against empty branch ( nullptr or default constructed unique
        // pointer )
        return pointIndexArr{};
    }
    if (not bool(*branch)) {
        // branch has no point, means it is a leaf,
        // no points to add
        return pointIndexArr();
    }

    const auto dim = pt.size();

    const auto r2 = rad * rad;

    const auto d = dist2(branch->x, pt);
    const auto dx = point_t(*branch).at(level) - pt.at(level);
    const auto dx2 = dx * dx;

    pointIndexArr nbh;
    pointIndexArr nbh_s;
    pointIndexArr nbh_o;
    if (d <= r2) {
        nbh.push_back(pointIndex(*branch));
    }

    const auto [section, other] = [&]() {
        if (dx > 0) {
            return std::make_pair(branch->left.get(), branch->right.get());
        }
        return std::make_pair(branch->right.get(), branch->left.get());
    }();

    nbh_s = neighborhood_(section, pt, rad, (level + 1) % dim);
    nbh.insert(nbh.end(), nbh_s.begin(), nbh_s.end());
    if (dx2 < r2) {
        nbh_o = neighborhood_(other, pt, rad, (level + 1) % dim);
        nbh.insert(nbh.end(), nbh_o.begin(), nbh_o.end());
    }

    return nbh;
};

pointIndexArr KDTree::neighborhood(const point_t& pt, double rad)
{
    return neighborhood_(root.get(), pt, rad, level0);
}

pointVec KDTree::neighborhood_points(const point_t& pt, double rad)
{
    const auto nbh = neighborhood_(root.get(), pt, rad, level0);
    pointVec nbhp;
    nbhp.reserve(nbh.size()); // allocate memory, do not create default values
    std::transform(nbh.begin(), nbh.end(), std::back_inserter(nbhp),
                   [](const auto& x) { return x.first; }); // use const ref here
    return nbhp;
}

indexArr KDTree::neighborhood_indices(const point_t& pt, double rad)
{
    const auto nbh = neighborhood_(root.get(), pt, rad, level0);
    indexArr nbhi;
    nbhi.reserve(nbh.size()); // allocate memory, do not create default values
    std::transform(
        nbh.begin(), nbh.end(), std::back_inserter(nbhi),
        [](const auto& x) { return x.second; }); // use const ref here
    return nbhi;
}

std::optional<std::size_t>
KDTree::firstNeighbor(const point_t& pt, double rad) const
{
    return firstNeighbor_(root.get(), pt, rad, level0);
}

std::optional<std::size_t>
KDTree::firstNeighbor_(const KDNode* branch, const point_t& pt, double rad,
                       size_t level) const
{
    if (not branch) {
        // check against empty branch ( nullptr or default constructed unique
        // pointer )
        return std::nullopt;
    }
    if (not bool(*branch)) {
        return std::nullopt;
    }

    const auto dim = pt.size();

    const auto r2 = rad * rad;

    const auto d = dist2(point_t(*branch), pt);
    const auto dx = point_t(*branch).at(level) - pt.at(level);
    const auto dx2 = dx * dx;

    if (d <= r2) {
        return (*branch).index;
    }

    const auto [section, other] = [&]() -> std::pair<KDNode*, KDNode*> {
        if (dx > 0) {
            return {branch->left.get(), branch->right.get()};
        }
        return {branch->right.get(), branch->left.get()};
    }();

    const auto index = firstNeighbor_(section, pt, rad, (level + 1) % dim);

    if (not index and dx2 < r2) {
        return firstNeighbor_(other, pt, rad, (level + 1) % dim);
    }

    return std::nullopt;
}