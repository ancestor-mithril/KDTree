#pragma once

/*
 * file: KDTree.hpp
 * author: J. Frederico Carvalho
 *
 * This is an adaptation of the KD-tree implementation in rosetta code
 *  https://rosettacode.org/wiki/K-d_tree
 * It is a reimplementation of the C code using C++.
 * It also includes a few more queries than the original
 *
 */

#include <algorithm>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

using point_t = std::vector<double>;
using indexArr = std::vector<size_t>;
using pointIndex = typename std::pair<std::vector<double>, size_t>;

class KDNode
{
  public:
    // why not unique pointers
    using KDNodePtr = std::unique_ptr<KDNode>;
    // const??
    point_t x;
    size_t index;
    KDNodePtr left;
    KDNodePtr right;

    // initializer
    KDNode() = delete;
    KDNode(const point_t&, size_t, KDNodePtr&&, KDNodePtr&&);
    KDNode(const pointIndex&, KDNodePtr&&, KDNodePtr&&);

    // getter
    double coord(size_t);

    // conversions
    explicit operator bool() const;
    explicit operator point_t() const;
    explicit operator size_t() const;
    explicit operator pointIndex() const;
};

using KDNodePtr = std::unique_ptr<KDNode>;

KDNodePtr NewKDNodePtr();

// square euclidean distance
inline double dist2(const point_t&, const point_t&);
inline double dist2(const KDNode*, const KDNode*);

// euclidean distance
inline double dist(const point_t&, const point_t&);
inline double dist(const KDNode*, const KDNode*);

// Need for sorting
class comparer
{
  public:
    size_t idx;
    explicit comparer(size_t idx_);
    inline bool compare_idx(const std::pair<std::vector<double>, size_t>&,
                            const std::pair<std::vector<double>, size_t>&);
};

using pointIndexArr = typename std::vector<pointIndex>;

inline void
sort_on_idx(pointIndexArr::iterator, pointIndexArr::iterator, size_t idx);

using pointVec = std::vector<point_t>;

class KDTree
{
    size_t array_size = 0;
    KDNodePtr root;
    KDNodePtr leaf;

    KDNodePtr
    make_tree(pointIndexArr::iterator begin, pointIndexArr::iterator end,
              size_t length, size_t level);

  public:
    KDTree() = default;
    explicit KDTree(const pointVec& point_array);

  private:
    const KDNode* nearest_(const KDNode* branch, const point_t& pt,
                           size_t level, const KDNode* best, double best_dist);

    // default caller
    const KDNode* nearest_(const point_t& pt);

  public:
    void insert_point(const point_t& pt);
    point_t nearest_point(const point_t& pt);
    size_t nearest_index(const point_t& pt);
    pointIndex nearest_pointIndex(const point_t& pt);

  private:
    pointIndexArr neighborhood_(const KDNode* branch, const point_t& pt,
                                double rad, size_t level);
    std::optional<std::size_t>
    firstNeighbor_(const KDNode* branch, const point_t& pt, double rad,
                   size_t level) const;

  public:
    pointIndexArr neighborhood(const point_t& pt, double rad);

    pointVec neighborhood_points(const point_t& pt, double rad);

    indexArr neighborhood_indices(const point_t& pt, double rad);

    std::optional<std::size_t>
    firstNeighbor(const point_t& pt, double rad) const;
};
