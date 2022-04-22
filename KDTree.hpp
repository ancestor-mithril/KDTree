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
#include <vector>

using point_t = std::vector< double >;
using indexArr = std::vector< size_t >;
using pointIndex = typename std::pair< std::vector< double >, size_t >;

class KDNode {
   public:
    // why not unique pointers 
    using KDNodePtr = std::shared_ptr< KDNode >;
    // const??
    point_t x;
    size_t index;
    KDNodePtr left;
    KDNodePtr right;

    // initializer
    KDNode() = default;
    KDNode(const point_t &, size_t, const KDNodePtr &,
           const KDNodePtr &);
    KDNode(const pointIndex &, const KDNodePtr &, const KDNodePtr &);

    // getter
    double coord(size_t);

    // conversions
    explicit operator bool();
    explicit operator point_t();
    explicit operator size_t();
    explicit operator pointIndex();
};

using KDNodePtr = std::shared_ptr< KDNode >;

KDNodePtr NewKDNodePtr();

// square euclidean distance
inline double dist2(const point_t &, const point_t &);
inline double dist2(const KDNodePtr &, const KDNodePtr &);

// euclidean distance
inline double dist(const point_t &, const point_t &);
inline double dist(const KDNodePtr &, const KDNodePtr &);

// Need for sorting
class comparer {
   public:
    size_t idx;
    explicit comparer(size_t idx_);
    inline bool compare_idx(
        const std::pair< std::vector< double >, size_t > &,
        const std::pair< std::vector< double >, size_t > &
    );
};

using pointIndexArr = typename std::vector< pointIndex >;

inline void sort_on_idx(pointIndexArr::iterator,
                        pointIndexArr::iterator,
                        size_t idx);

using pointVec = std::vector< point_t >;

class KDTree {
    KDNodePtr root;
    KDNodePtr leaf;

    KDNodePtr make_tree(pointIndexArr::iterator begin,
                        pointIndexArr::iterator end,
                        size_t length,
                        size_t level
    );

   public:
    KDTree() = default;
    explicit KDTree(const pointVec &point_array);

   private:
    KDNodePtr nearest_(
        const KDNodePtr &branch,
        const point_t &pt,
        size_t level,
        const KDNodePtr &best,
        double best_dist
    );

    // default caller
    KDNodePtr nearest_(const point_t &pt);

   public:
    point_t nearest_point(const point_t &pt);
    size_t nearest_index(const point_t &pt);
    pointIndex nearest_pointIndex(const point_t &pt);

   private:
    pointIndexArr neighborhood_(
        const KDNodePtr &branch,
        const point_t &pt,
        double rad,
        size_t level
    );

   public:
    pointIndexArr neighborhood(
        const point_t &pt,
        double rad);

    pointVec neighborhood_points(
        const point_t &pt,
        double rad);

    indexArr neighborhood_indices(
        const point_t &pt,
        double rad);
};
