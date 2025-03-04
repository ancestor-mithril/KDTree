# KDTree

Simple C++ KD-Tree implementation with minimal functionality.
Adapted from [original](https://github.com/crvs/KDTree) and improved.
- points are given as STL vectors (and inserted in their own STL vector) so supports n-dimensional points for any n
- makes full trees, (i.e. does not cut-off the branching at some arbitrary level) giving the nearest neighbor query have (strong) logarithmic complexity.
- builds the tree in one go ~~(does not support adding nodes, the tree is built from a list of points and cannot be altered afterwards)~~
- supports adding nodes and ~~rebalancing~~
- points are assumed to be STL vectors
- it provides the following queries:
	- nearest neighbor
	- neighbors within a given distance
	- first neighbor within a given distance
	- ~~nearest neighbor within a given distance~~

## License and copyright

© J. Frederico Carvalho (2018 - 2021)
Licensed under the [BSD3 License](ORIGINAL_LICENSE)

© G. Stoica and S. Amihaiesei (2022 - )
Licensed under the [GPL3 License](LICENSE)
