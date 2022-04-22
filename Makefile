opt=-O2
compiler=g++-11
flags=-std=c++20 ${OPT} -Wall -Wextra -Wpedantic

all: KDTree

KDTree:
	${compiler} -c KDTree.cpp ${flags}

clean:
	rm *.o
