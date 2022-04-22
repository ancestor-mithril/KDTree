BUILDDIR=./build

builddir:
	if [ -d "${BUILDDIR}" ]; then echo "Directory ${BUILDDIR} already exists"; else echo "Creating ${BUILDDIR}"; mkdir ${BUILDDIR}; fi;

clean:
	rm -rf ${BUILDDIR}

all: builddir
	cmake -S . -B ./build
	cmake --build ./build

FLAGS=-Wall -std=c++20 -Wextra -Wpedantic
COMPILER=g++-11
RELEASE=-O3 -DNDEBUG

compile: builddir
	cd ${BUILDDIR} \
	&& ${COMPILER} ${FLAGS} ${RELEASE} -c ../KDTree.cpp \
	&& ${COMPILER} ${FLAGS} ${RELEASE} -c ../tests/construction_time.cpp \
	&& ${COMPILER} ${FLAGS} ${RELEASE} KDTree.o construction_time.o -o construction_time.exe

run: compile
	./build/construction_time.exe
