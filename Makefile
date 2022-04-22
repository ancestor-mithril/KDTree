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
	&& ${COMPILER} ${FLAGS} ${RELEASE} -c ../tests/toy_test.cpp \
	&& ${COMPILER} ${FLAGS} ${RELEASE} -c ../tests/error_test.cpp \
	&& ${COMPILER} ${FLAGS} ${RELEASE} KDTree.o construction_time.o -o construction_time.exe \
	&& ${COMPILER} ${FLAGS} ${RELEASE} KDTree.o toy_test.o -o toy_test.exe \
	&& ${COMPILER} ${FLAGS} ${RELEASE} KDTree.o error_test.o -o error_test.exe

run: compile
	./build/error_test.exe
