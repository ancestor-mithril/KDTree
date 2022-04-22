opt=-O2
compiler=g++-11
flags=-std=c++20 ${OPT} -Wall -Wextra -Wpedantic

all: toy_test error_test_2 construction_time

run:
	./app1.exe
	./app2.exe
	./app3.exe

toy_test:
	${compiler} -c toy_test.cpp ${flags}
	${compiler} ${flags} ../KDTree.o toy_test.o -o app1.exe

error_test_2:
	${compiler} -c error_test_2.cpp ${flags}
	${compiler} ${flags} ../KDTree.o error_test_2.o -o app2.exe

construction_time:
	${compiler} -c construction_time.cpp ${flags}
	${compiler} ${flags} ../KDTree.o construction_time.o -o app3.exe

clean:
	rm *.o
