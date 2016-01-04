all: main

main: main.cpp
	mpicxx --std=c++11 main.cpp grid.cpp matrix.cpp fox.cpp -o main

clean:
	rm ./main

run: main
	mpirun -n 4 ./main input.txt
