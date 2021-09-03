OPT= -pthread -std=c++14 -Wall  -O3 -ffast-math -march=native  -mfpmath=sse  
a.out: main.o jordan.o time.o
	g++ $(OPT) jordan.o main.o time.o
main.o: main.cpp jordan.h
	g++ $(OPT) -c main.cpp
jordan.o: jordan.cpp jordan.h
	g++ $(OPT) -c jordan.cpp
time.o: time.cpp jordan.h
	g++ $(OPT) -c time.cpp
