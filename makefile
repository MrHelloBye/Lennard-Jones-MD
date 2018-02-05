main.o: main.cpp
	clang++ main.cpp -c -std=c++14

Particle.o: Particle.cpp
	clang++ Particle.cpp -c -std=c++14

vose.o: vose.cpp
	clang++ vose.cpp -c std=c++14

ljmd: main.o Particle.o vose.o
	clang++ *.o -o ljmd
