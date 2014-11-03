#g++ -fPIC -c lib.cpp
g++ -std=c++11 -fopenmp -O3 lib.o main.cpp Potential.cpp System.cpp AlphaHarmonicOscillator.cpp BasisFunctions.cpp Jastrow.cpp -larmadillo -llapack -lblas
