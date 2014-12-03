#g++ -fPIC -c lib.cpp
g++ -std=c++11 -fopenmp -O2 lib.o main.cpp Potential.cpp System.cpp AlphaHarmonicOscillator.cpp BasisFunctions.cpp Jastrow.cpp GaussianDeviate.cpp -larmadillo -llapack -lblas
