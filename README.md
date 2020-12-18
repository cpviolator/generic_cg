# generic_cg

A generic CG linear operator

## Build

In case you didn't already, please ensure that you clone this repository using
the --recursive flag:
```
git clone --recursive https://github.com/cpviolator/generic_cg.git
```

Usual CMake build procedure. Create a build directory and invoke CMake via

```
mkdir build
cd build
cmake <path/to/generic_cg>
make
```
or, if you prefer the GUI
```
mkdir build
cd build
ccmake <path/to/generic_cg>
make
```

## Running

The library will create a test execuatble in `build/test`. You should get something
similar to this:
```
~/sambanova/build/test$ ./test 
CG iter 1, rsq = 0.0778255
CG iter 2, rsq = 0.0119006
CG iter 3, rsq = 0.00262384
CG iter 4, rsq = 0.000575028
CG iter 5, rsq = 0.000141541
CG iter 6, rsq = 3.16102e-05
CG iter 7, rsq = 7.8289e-06
CG iter 8, rsq = 1.65913e-06
CG iter 9, rsq = 3.85137e-07
CG iter 10, rsq = 7.32421e-08
CG iter 11, rsq = 1.44106e-08
CG iter 12, rsq = 2.36176e-09
CG iter 13, rsq = 3.89963e-10
CG iter 14, rsq = 5.5115e-11
CG iter 15, rsq = 7.81776e-12
CG iter 16, rsq = 9.66721e-13
CG iter 17, rsq = 1.19105e-13
CG iter 18, rsq = 1.31638e-14
CG iter 19, rsq = 1.50224e-15
CG iter 20, rsq = 1.51587e-16
CG iter 21, rsq = 1.47745e-17
CG iter 22, rsq = 1.41905e-18
CG iter 23, rsq = 1.7268e-19
CG iter 24, rsq = 1.51408e-20
CG iter 25, rsq = 4.12606e-22
source norm = 1
sol norm = 0.035301
CG: Converged iter = 25, rsq = 4.1260594705864369e-22, truersq = 4.1261018512255343e-22
test routine source norm = 1
test routine residual norm = 2.03128e-11
test routine solution norm = 0.035301
```
The numerical values may differ due to the RNG.

The code has plenty of comments detailing the CG process

Happy hacking!
