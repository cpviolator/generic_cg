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
./test
CG iter 1, rsq = 0.22321
CG iter 2, rsq = 0.0731392
CG iter 3, rsq = 0.0475098
CG iter 4, rsq = 0.0191602
CG iter 5, rsq = 0.0134318
CG iter 6, rsq = 0.00405207
CG iter 7, rsq = 0.00614416
CG iter 8, rsq = 0.00242465
CG iter 9, rsq = 0.000982737
CG iter 10, rsq = 0.00153331
CG iter 11, rsq = 0.000634916
CG iter 12, rsq = 0.000357081
CG iter 13, rsq = 0.000369652
CG iter 14, rsq = 0.000171095
CG iter 15, rsq = 4.01766e-05
CG iter 16, rsq = 2.32297e-05
CG iter 17, rsq = 1.90525e-05
CG iter 18, rsq = 1.22e-05
CG iter 19, rsq = 1.51786e-05
CG iter 20, rsq = 6.5849e-06
CG iter 21, rsq = 1.19161e-05
CG iter 22, rsq = 5.33034e-06
CG iter 23, rsq = 2.42755e-06
CG iter 24, rsq = 4.23812e-06
CG iter 25, rsq = 1.88424e-06
CG iter 26, rsq = 1.05526e-06
CG iter 27, rsq = 1.2778e-06
CG iter 28, rsq = 5.68012e-07
CG iter 29, rsq = 1.06401e-06
CG iter 30, rsq = 4.76929e-07
CG iter 31, rsq = 1.62974e-07
CG iter 32, rsq = 2.24951e-07
CG iter 33, rsq = 7.99395e-08
CG iter 34, rsq = 2.07752e-08
CG iter 35, rsq = 7.47581e-09
CG iter 36, rsq = 1.18944e-09
CG iter 37, rsq = 3.31833e-09
CG iter 38, rsq = 1.83495e-09
CG iter 39, rsq = 4.72716e-10
CG iter 40, rsq = 1.28103e-09
CG iter 41, rsq = 6.95287e-10
CG iter 42, rsq = 1.59038e-10
CG iter 43, rsq = 3.80447e-10
CG iter 44, rsq = 1.96258e-10
CG iter 45, rsq = 6.42723e-11
CG iter 46, rsq = 1.93356e-11
CG iter 47, rsq = 4.17772e-12
CG iter 48, rsq = 5.47976e-12
CG iter 49, rsq = 1.82868e-12
CG iter 50, rsq = 3.10877e-13
CG iter 51, rsq = 8.64289e-13
CG iter 52, rsq = 4.74599e-13
CG iter 53, rsq = 1.03973e-13
CG iter 54, rsq = 1.98465e-13
CG iter 55, rsq = 9.24602e-14
CG iter 56, rsq = 4.64662e-14
CG iter 57, rsq = 7.17611e-14
CG iter 58, rsq = 3.12874e-14
CG iter 59, rsq = 3.0569e-14
CG iter 60, rsq = 1.54045e-14
CG iter 61, rsq = 1.38146e-14
CG iter 62, rsq = 3.88254e-15
CG iter 63, rsq = 3.21085e-15
CG iter 64, rsq = 1.79058e-15
CG iter 65, rsq = 3.23747e-15
CG iter 66, rsq = 1.45667e-15
CG iter 67, rsq = 7.24703e-16
CG iter 68, rsq = 1.19697e-15
CG iter 69, rsq = 5.34587e-16
CG iter 70, rsq = 4.51721e-16
CG iter 71, rsq = 2.77323e-16
CG iter 72, rsq = 4.43518e-16
CG iter 73, rsq = 1.96571e-16
CG iter 74, rsq = 1.86031e-16
CG iter 75, rsq = 1.00107e-16
CG iter 76, rsq = 1.81349e-16
CG iter 77, rsq = 7.94729e-17
CG iter 78, rsq = 2.67478e-17
CG iter 79, rsq = 4.73908e-17
CG iter 80, rsq = 1.99159e-17
CG iter 81, rsq = 4.02455e-18
CG iter 82, rsq = 1.26142e-17
CG iter 83, rsq = 7.33319e-18
CG iter 84, rsq = 1.58941e-18
CG iter 85, rsq = 2.75562e-18
CG iter 86, rsq = 1.2567e-18
CG iter 87, rsq = 9.91884e-19
CG iter 88, rsq = 7.00708e-19
CG iter 89, rsq = 7.99653e-19
CG iter 90, rsq = 3.93433e-19
CG iter 91, rsq = 1.28058e-19
CG iter 92, rsq = 4.43932e-20
CG iter 93, rsq = 1.1473e-19
CG iter 94, rsq = 6.09744e-20
CG iter 95, rsq = 1.53772e-20
CG iter 96, rsq = 4.5342e-20
CG iter 97, rsq = 2.56199e-20
CG iter 98, rsq = 5.6015e-21
source norm = 1
sol norm = 0.025282
CG: Converged iter = 98, rsq = 5.6014954288578712e-21, truersq = 5.6014978461128073e-21
test routine source norm = 1
test routine residual norm = 7.48432e-11
test routine solution norm = 0.025282
```
The numerical values may differ due to the RNG.

The code has plenty of comments detailing the CG process

Happy hacking!
