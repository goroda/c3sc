# Compressed Continuous Computation for Stochastic Optimal Control (C3SC)
Perform stochastic optimal control using the function-train decomposition from the C3 package

## Prerequisites
  * [Compressed Continuous Computation (C3)](https://github.com/goroda/Compressed-Continuous-Computation) 
  * [Dynamical Systems in C (cdyn)](https://github.com/goroda/cdyn)

This library is used for these papers

1. Alex Gorodetsky, Sertac Karaman, Youssef Marzouk: [High-Dimensional Stochastic Optimal Control using Continuous Tensor Decompositions](https://alexgorodetsky.com/wp-content/uploads/2018/02/1611.04706.pdf). In: International Journal of Robotics Research, Accepted 2018.
2. Ezra Tal, Alex Gorodetsky, Sertac Karaman: Continuous Tensor Train-Based Dynamic Programming for High-Dimensional Zero-Sum Differential Games. In: American Control Conference (ACC), Milwaukee, WI, Accepted 2018.

## Installation Instructions

We will install the prerequisite packages (C3 and CDYN) into a generic directory denoted by `<c3sc-prereq>`. Replace this flag with what is appropriate for your system. Then we will tell C3SC where to find the directory. We will install C3SC into the directory denoted by `<c3sc-installed>.`

### Install C3 

```
git clone https://github.com/goroda/Compressed-Continuous-Computation.git c3
cd c3
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=<c3sc-prereq> ..
make
make install
```

### Install CDYN

```shell
git clone https://github.com/goroda/cdyn.git cdyn
cd cdyn
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=<c3sc-prereq> ..
make
make install
```

### Install C3SC

```shell
git clone https://github.com/goroda/Compressed-Continuous-Computation.git c3
cd c3
mkdir build
cd build
cmake -DC3_INCLUDE_DIR=<c3sc-prereq>/include -DCDYN_INCLUDE_DIR=<c3sc-prereq>/include -DC3_LIB_PATH=<c3sc-prereq>/lib -DCDYN_LIB_PATH=<c3sc-prereq>/lib -DCMAKE_INSTALL_PREFIX=<c3sc-installed> ..
make
make install
```


<!-- http://www.alexgorodetsky.com/c3/html/ -->

Author: [Alex A. Gorodetsky](https://www.alexgorodetsky.com)  
Contact: [goroda@umich.edu](mailto:goroda@umich.edu)  
Copyright (c) 2015-2016, Massachusetts Institute of Technology  
Copyright (c) 2018, University of Michigan  
License: BSD
