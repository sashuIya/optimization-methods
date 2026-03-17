# Optimization Methods: Lennard-Jones Potential Minimizer

This project implements numerical optimization methods to find the minimum energy configuration of $n$ points in $m$-dimensional space, governed by the Lennard-Jones potential.

## Problem Overview

The potential energy between a pair of points $p_i$ and $p_j$ is defined as:
$$U_{ij} = \left(\frac{1}{r_{ij}}\right)^{12} - \left(\frac{1}{r_{ij}}\right)^6$$
where $r_{ij}$ is the Euclidean distance between points.

The goal is to minimize the average energy per point:
$$U_{at} = \frac{1}{n} \sum_{i < j} U_{ij}$$

## Optimization Methods

The solver utilizes two complementary methods to navigate the potential landscape and avoid local minima:

1.  **Dynamic Method (Gradient Descent with Momentum)**: Uses analytical gradients of the potential function.
2.  **Powell's Method**: A derivative-free optimization algorithm that performs line searches along a set of conjugate directions.
3.  **Golden Section Search**: Used for robust 1D line minimization.

## Quick Start

### Prerequisites
- CMake >= 3.14
- C++17 compliant compiler (GCC, Clang, MSVC)

### Build
```bash
mkdir build && cd build
cmake ..
make
```

### Usage
Run the solver for a range of points ($n$) and dimensions ($m$):
```bash
./task3 [n_start [n_end [m_start [m_end]]]]
# Example: solve only for n=5, m=3
./task3 5 5 3 3
```

### Testing
```bash
./unit_tests
```

## Project Structure
- `src/`: Core implementation and mathematical primitives.
- `tests/`: Unit tests using Google Test.
- `CMakeLists.txt`: Build configuration.

## License
Licensed under the [GNU GPL v3](LICENSE).
