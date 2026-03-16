# Modernization Plan: Optimization Methods Solver

This plan outlines the strategy for transforming the legacy `main.cpp` optimization solver into a professional, modular, and fully verified C++ application.

## 1. Analysis & Baseline
- **Current State**: Monolithic `main.cpp` with mixed responsibilities (math, physics, optimizers, I/O).
- **Goal**: Establish a numerical baseline and maintain parity throughout refactoring.
- **Tasks**:
  - [ ] Capture the original output of `main.cpp` as a reference (`baseline_output.txt`).
  - [ ] Create a small regression script to compare new outputs with the baseline.

## 2. Infrastructure Setup
- **Goal**: Professional build and development environment.
- **Tasks**:
  - [ ] Migrate to **CMake (>= 3.14)**.
  - [ ] Integrate **Google Test** for unit testing.
  - [ ] Configure **.clang-format** (Google style).
  - [ ] Set up a modern directory structure (`src/`, `include/`, `tests/`).

## 3. Modularization & Refactoring
The code will be decomposed into logical components:

| Module | Responsibility | Planned Improvements |
| :--- | :--- | :--- |
| **`Geometry`** | Vector/Matrix math | Create a `Vector` class with operator overloads (`+`, `-`, `*`). Replace `vector<double>`. |
| **`Physics`** | Lennard-Jones Model | Extract `get_U_tot`, `get_U_at`, and `get_gradient` into a `LennardJonesSystem` class. |
| **`Optimizers`** | Solver Algorithms | Implement `PowellOptimizer` and `DynamicOptimizer` (Gradient Descent) as separate classes. |
| **`App`** | CLI & Orchestration | Modernize `main.cc` to handle parameters and output formatting using standard C++ streams. |

## 4. Documentation & UX
- **README.md**: Rewrite with LaTeX formulas and clear usage instructions.
- **Logging**: Add progress indicators for the optimization process.

## 5. Verification & Quality Assurance
- [ ] Implement unit tests for all mathematical primitives.
- [ ] Verify optimization convergence for standard cases.
- [ ] Replace `rand()` with `<random>` for reproducible benchmarks.
- [ ] Run `clang-tidy` to identify and fix legacy code smells.

## Final State
A clean, modularized codebase where each component is independently testable and the system is easy to extend with new optimization methods or physical models.
