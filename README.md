# EE342 Numerical Linear Algebra Projects

This repository contains two C++ projects from my EE342 coursework. I worked on numerical linear algebra problems such as Gaussian elimination, matrix inversion, and eigenvalue-related calculations.

The repository has 2 project folders and 3 C++ source files. I also included sample input and output files so the results can be checked easily.

## What I Did

- Implemented Gaussian elimination with partial pivoting.
- Solved a linear system using matrix and vector input files.
- Implemented a matrix inverse workflow with a reusable matrix class.
- Compared outputs with the provided sample result files.
- Kept one AI-assisted version of the first solver as a separate file.

## Repository Structure

```text
project1/    Gaussian elimination solver
project2/    Matrix inverse and eigenvalue-related workflow
```

## How to Run
```
g++ project1/gaussian_elimination_solver.cpp -o project1_solver
./project1_solver project1/A.txt project1/B.txt

g++ project2/matrix_inverse_solver.cpp -o project2_solver
./project2_solver project2/A.txt
```
