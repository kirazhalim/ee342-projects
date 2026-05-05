# EE342 Projects

This repository contains two C++ numerical linear algebra projects developed for coursework. The projects include matrix/vector input files and sample outputs.

## Projects

- `project1`: matrix/vector operations and solver implementation with sample input files `A.txt` and `B.txt`.
- `project2`: matrix inverse/solver workflow with sample input files `A.txt` and `A_inv.txt`.

## What Is Excluded

Compiled executables, ZIP submissions, editor settings, assignment prompts, and docx reports are excluded.

## Build and Run

Compile with a standard C++ compiler:

```bash
g++ project1/abdulhalim_kiraz.cpp -o project1_solver
./project1_solver project1/A.txt project1/B.txt

g++ project2/abdulhalim_kiraz.cpp -o project2_solver
./project2_solver project2/A.txt
```

Adjust command-line arguments if you use the alternative implementation in `project1/abdulhalim_kiraz_ai.cpp`.
