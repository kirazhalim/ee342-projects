/**
 * Gaussian Elimination with Partial Pivoting
 * 
 * This program implements Gaussian elimination with partial pivoting and backward
 * substitution to solve a linear system Ax = b, where A is an n×n square matrix.
 * 
 * The program reads matrix A and vector b from input files, solves the system, and
 * outputs the solution vector x to a text file.
 */

 #include <iostream>
 #include <fstream>
 #include <cmath>
 #include <cstring>
 #include <limits>
 
 // Machine epsilon (precision) for floating-point comparison
 const double EPSILON = std::numeric_limits<double>::epsilon() * 100;
 
 /**
  * Count the number of rows in a file
  * @param filename The file to read from
  * @return The number of rows in the file
  */
 int countRows(const char* filename) {
     std::ifstream file(filename);
     if (!file.is_open()) {
         std::cerr << "Error: Could not open file " << filename << std::endl;
         exit(1);
     }
     
     int rows = 0;
     char buffer[1024];
     
     while (file.getline(buffer, 1024)) {
         rows++;
     }
     
     file.close();
     return rows;
 }
 
 /**
  * Count the number of columns in the first row of a file
  * @param filename The file to read from
  * @return The number of columns in the first row
  */
 int countCols(const char* filename) {
     std::ifstream file(filename);
     if (!file.is_open()) {
         std::cerr << "Error: Could not open file " << filename << std::endl;
         exit(1);
     }
     
     int cols = 0;
     char buffer[1024];
     
     if (file.getline(buffer, 1024)) {
         char* token = strtok(buffer, " \t");
         while (token != NULL) {
             cols++;
             token = strtok(NULL, " \t");
         }
     }
     
     file.close();
     return cols;
 }
 
 /**
  * Reads matrix data from a file into a dynamically allocated array
  * @param filename The file to read from
  * @param rows Reference to store the number of rows
  * @param cols Reference to store the number of columns
  * @return A dynamically allocated 2D array containing the matrix data
  */
 double** readMatrixFromFile(const char* filename, int& rows, int& cols) {
     // Count rows and columns
     rows = countRows(filename);
     cols = countCols(filename);
     
     // Allocate memory for the matrix
     double** matrix = new double*[rows];
     for (int i = 0; i < rows; i++) {
         matrix[i] = new double[cols];
     }
     
     // Read the data from the file
     std::ifstream file(filename);
     if (!file.is_open()) {
         std::cerr << "Error: Could not open file " << filename << std::endl;
         exit(1);
     }
     
     for (int i = 0; i < rows; i++) {
         for (int j = 0; j < cols; j++) {
             if (!(file >> matrix[i][j])) {
                 std::cerr << "Error: Invalid data format in file " << filename << std::endl;
                 exit(1);
             }
         }
     }
     
     file.close();
     return matrix;
 }
 
 /**
  * Reads vector data from a file into a dynamically allocated array
  * @param filename The file to read from
  * @param size Reference to store the size of the vector
  * @return A dynamically allocated array containing the vector data
  */
 double* readVectorFromFile(const char* filename, int& size) {
     // Count rows (size of the vector)
     size = countRows(filename);
     
     // Allocate memory for the vector
     double* vector = new double[size];
     
     // Read the data from the file
     std::ifstream file(filename);
     if (!file.is_open()) {
         std::cerr << "Error: Could not open file " << filename << std::endl;
         exit(1);
     }
     
     for (int i = 0; i < size; i++) {
         if (!(file >> vector[i])) {
             std::cerr << "Error: Invalid data format in file " << filename << std::endl;
             exit(1);
         }
     }
     
     file.close();
     return vector;
 }
 
 /**
  * Free the memory allocated for a 2D matrix
  * @param matrix The matrix to free
  * @param rows The number of rows in the matrix
  */
 void freeMatrix(double** matrix, int rows) {
     for (int i = 0; i < rows; i++) {
         delete[] matrix[i];
     }
     delete[] matrix;
 }
 
 /**
  * Calculates the condition number of a 2×2 matrix using norm-1
  * @param A The 2×2 matrix
  * @return The condition number under norm-1
  */
 double conditionNumber1(double** A) {
     // Calculate determinant
     double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
     
     if (std::abs(det) < EPSILON) {
         return std::numeric_limits<double>::infinity(); // Singular matrix
     }
     
     // Calculate inverse
     double inv[2][2];
     inv[0][0] = A[1][1] / det;
     inv[0][1] = -A[0][1] / det;
     inv[1][0] = -A[1][0] / det;
     inv[1][1] = A[0][0] / det;
     
     // Calculate norm-1 for A (maximum absolute column sum)
     double normA = std::max(std::abs(A[0][0]) + std::abs(A[1][0]), 
                            std::abs(A[0][1]) + std::abs(A[1][1]));
     
     // Calculate norm-1 for A inverse
     double normAInv = std::max(std::abs(inv[0][0]) + std::abs(inv[1][0]), 
                               std::abs(inv[0][1]) + std::abs(inv[1][1]));
     
     return normA * normAInv;
 }
 
 /**
  * Calculates the condition number of a 2×2 matrix using norm-infinity
  * @param A The 2×2 matrix
  * @return The condition number under norm-infinity
  */
 double conditionNumberInf(double** A) {
     // Calculate determinant
     double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
     
     if (std::abs(det) < EPSILON) {
         return std::numeric_limits<double>::infinity(); // Singular matrix
     }
     
     // Calculate inverse
     double inv[2][2];
     inv[0][0] = A[1][1] / det;
     inv[0][1] = -A[0][1] / det;
     inv[1][0] = -A[1][0] / det;
     inv[1][1] = A[0][0] / det;
     
     // Calculate norm-infinity for A (maximum absolute row sum)
     double normA = std::max(std::abs(A[0][0]) + std::abs(A[0][1]), 
                            std::abs(A[1][0]) + std::abs(A[1][1]));
     
     // Calculate norm-infinity for A inverse
     double normAInv = std::max(std::abs(inv[0][0]) + std::abs(inv[0][1]), 
                               std::abs(inv[1][0]) + std::abs(inv[1][1]));
     
     return normA * normAInv;
 }
 
 /**
  * Performs Gaussian elimination with partial pivoting on the augmented matrix [A|b]
  * @param A The coefficient matrix
  * @param b The right-hand side vector
  * @param n The size of the system (n×n)
  * @return True if successful, false if matrix is singular
  */
 bool gaussianEliminationWithPivoting(double** A, double* b, int n) {
     // Forward elimination with partial pivoting
     for (int k = 0; k < n-1; k++) {
         // Find pivot (maximum element in current column)
         int maxRow = k;
         double maxVal = std::abs(A[k][k]);
         
         for (int i = k+1; i < n; i++) {
             if (std::abs(A[i][k]) > maxVal) {
                 maxVal = std::abs(A[i][k]);
                 maxRow = i;
             }
         }
         
         // Check for singularity (pivot too close to zero)
         if (std::abs(maxVal) < EPSILON) {
             std::cerr << "Error: Matrix is singular or nearly singular." << std::endl;
             return false;
         }
         
         // Swap rows if needed
         if (maxRow != k) {
             // Swap rows in matrix A
             double* temp = A[k];
             A[k] = A[maxRow];
             A[maxRow] = temp;
             
             // Swap elements in vector b
             double tempB = b[k];
             b[k] = b[maxRow];
             b[maxRow] = tempB;
         }
         
         // Eliminate below the pivot
         for (int i = k+1; i < n; i++) {
             double factor = A[i][k] / A[k][k];
             b[i] -= factor * b[k];
             
             for (int j = k; j < n; j++) {
                 A[i][j] -= factor * A[k][j];
             }
         }
     }
     
     // Check if the last diagonal element is close to zero
     if (std::abs(A[n-1][n-1]) < EPSILON) {
         std::cerr << "Error: Matrix is singular or nearly singular." << std::endl;
         return false;
     }
     
     return true;
 }
 
 /**
  * Performs backward substitution to find the solution vector x
  * @param A The upper triangular matrix after Gaussian elimination
  * @param b The modified right-hand side vector
  * @param n The size of the system (n×n)
  * @param x Output parameter for the solution vector
  */
 void backwardSubstitution(double** A, double* b, int n, double* x) {
     // Solve for x from bottom to top
     for (int i = n-1; i >= 0; i--) {
         double sum = 0.0;
         
         // Sum up the known terms
         for (int j = i+1; j < n; j++) {
             sum += A[i][j] * x[j];
         }
         
         // Solve for x[i]
         x[i] = (b[i] - sum) / A[i][i];
     }
 }
 
 /**
  * Writes the solution vector to a file
  * @param x The solution vector
  * @param n The size of the vector
  * @param filename The output file name
  */
 void writeOutputToFile(double* x, int n, const char* filename) {
     std::ofstream file(filename);
     
     if (!file.is_open()) {
         std::cerr << "Error: Could not open output file " << filename << std::endl;
         return;
     }
     
     file.precision(10);
     for (int i = 0; i < n; i++) {
         file << x[i] << std::endl;
     }
     
     file.close();
     std::cout << "Solution written to " << filename << std::endl;
 }
 
 /**
  * Main function that orchestrates the solution process
  * @param argc Number of command-line arguments
  * @param argv Array of command-line arguments
  * @return Exit status
  */
 int main(int argc, char* argv[]) {
     // Check if correct number of arguments provided
     if (argc != 3) {
         std::cerr << "Usage: " << argv[0] << " <matrix_A_file> <vector_b_file>" << std::endl;
         return 1;
     }
     
     // Parse command-line arguments
     const char* matrixFile = argv[1];
     const char* vectorFile = argv[2];
     
     // Read matrix A
     int rowsA, colsA;
     double** A = readMatrixFromFile(matrixFile, rowsA, colsA);
     
     // Read vector b
     int rowsB;
     double* b = readVectorFromFile(vectorFile, rowsB);
     
     // Check if A is square
     if (rowsA != colsA) {
         std::cerr << "Error: Matrix A is not square." << std::endl;
         freeMatrix(A, rowsA);
         delete[] b;
         return 1;
     }
     
     // Check if dimensions are compatible
     if (rowsA != rowsB) {
         std::cerr << "Error: Dimensions of A and b are not compatible." << std::endl;
         freeMatrix(A, rowsA);
         delete[] b;
         return 1;
     }
     
     int n = rowsA; // Size of the system
     
     // Calculate condition numbers for 2×2 matrices
     if (n == 2) {
         double cond1 = conditionNumber1(A);
         double condInf = conditionNumberInf(A);
         
         std::cout << "Matrix size: 2×2" << std::endl;
         std::cout << "Condition number (norm-1): " << cond1 << std::endl;
         std::cout << "Condition number (norm-infinity): " << condInf << std::endl;
         
         if (cond1 > 1e6 || condInf > 1e6) {
             std::cout << "Warning: High condition number detected. Results may be inaccurate." << std::endl;
         }
     }
     
     // Make copies of A and b for the actual computation (to preserve originals)
     double** A_copy = new double*[n];
     for (int i = 0; i < n; i++) {
         A_copy[i] = new double[n];
         for (int j = 0; j < n; j++) {
             A_copy[i][j] = A[i][j];
         }
     }
     
     double* b_copy = new double[n];
     for (int i = 0; i < n; i++) {
         b_copy[i] = b[i];
     }
     
     // Perform Gaussian elimination with partial pivoting
     bool success = gaussianEliminationWithPivoting(A_copy, b_copy, n);
     
     if (!success) {
         std::cerr << "Failed to solve the system due to singularity." << std::endl;
         freeMatrix(A, n);
         freeMatrix(A_copy, n);
         delete[] b;
         delete[] b_copy;
         return 1;
     }
     
     // Perform backward substitution to find the solution x
     double* x = new double[n];
     backwardSubstitution(A_copy, b_copy, n, x);
     
     // Print the solution
     std::cout << "Solution vector x:" << std::endl;
     for (int i = 0; i < n; i++) {
         std::cout << "x[" << i << "] = " << x[i] << std::endl;
     }
     
     // Write the solution to a file
     writeOutputToFile(x, n, "solution.txt");
     
     // Free all dynamically allocated memory
     freeMatrix(A, n);
     freeMatrix(A_copy, n);
     delete[] b;
     delete[] b_copy;
     delete[] x;
     
     return 0;
 }