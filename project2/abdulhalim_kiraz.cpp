#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

template<typename T>
class Matrix {
    public:
    int rows, cols;
    T* matrix;
    

    Matrix() : rows(0), cols(0), matrix(nullptr) {}
    // create default constructor for empty matrix
    Matrix(int rows, int cols)
    // create matrix array according to the size of matrix 
    : rows(rows), cols(cols), matrix(nullptr) {
        if (rows > 0 && cols > 0) {
            matrix = new T[rows * cols];
        }
    }

    ~Matrix() {
        // create destructor to free the storage
        delete[] matrix;
    }

    int getMatrixFromFile(const string& filename) {
        rows = 0, cols = 0;
        delete[] matrix;
        matrix = nullptr;
        // free the memory if the matrix is not empty

        ifstream file(filename);
        // open the file
        if (!file) {
            // raise an error message if there is something wrong with the file
            cout << filename << " couldn't open." << endl;
            delete[] matrix;
            return 1;
        }

        string line;
        bool first_row = true;
        // check if the first row is empty

        while (getline(file, line)) {
            // read the lines of the file
            stringstream ss(line);
            double value;
            int cols_first_row = 0;

            while (ss >> value) {
                ++cols_first_row;
            }

            if (cols_first_row == 0) {
                continue;
            }
            if (first_row) {
                cols   = cols_first_row;
                first_row = false;
            } else if (cols_first_row != cols) {
                cout << filename << " is not a rectangular matrix.\n";
                return 1;
            }
            ++rows; // count the rows
        }
        if (rows == 0 || cols == 0) {
            // check if the matrix is valid
            cout << filename << " is not a valid matrix." << endl;
            return 1;
        }

        matrix = new T[rows * cols];
        // allocate memory for the matrix A

        file.clear();
        file.seekg(0);
        // reset the pointer to the beginning of the file

        for (int i = 0; i < rows; i++) {
            // iterate through rows of matrix
            getline(file, line);
            // seperate the lines
            stringstream ss(line);
            // turn the lines into streams to read each value
            for (int j = 0; j < cols; j++) {
                // iterate through columns of matrix
                ss >> matrix[i * cols + j];
                // read the values into the matrix
            }
        }
        file.close();
        // close the file
        return 0;    
    }
};

template<typename T>
void find_upper_of_a(T* matrix_a, T* matrix_b, int rows_a, int columns_a, int rows_b) {

    for (int i = 0; i < rows_a -1; i++) {
        //iterate through the rows of A
        int maxRow = i;
        for (int j = i + 1; j < rows_a; j++) {
            // iterate through the rows to find the maximum element
            if (abs(matrix_a[j * columns_a + i]) > abs(matrix_a[maxRow * columns_a + i])) {
                maxRow = j;
            }
        }

        if (fabs(matrix_a[maxRow * columns_a + i]) < 1e-9) {
            // check if the matrix is singular considering machine precision
            cout << "Matrix A is singular." << endl;
            return;
        }
        if (maxRow != i) {
            // check if the maximum row is not the current row, swap them
            for (int k = i; k < columns_a; k++) {
                T temp = matrix_a[i * columns_a + k];
                matrix_a[i * columns_a + k] = matrix_a[maxRow * columns_a + k];
                matrix_a[maxRow * columns_a + k] = temp;
            }
            T temp = matrix_b[i];
            matrix_b[i] = matrix_b[maxRow];
            matrix_b[maxRow] = temp;
            // swap the rows of B
        }

        for (int j = i + 1; j < rows_a; j++) {
            if (matrix_a[j * columns_a + i] != 0) {
                T multiplier = matrix_a[j * columns_a + i]/matrix_a[i * columns_a + i];
                //find the factor for subtraction
                for (int k = i; k < columns_a; k++) {
                    matrix_a[j * columns_a + k] = matrix_a[j * columns_a + k] - multiplier*matrix_a[i * columns_a + k];
                //perform subtraction to make the pivot zero
                }

            matrix_b[j] = matrix_b[j] - multiplier*matrix_b[i];
            // update the right-hand side
            }
        }
    } 
}

template<typename T>
void find_x(T* matrix_a, T* matrix_b,T* matrix_x, int rows_a, int columns_a, int rows_b) {
    // this function is where we applied back substitution
    for (int j = rows_b -1 ; j >= 0; j--) {
        //loop backwards over columns
        matrix_x[j] = matrix_b[j]/matrix_a[j * columns_a + j];
        // compute solution
        for (int i = j - 1; i >= 0; i--) {
            matrix_b[i] = matrix_b[i] - matrix_a[i * columns_a + j]*matrix_x[j];
            // update right-hand side
        }
    }
}

// create a function to perform inverse iteration algorithm to compute the smallest eigenvalue
template<typename T>
void inverseIteration(const Matrix<T>& matrix, double tolerance, T& lambda, T* eigenvector) {
    if (matrix.rows != matrix.cols) {
        cout << "Matrix is not square" << endl;
        return;
    }
    // check if the matrix is square
    if (matrix.rows == 0) {
        cout << "Matrix is empty" << endl;
        return;
    }

    int n = matrix.rows;

    T* x = new T[n];
    T* y = new T[n];
    T* y_next = new T[n];
    T* a_copy = new T[n * n];
    T* b = new T[n];
    // create pointers

    for (int i = 0; i < n; i++) {
        x[i] = 0.0;
    }
    x[n-1] = 1.0;
    // choose the x vector

    int count = 0;
    while (true) {

        copy(matrix.matrix, matrix.matrix + n * n, a_copy);
        // copy the matrix to a new matrix
        copy(x, x + n, b);
        // copy the x vector to a new vector

        find_upper_of_a(a_copy, b, n, n, n);
        find_x(a_copy, b, y, n, n, n);
        
        T norm = 0.0;
        for (int i = 0; i < n; i++) {
            if (abs(y[i]) > norm) {
                norm = abs(y[i]);
            }
        }
        // calculate the inf-norm of the y vector

        if (norm == 0.0) {
            cout << "Matrix is singular" << endl;
            delete[] x;
            delete[] y;
            delete[] y_next;
            delete[] a_copy;
            delete[] b;
            return;
        }
        // check if the matrix is singular

        for (int i = 0; i < n; i++) {
            y_next[i] = y[i] / norm;
        }
        // normalize the y vector

        T diff = 0.0;
        for (int i = 0; i < n; i++) {
            diff = max(diff, abs(y_next[i]-x[i]));
        }
        // calculate the difference between the two vectors

        if (diff < tolerance) {
            break;
        }
        // check if the difference is less than the tolerance

        if (count > 1000) {
            cout << "Maximum iterations reached" << endl;
            break;
        }
        count++;
        // if the program gets into a loop, break the loop

        for (int i = 0; i < n; i++) {
            x[i] = y_next[i];
        }
        // update the x vector


    }

    T numerator = 0.0, denominator = 0.0;
    for (int i = 0; i < n; i++) {
        T Ax_i = 0.0;
        for (int j = 0; j < n; j++) {
            Ax_i += matrix.matrix[i * n + j] * y_next[j];
        }
        numerator += y_next[i] * Ax_i;
        denominator += y_next[i] * y_next[i];
    }
    lambda = numerator / denominator;
    // calculate the eigenvalue using the rayleigh quotient

    T factor = y_next[n-1];
    for (int i = 0; i < n; i++) {
        y_next[i] /= factor;
    }
    // normalize the eigenvector

    for (int i = 0; i < n; i++) {
        eigenvector[i] = y_next[i];
    }

    delete[] x;
    delete[] y;
    delete[] y_next;
    delete[] a_copy;
    delete[] b;
    // free the memory
}

template<typename T>
// create a function to apply normalized power iteration
void normalizedPowerIteration(const Matrix<T>& matrix, double tolerance, T& lambda, T* eigenvector) {
    
    if (matrix.rows != matrix.cols) {
        cout << "Matrix is not square" << endl;
        return;
    }
    // check if the matrix is square

    int n = matrix.rows;

    T* x = new T[n];
    T* y = new T[n];
    T* y_next = new T[n];

    for (int i = 0; i < n -1; i++) {
        x[i] = 0.0;
    }
    x[n-1] = 1.0;
    // choose the x vector

    int count = 0;
    while (true) {
        for (int i = 0; i < n; i++) {
            y[i] = 0.0;
        }
        // set the y vector to be 0.0

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                y[i] += matrix.matrix[i * n + j] * x[j];
            }
        }
        // multiply the matrix with the x vector

        T norm = 0.0;
        for (int i = 0; i < n; i++) {
            if (abs(y[i]) > norm) {
                norm = abs(y[i]);
            }
        }
        // calculate the inf-norm of the y vector

        if (norm == 0.0) {
            cout << "Matrix is singular" << endl;
            delete[] x;
            delete[] y;
            delete[] y_next;
            return;
        }
        // check if the matrix is singular

        for (int i = 0; i < n; i++) {
            y_next[i] = y[i] / norm;
        }
        // normalize the y vector

        T diff = 0.0;
        for (int i = 0; i < n; i++) {
            diff = max(diff, abs(y_next[i]-x[i]));
        }
        // calculate the difference between the two vectors

        if (diff < tolerance) {
            break;
        }

        if (count > 1000) {
            cout << "Maximum iterations reached" << endl;
            break;
        }
        // if the program gets into a loop, break the loop
        count++;


        for (int i = 0; i < n; i++) {
            x[i] = y_next[i];
        }
        // update the x vector


    }

    // using the rayleigh quotient to calculate the eigenvalue
    T numerator = 0.0;
    T denominator = 0.0;
    for (int i = 0; i < n; i++) {
        T Ax_i = 0.0;
        for (int j = 0; j < n; j++) {
            Ax_i += matrix.matrix[i * n + j] * y_next[j];
        }
        numerator += y_next[i] * Ax_i;
        denominator += y_next[i] * y_next[i];
    }
    lambda = numerator / denominator;

    T factor = y_next[n-1];
    for (int i = 0; i < n; i++) {
        y_next[i] /= factor;
    }
    // normalize the eigenvector

    for (int i = 0; i < n; i++) {
        eigenvector[i] = y_next[i];
    }
    
    delete[] x;
    delete[] y;
    delete[] y_next;
    // free the memory
    
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        cout << "NOT MATCHED ARGUMENT NUMBERS" << endl;
        return 1;
        // checks if the number of arguments is correct 
    }

    
    Matrix<double> matrix_a;
    Matrix<double> matrix_a_inv;
    // create matrix objects
    matrix_a.getMatrixFromFile(argv[1]);
    matrix_a_inv.getMatrixFromFile(argv[2]);
    // read the matrices from the files

    if (matrix_a.rows != matrix_a_inv.rows || matrix_a.cols != matrix_a_inv.cols) {
        // check if the matrices are the same size
        cout << "Matrices are not the same size" << endl;
        return 1;
    }

    istringstream ss(argv[3]);
    // read tolerance from argv[3]
    double tolerance;
    if (!(ss >> tolerance)) {
        // check if the input tolerance is valid
        cout << "Invalid tolerance" << endl;
        return 1;
    }
    if (tolerance < 0) {
        // check if the tolerance is negative
        cout << "Tolerance must be a positive number" << endl;
        return 1;
    }

    ofstream outFile(argv[4]);
    // open the output file
    if (!outFile) {
        // check if the output file is valid
        cout << "Output file couldn't open." << endl;
        return 1;
    }

    double lambda1;
    double* eigenvector_A = new double[matrix_a.rows];

    inverseIteration(matrix_a, tolerance, lambda1, eigenvector_A);
    // compute the smallest eigenvalue of A and its eigenvector with inverse iteration
    outFile << "Smallest eigenvalue of A in magnitude: " << lambda1 << endl;
    outFile << "Corresponding eigenvector: " << endl;
    for (int i = 0; i < matrix_a.rows; i++) {
        outFile << eigenvector_A[i] << endl;
    }

    double lambda2;
    double* eigenvector_A_inv = new double[matrix_a_inv.rows];
    normalizedPowerIteration(matrix_a_inv, tolerance, lambda2, eigenvector_A_inv);
    // compute the largest eigenvalue of A and its eigenvector with normalized power iteration
    outFile << "Dominant eigenvalue of A_inv in magnitude: " << lambda2 << endl;
    outFile << "Corresponding eigenvector: " << endl;
    for (int i = 0; i < matrix_a_inv.rows; i++) {
        outFile << eigenvector_A_inv[i] << endl;
    }

    outFile << "1.0/" << lambda2 << "=" << 1.0/lambda2 << endl;
    // verify that the dominant eigenvalue of A inverse is the smallest eigenvalue of A

    delete[] eigenvector_A;
    delete[] eigenvector_A_inv;
    outFile.close();
    // free the memory

    return 0;

}
