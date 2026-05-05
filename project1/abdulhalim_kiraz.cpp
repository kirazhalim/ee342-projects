#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
using namespace std;

void find_upper_of_a(double* matrix_a, double* matrix_b, int rows_a, int columns_a, int rows_b) {
    
    for (int i = 0; i < rows_a -1; i++) {
        //iterate through the rows of A
        int maxRow = i;
        for (int j = i + 1; j < rows_a; j++) {
            // iterate through the rows to find the maximum element
            if (abs(matrix_a[j * columns_a + i]) > abs(matrix_a[maxRow * columns_a + i])) {
                maxRow = j;
            }
        }

        if (matrix_a[maxRow * columns_a + i] < std::numeric_limits<double>::epsilon()) {
            // check if the matrix is singular considering machine precision
            cout << "Matrix A is singular." << endl;
            delete[] matrix_a;
            delete[] matrix_b;
            return;
        }
        if (maxRow != i) {
            // check if the maximum row is not the current row, swap them
            for (int k = i; k < columns_a; k++) {
                double temp = matrix_a[i * columns_a + k];
                matrix_a[i * columns_a + k] = matrix_a[maxRow * columns_a + k];
                matrix_a[maxRow * columns_a + k] = temp;
            }
            double temp = matrix_b[i];
            matrix_b[i] = matrix_b[maxRow];
            matrix_b[maxRow] = temp;
            // swap the rows of B
        }

        for (int j = i + 1; j < rows_a; j++) {
            if (matrix_a[j * columns_a + i] != 0) {
                double multiplier = matrix_a[j * columns_a + i]/matrix_a[i * columns_a + i];
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
    cout << " Result matrix A is: " << endl;
    for (int i = 0; i < rows_a; i++) {
        for (int j = 0; j < columns_a; j++) {
            // print the matrix A
            cout << matrix_a[i * columns_a + j] << " ";
        }
        cout << endl;
    }

    cout << " Result matrix B is: " << endl;
    for (int i = 0; i < rows_a; i++) {
        // print the matrix B
            cout << matrix_b[i] << endl;

    }    
}


void find_x(double* matrix_a, double* matrix_b,double* matrix_x, int rows_a, int columns_a, int rows_b) {
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

    cout << " Result vector x is: " << endl;
    for (int i = 0; i < rows_b; i++) {
         cout << matrix_x[i] << endl;
    } 
}

void calculate_condition_numbers(double* matrix_a, int rows_a, int columns_a) {

    double norm_one = std::max(abs(matrix_a[0]) + abs(matrix_a[2]),
                                abs(matrix_a[1]) + abs(matrix_a[3]));
    // find the norm-1 of the matrix A
    double norm_inf = std::max(abs(matrix_a[0]) + abs(matrix_a[1]),
                                abs(matrix_a[2]) + abs(matrix_a[3]));
    // find the norm-infinity of the matrix A

    double det = (matrix_a[0] * matrix_a[3]) - (matrix_a[1] * matrix_a[2]);
    // calculate the determinant of the matrix A
    double inverse_det = 1.0 / det;
    // calculate the inverse of the determinant
    double inv_matrix_a[4];
    // calculate the inverse of the matrix A
    inv_matrix_a[0] = matrix_a[3] * inverse_det;
    inv_matrix_a[1] = -matrix_a[1] * inverse_det;
    inv_matrix_a[2] = -matrix_a[2] * inverse_det;
    inv_matrix_a[3] = matrix_a[0] * inverse_det;

    double norm_one_inv = std::max(abs(inv_matrix_a[0]) + abs(inv_matrix_a[2]),
                                    abs(inv_matrix_a[1]) + abs(inv_matrix_a[3]));
    double norm_inf_inv = std::max(abs(inv_matrix_a[0]) + abs(inv_matrix_a[1]),
                                    abs(inv_matrix_a[2]) + abs(inv_matrix_a[3]));
    // calculate norm-1 and norm-infinity of the inverse matrix A
    double condition_number_one = norm_one * norm_one_inv;
    double condition_number_inf = norm_inf * norm_inf_inv;
    // calculate the condition numbers
    cout << "Condition number at norm-1: " << condition_number_one << endl;
    cout << "Condition number at norm-infinity: " << condition_number_inf << endl;
}


int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "NOT MATCHED ARGUMENT NUMBERS" << endl;
        return 1;
        // checks if the number of arguments is correct
        }

    int rows_a = 0, columns_a = 0;
    int rows_b = 0;
    string line;
    double* matrix_a = nullptr;
    double* matrix_b = nullptr;
    double* matrix_x = nullptr;
    // create pointers for the matrices, since we don't know the size of the matrices
    // create null pointers for them

    ifstream matrix_a_file(argv[1]);
    // open the first input file
    ifstream matrix_b_file(argv[2]);
    // open the second input file
    ofstream matrix_x_file("result.txt");
    // open the output file


    if (matrix_a_file.is_open()) {
        // read the first file
        while (getline(matrix_a_file, line)) {
            // iterate through each line 
                stringstream ss(line);
                // read the line into a string stream to find the number of columns
                double value;
                while (ss >> value) {
                    // if ss can read a value, it means there is a column
                    columns_a++;
                    // count the number of columns in the first file
                }
            rows_a++;
            // first row counted
            break;
        }
        while (getline(matrix_a_file, line)) {
            ++rows_a;
            // count the number of remaining lines
        }
        if (rows_a == 0 || columns_a == 0) {
            // check if the file is appropriately formatted
            cout << "File " << argv[1] << " is not properly formatted." << endl;
            return 1;
        }
        
        matrix_a = new double[rows_a * columns_a];
        // allocate memory for the matrix A

        matrix_a_file.clear();
        matrix_a_file.seekg(0);
        // reset the pointer to the beginning of the file
        
        for (int i = 0; i < rows_a; i++) {
            getline(matrix_a_file, line);
            stringstream ss(line);
            for (int j = 0; j < columns_a; j++) {
                ss >> matrix_a[i * columns_a + j];
                // read the values into the matrix
            }
        }

        matrix_a_file.close();
    } else {
        cout << "Unable to open file " << argv[1] << endl;
        delete[] matrix_a;
        // free the memory if the file cannot be opened
        return 1;
    }



    if (matrix_b_file.is_open()) {
        while (getline(matrix_b_file, line)) {
                ++rows_b;
            // count the number of lines in the second file
        }
        if (rows_b == 0) {
            // check if the file is appropriately formatted
            cout << "File " << argv[2] << " is not properly formatted." << endl;
            return 1;
        }

        matrix_b = new double[rows_b];

        matrix_b_file.clear();
        matrix_b_file.seekg(0);
        // reset the pointer to the beginning of the file

        for (int i = 0; i < rows_b; i++) {
            getline(matrix_b_file, line);
            stringstream ss(line);
            ss >> matrix_b[i];
            // read the values into the matrix
        }

        matrix_b_file.close();
        
    } else {
        cout << "Unable to open file " << argv[2] << endl;
        delete[] matrix_a;
        delete[] matrix_b;
        return 1;
    }

    cout << "Matrix A is " << rows_a << "x" << columns_a << endl;
    cout << "Matrix B is " << rows_b << "x1" << endl;


    if (matrix_x_file.is_open()) {

        matrix_x = new double[rows_b];

        for(int i = 0; i < rows_b; i++) {
            matrix_x[i] = 0.0;
            // initialize the matrix x to zero
        }
    } else {
        cout << "Unable to create result.txt file" << endl;
        delete[] matrix_a;
        delete[] matrix_b;
        delete[] matrix_x;
        return 1;
    }


    if (rows_a == 2 && columns_a == 2) {
        // if the matrix is 2x2, calculate the condition numbers
        calculate_condition_numbers(matrix_a,rows_a,columns_a);
    } else {
        cout << "The matrix is not 2x2, thus the condition number is not calculated." << endl;
    }

    
    find_upper_of_a(matrix_a, matrix_b, rows_a, columns_a, rows_b); 
    // perform Gaussian Elimination with partial pivoting
    find_x(matrix_a, matrix_b, matrix_x, rows_a, columns_a, rows_b);
    // perform back substitution to find the solution

    for(int i = 0; i < rows_b; i++) {
        // write the solution x to the output file
        matrix_x_file << matrix_x[i] << endl;
    }
    matrix_x_file.close();

    delete[] matrix_a;
    delete[] matrix_b;
    delete[] matrix_x;
    //free the memory

    return 0;
}