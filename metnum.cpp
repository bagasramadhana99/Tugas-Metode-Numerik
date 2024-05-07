// Bagas Ramadhana
// 21120122120001
// kelas D
// Tugas Implementasi Sistem Persamaan Linier

#include <iostream>
#include <vector>

using namespace std;

// Fungsi untuk mentransposisi matriks
vector<vector<double>> transpose(const vector<vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    vector<vector<double>> result(cols, vector<double>(rows));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[j][i] = matrix[i][j];
        }
    }
    return result;
}

// Fungsi untuk perkalian matriks
vector<vector<double>> matrix_multiply(const vector<vector<double>>& matrix1, const vector<vector<double>>& matrix2) {
    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int cols2 = matrix2[0].size();
    vector<vector<double>> result(rows1, vector<double>(cols2, 0.0));
    for (int i = 0; i < rows1; ++i) {
        for (int j = 0; j < cols2; ++j) {
            for (int k = 0; k < cols1; ++k) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return result;
}

// Fungsi untuk menghitung matriks balikan
vector<vector<double>> matrix_inverse(const vector<vector<double>>& matrix) {
    double det = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
                 matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
                 matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
    double invdet = 1 / det;
    vector<vector<double>> inverse(3, vector<double>(3, 0.0));
    inverse[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) * invdet;
    inverse[0][1] = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) * invdet;
    inverse[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) * invdet;
    inverse[1][0] = (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) * invdet;
    inverse[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) * invdet;
    inverse[1][2] = (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) * invdet;
    inverse[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) * invdet;
    inverse[2][1] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) * invdet;
    inverse[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) * invdet;
    return inverse;
}

// Fungsi untuk melakukan dekomposisi LU
pair<vector<vector<double>>, vector<vector<double>>> lu_decomposition(const vector<vector<double>>& matrix) {
    int size = matrix.size();
    vector<vector<double>> L(size, vector<double>(size, 0.0));
    vector<vector<double>> U(size, vector<double>(size, 0.0));
    for (int i = 0; i < size; ++i) {
        L[i][i] = 1.0;
        for (int j = i; j < size; ++j) {
            double total = 0.0;
            for (int k = 0; k < i; ++k) {
                total += L[i][k] * U[k][j];
            }
            U[i][j] = matrix[i][j] - total;
        }
        for (int j = i + 1; j < size; ++j) {
            double total = 0.0;
            for (int k = 0; k < i; ++k) {
                total += L[j][k] * U[k][i];
            }
            L[j][i] = (matrix[j][i] - total) / U[i][i];
        }
    }
    return make_pair(L, U);
}

// Fungsi untuk melakukan dekomposisi Crout
pair<vector<vector<double>>, vector<vector<double>>> crout_decomposition(const vector<vector<double>>& matrix) {
    int size = matrix.size();
    vector<vector<double>> L(size, vector<double>(size, 0.0));
    vector<vector<double>> U(size, vector<double>(size, 0.0));
    for (int i = 0; i < size; ++i) {
        U[i][i] = 1.0;
        for (int j = i; j < size; ++j) {
            double sum = 0.0;
            for (int k = 0; k < i; ++k) {
                sum += L[i][k] * U[k][j];
            }
            L[j][i] = matrix[j][i] - sum;
        }
        for (int j = i + 1; j < size; ++j) {
            double sum = 0.0;
            for (int k = 0; k < i; ++k) {
                sum += L[j][k] * U[k][i];
            }
            U[i][j] = (matrix[i][j] - sum) / L[i][i];
        }
    }
    return make_pair(L, U);
}

// Fungsi untuk menyelesaikan sistem persamaan linear menggunakan matriks balikan
vector<double> solve_using_inverse(const vector<vector<double>>& A, const vector<double>& b) {
    vector<vector<double>> A_inv = matrix_inverse(A);
    vector<vector<double>> b_matrix = {b};
    vector<vector<double>> x_matrix = matrix_multiply(A_inv, transpose(b_matrix));
    vector<double> x = {x_matrix[0][0], x_matrix[1][0], x_matrix[2][0]};
    return x;
}

// Fungsi untuk menyelesaikan sistem persamaan linear menggunakan dekomposisi LU
vector<double> solve_using_lu(const vector<vector<double>>& A, const vector<double>& b) {
    auto LU = lu_decomposition(A);
    const vector<vector<double>>& L = LU.first;
    const vector<vector<double>>& U = LU.second;
    int size = A.size();
    vector<double> y(size, 0.0);
    vector<double> x(size, 0.0);

    // Solve Ly = b
    for (int i = 0; i < size; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }

    // Solve Ux = y
    for (int i = size - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < size; ++j) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }
    return x;
}

// Fungsi untuk menyelesaikan sistem persamaan linear menggunakan dekomposisi Crout
vector<double> solve_using_crout(const vector<vector<double>>& A, const vector<double>& b) {
    auto LU = crout_decomposition(A);
    const vector<vector<double>>& L = LU.first;
    const vector<vector<double>>& U = LU.second;
    int size = A.size();
    vector<double> y(size, 0.0);
    vector<double> x(size, 0.0);

    // Solve Ly = b
    for (int i = 0; i < size; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }

    // Solve Ux = y
    for (int i = size - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < size; ++j) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum);
    }
    return x;
}

// Fungsi untuk memasukkan matriks dari stdin
vector<vector<double>> input_matrix(int size) {
    vector<vector<double>> matrix(size, vector<double>(size));
    cout << "Masukkan elemen-elemen matriks A:" << endl;
    for (int i = 0; i < size; ++i) {
        cout << "Baris ke-" << i + 1 << ": ";
        for (int j = 0; j < size; ++j) {
            cin >> matrix[i][j];
        }
    }
    return matrix;
}

// Fungsi untuk memasukkan vektor dari stdin
vector<double> input_vector(int size) {
    vector<double> vector(size);
    cout << "Masukkan elemen-elemen vektor b:" << endl;
    for (int i = 0; i < size; ++i) {
        cin >> vector[i];
    }
    return vector;
}

// Fungsi utama
int main() {
    cout << "Masukkan ukuran matriks A: ";
    int n;
    cin >> n;

    // Masukkan matriks A
    vector<vector<double>> A = input_matrix(n);

    // Masukkan vektor b
    cout << "Masukkan ukuran vektor b: ";
    vector<double> b = input_vector(n);

    // Solusi menggunakan matriks balikan
    vector<double> x_inverse = solve_using_inverse(A, b);
    cout << "\nSolusi menggunakan metode matriks balikan:" << endl;
    for (int i = 0; i < x_inverse.size(); ++i) {
        cout << "x[" << i + 1 << "] = " << x_inverse[i] << endl;
    }
    cout << endl;

    // Solusi menggunakan dekomposisi LU
    vector<double> x_lu = solve_using_lu(A, b);
    cout << "Solusi menggunakan metode dekomposisi LU:" << endl;
    for (int i = 0; i < x_lu.size(); ++i) {
        cout << "x[" << i + 1 << "] = " << x_lu[i] << endl;
    }
    cout << endl;

    // Solusi menggunakan dekomposisi Crout
    vector<double> x_crout = solve_using_crout(A, b);
    cout << "Solusi menggunakan metode dekomposisi Crout:" << endl;
    for (int i = 0; i < x_crout.size(); ++i) {
        cout << "x[" << i + 1 << "] = " << x_crout[i] << endl;
    }

    return 0;
}
