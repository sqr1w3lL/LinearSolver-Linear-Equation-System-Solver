#include <iostream>
#include <locale>
#include <windows.h>
#include <iomanip>
#include <cmath>
using namespace std;

// ---------- Метод Гауса ----------
void gauss(double** A, double* B, double* X, int n) {
    for (int i = 0; i < n; i++) {
        double diag = A[i][i];
        for (int j = i; j < n; j++) A[i][j] /= diag;
        B[i] /= diag;

        for (int k = 0; k < n; k++) {
            if (k == i) continue;
            double factor = A[k][i];
            for (int j = i; j < n; j++) A[k][j] -= factor * A[i][j];
            B[k] -= factor * B[i];
        }
    }
    for (int i = 0; i < n; i++) X[i] = B[i];
}

// ---------- Обчислення визначника ----------
double determinant(double** M, int n) {
    if (n == 1) return M[0][0];
    if (n == 2) return M[0][0] * M[1][1] - M[0][1] * M[1][0];
    double det = 0;
    for (int p = 0; p < n; p++) {
        double** sub = new double* [n - 1];
        for (int i = 0; i < n - 1; i++) sub[i] = new double[n - 1];

        for (int i = 1; i < n; i++) {
            int colIndex = 0;
            for (int j = 0; j < n; j++) {
                if (j == p) continue;
                sub[i - 1][colIndex++] = M[i][j];
            }
        }
        det += M[0][p] * pow(-1, p) * determinant(sub, n - 1);

        for (int i = 0; i < n - 1; i++) delete[] sub[i];
        delete[] sub;
    }
    return det;
}

// ---------- Метод Крамера ----------
void cramer(double** A, double* B, double* X, int n) {
    double detA = determinant(A, n);
    for (int i = 0; i < n; i++) {
        double** Ai = new double* [n];
        for (int j = 0; j < n; j++) {
            Ai[j] = new double[n];
            for (int k = 0; k < n; k++) Ai[j][k] = A[j][k];
            Ai[j][i] = B[j];
        }
        X[i] = determinant(Ai, n) / detA;
        for (int j = 0; j < n; j++) delete[] Ai[j];
        delete[] Ai;
    }
}

// ---------- Матричний метод ----------
void matrixInverseSolve(double** A, double* B, double* X, int n) {
    double** aug = new double* [n];
    for (int i = 0; i < n; i++) {
        aug[i] = new double[2 * n];
        for (int j = 0; j < n; j++) aug[i][j] = A[i][j];
        for (int j = n; j < 2 * n; j++) aug[i][j] = (j - n == i) ? 1 : 0;
    }

    for (int i = 0; i < n; i++) {
        double diag = aug[i][i];
        for (int j = 0; j < 2 * n; j++) aug[i][j] /= diag;
        for (int k = 0; k < n; k++) {
            if (k == i) continue;
            double factor = aug[k][i];
            for (int j = 0; j < 2 * n; j++) aug[k][j] -= factor * aug[i][j];
        }
    }

    for (int i = 0; i < n; i++) {
        X[i] = 0;
        for (int j = 0; j < n; j++) X[i] += aug[i][n + j] * B[j];
    }

    for (int i = 0; i < n; i++) delete[] aug[i];
    delete[] aug;
}

// ---------- Вивід системи ----------
void printSystem(double** A, double* B, int n) {
    cout << "Система рівнянь:" << endl;
    for (int i = 0; i < n; i++) {
        cout << fixed << setprecision(2);
        for (int j = 0; j < n; j++) {
            if (j > 0 && A[i][j] >= 0) cout << " + ";
            else if (A[i][j] < 0) cout << " - ";
            cout << fabs(A[i][j]) << "x" << j + 1;
        }
        cout << " = " << B[i] << endl;
    }
    cout << endl;
}

int main() {
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
    setlocale(LC_ALL, "ukr");

    int n = 4;

    double initA[4][4] = {
        {0.63, 1, 0.11, 0.34},
        {0.17, 1.18, -0.45, 0.11},
        {0.31, -0.15, 1.17, -2.35},
        {0.58, 0.21, -3.45, -1.18}
    };
    double initB[4] = { 2.08, 0.17, 1.28, 0.05 };


    double** Ashow = new double* [n];
    for (int i = 0; i < n; i++) {
        Ashow[i] = new double[n];
        for (int j = 0; j < n; j++) Ashow[i][j] = initA[i][j];
    }
    printSystem(Ashow, initB, n);

    bool again = true;
    while (again) {
        cout << "Оберіть метод розв’язання:\n";
        cout << "1. Метод Гауса\n";
        cout << "2. Метод Крамера\n";
        cout << "3. Матричний метод\n";
        cout << "Ваш вибір: ";
        int choice;
        cin >> choice;

        double* X = new double[n];

        if (choice == 1) {
            double** A1 = new double* [n];
            double* B1 = new double[n];
            for (int i = 0; i < n; i++) {
                A1[i] = new double[n];
                for (int j = 0; j < n; j++) A1[i][j] = initA[i][j];
                B1[i] = initB[i];
            }
            gauss(A1, B1, X, n);
            cout << "Розв’язок (Гаус): ";
            for (int i = 0; i < n; i++) cout << fixed << setprecision(3) << X[i] << " ";
            cout << endl;
            for (int i = 0; i < n; i++) delete[] A1[i];
            delete[] A1; delete[] B1;
        }
        else if (choice == 2) {
            double** A2 = new double* [n];
            for (int i = 0; i < n; i++) {
                A2[i] = new double[n];
                for (int j = 0; j < n; j++) A2[i][j] = initA[i][j];
            }
            cramer(A2, initB, X, n);
            cout << "Розв’язок (Крамер): ";
            for (int i = 0; i < n; i++) cout << fixed << setprecision(3) << X[i] << " ";
            cout << endl;
            for (int i = 0; i < n; i++) delete[] A2[i];
            delete[] A2;
        }
        else if (choice == 3) {
            double** A3 = new double* [n];
            for (int i = 0; i < n; i++) {
                A3[i] = new double[n];
                for (int j = 0; j < n; j++) A3[i][j] = initA[i][j];
            }
            matrixInverseSolve(A3, initB, X, n);
            cout << "Розв’язок (Матричний метод): ";
            for (int i = 0; i < n; i++) cout << fixed << setprecision(3) << X[i] << " ";
            cout << endl;
            for (int i = 0; i < n; i++) delete[] A3[i];
            delete[] A3;
        }
        else {
            cout << "Невірний вибір!" << endl;
        }

        delete[] X;

        cout << "\nБажаєте вирішити іншим методом? (1 - так, 0 - ні): ";
        int ans;
        cin >> ans;
        again = (ans == 1);
    }

    for (int i = 0; i < n; i++) delete[] Ashow[i];
    delete[] Ashow;

    cout << "Програма завершена." << endl;
    return 0;
}
