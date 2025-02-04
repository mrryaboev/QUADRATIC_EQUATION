#ifndef QUADRATIC_EQUATION_H
#define QUADRATIC_EQUATION_H

#include <vector>
#include <complex>

// Структура для хранения результатов
struct QuadraticRoots {
    std::vector<std::complex<double>> roots;
};

// Функция для решения квадратного уравнения
QuadraticRoots solveQuadraticEquation(double a, double b, double c);

#endif
