#include "quadratic_equation.h"
#include <cmath>
#include <iostream>

QuadraticRoots solveQuadraticEquation(double a, double b, double c) {
    QuadraticRoots result;

    if (a == 0.0) {
        // Линейное уравнение bx + c = 0
        if (b == 0.0) {
            // Если b == 0, то либо бесконечное количество решений (c == 0), либо нет решений (c != 0)
            // В данном случае, возвращаем пустой вектор, сигнализируя об отсутствии решения
            return result;
        } else {
            result.roots.push_back(-c / b);
            return result;
        }
    }

    double discriminant = b * b - 4 * a * c;

    if (discriminant > 0) {
        // Два вещественных корня
        result.roots.push_back((-b + std::sqrt(discriminant)) / (2 * a));
        result.roots.push_back((-b - std::sqrt(discriminant)) / (2 * a));
    } else if (discriminant == 0) {
        // Один вещественный корень (кратный)
        result.roots.push_back(-b / (2 * a));
    } else {
        // Два комплексных корня
        double realPart = -b / (2 * a);
        double imaginaryPart = std::sqrt(-discriminant) / (2 * a);
        result.roots.push_back(std::complex<double>(realPart, imaginaryPart));
        result.roots.push_back(std::complex<double>(realPart, -imaginaryPart));
    }

    return result;
}
