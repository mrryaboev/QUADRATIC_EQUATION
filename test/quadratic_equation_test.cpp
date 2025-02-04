#include "gtest/gtest.h"
#include "quadratic_equation.h"
#include <limits>
#include <cmath>

// Функция для сравнения double с учетом погрешности
bool compareDouble(double a, double b, double epsilon = 1e-6) {
    return std::abs(a - b) < epsilon;
}

// Функция для сравнения комплексных чисел
bool compareComplex(const std::complex<double>& a, const std::complex<double>& b, double epsilon = 1e-6) {
    return compareDouble(a.real(), b.real(), epsilon) && compareDouble(a.imag(), b.imag(), epsilon);
}

TEST(QuadraticEquationTest, TwoRealRoots) {
    QuadraticRoots roots = solveQuadraticEquation(1, -3, 2);
    ASSERT_EQ(roots.roots.size(), 2);
    EXPECT_TRUE(compareDouble(roots.roots[0].real(), 2.0));
    EXPECT_TRUE(compareDouble(roots.roots[1].real(), 1.0));
}

TEST(QuadraticEquationTest, OneRealRoot) {
    QuadraticRoots roots = solveQuadraticEquation(1, -2, 1);
    ASSERT_EQ(roots.roots.size(), 1);
    EXPECT_TRUE(compareDouble(roots.roots[0].real(), 1.0));
}

TEST(QuadraticEquationTest, ComplexRoots) {
    QuadraticRoots roots = solveQuadraticEquation(1, 2, 5);
    ASSERT_EQ(roots.roots.size(), 2);
    std::complex<double> expectedRoot1(-1.0, 2.0);
    std::complex<double> expectedRoot2(-1.0, -2.0);
    EXPECT_TRUE(compareComplex(roots.roots[0], expectedRoot1));
    EXPECT_TRUE(compareComplex(roots.roots[1], expectedRoot2));
}

TEST(QuadraticEquationTest, LinearEquation) {
    QuadraticRoots roots = solveQuadraticEquation(0, 2, -4);
    ASSERT_EQ(roots.roots.size(), 1);
    EXPECT_TRUE(compareDouble(roots.roots[0].real(), 2.0));
}

TEST(QuadraticEquationTest, NoRoots) {
    QuadraticRoots roots = solveQuadraticEquation(0, 0, 5);
    ASSERT_EQ(roots.roots.size(), 0);
}

TEST(QuadraticEquationTest, ZeroCoefficientA) {
    QuadraticRoots roots = solveQuadraticEquation(0, 5, 2);
    ASSERT_EQ(roots.roots.size(), 1);
    EXPECT_TRUE(compareDouble(roots.roots[0].real(), -0.4));
}

TEST(QuadraticEquationTest, LargeCoefficients) {
    double a = 1e10;
    double b = 2e10;
    double c = 1e10;
    QuadraticRoots roots = solveQuadraticEquation(a, b, c);
    ASSERT_EQ(roots.roots.size(), 1);
    EXPECT_TRUE(compareDouble(roots.roots[0].real(), -1.0));
}

TEST(QuadraticEquationTest, SmallCoefficients) {
    double a = 1e-10;
    double b = 2e-10;
    double c = 1e-10;
    QuadraticRoots roots = solveQuadraticEquation(a, b, c);
    ASSERT_EQ(roots.roots.size(), 1);
    EXPECT_TRUE(compareDouble(roots.roots[0].real(), -1.0));
}

TEST(QuadraticEquationTest, ExtremeValues) {
    double a = std::numeric_limits<double>::max();
    double b = std::numeric_limits<double>::max();
    double c = std::numeric_limits<double>::max();
    QuadraticRoots roots = solveQuadraticEquation(a, b, c);
    ASSERT_EQ(roots.roots.size(), 1);
    EXPECT_TRUE(compareDouble(roots.roots[0].real(), -1.0));
}

TEST(QuadraticEquationTest, NegativeDiscriminant) {
    QuadraticRoots roots = solveQuadraticEquation(1, 1, 1);
    ASSERT_EQ(roots.roots.size(), 2);
    EXPECT_TRUE(compareComplex(roots.roots[0], std::complex<double>(-0.5, std::sqrt(3.0)/2.0)));
    EXPECT_TRUE(compareComplex(roots.roots[1], std::complex<double>(-0.5, -std::sqrt(3.0)/2.0)));
}

TEST(QuadraticEquationTest, ZeroCoefficients) {
    QuadraticRoots roots = solveQuadraticEquation(0, 0, 0);
    ASSERT_EQ(roots.roots.size(), 0);
}
