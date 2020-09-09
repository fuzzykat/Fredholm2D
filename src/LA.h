#ifndef LA_H_
#define LA_H_
#include "Array.h"

/*
 * Пространство имен LA (Linear Algebra)
 * Здесь определены классы матриц и векторов и
 * операции с ними
 */
namespace LA
{
    typedef array_1d Vector;
    typedef array_2d Matrix;
    
    // Прибавить константу к главной диагонали
    void cmAdd(Matrix &A, double alpha);
    
    // Умножение вектора на константу
    void cvMult(Vector &A, double alpha);
    
    // Умножение матрицы на константу
    void cmMult(Matrix &A, double alpha);
    
    // Умножение матрицы на вектор
    void mvMult(const Matrix &A, Vector &V);
    
    // Умножение матриц
    void mmMult(const Matrix &A, Matrix &B);
    
    // Сумма векторов
    void vvAdd(const Vector &A, Vector &V);
    
    // Разность векторов
    void vvSub(const Vector &A, Vector &V);
    
    // Транспонирование матрицы A: A = A^T
    void mTranspose(Matrix &A);
    
    // Вычисляет норму матрицы A
    double mNorm(const Matrix &A);
    
    // Вычисляет норму вектора A
    double vNorm(const Vector &A);
    
    // Вычисляет скалярное произведение
    double vvProd(const Vector &A, const Vector &B);
    
    // Вычисляет максимальный модуль разности элементов векторов 
    inline double vvMaxAbsDiff(const Vector &A, const Vector &B)
    {
        return std::abs(A - B).max();
    }
    
} // namespace LA

#endif