#ifndef SPLINE_2D_H_
#define SPLINE_2D_H_
#include "Array.h"
#include "GaussNodes.h"

/*
 * Класс, реализующий двумерный интерполяционный полином Лагранжа
 */
class LagrangePoly2D
{
public:

    LagrangePoly2D() {}

    LagrangePoly2D(
        const array_1d &_xnodes,
        const array_1d &_ynodes,
        const array_2d &_values)
    : xnodes(_xnodes), ynodes(_ynodes), values(_values) {}
        
    ~LagrangePoly2D() {}
    
    const array_1d &getXNodes() const
    {
        return xnodes;
    }
    
    const array_1d &getYNodes() const
    {
        return ynodes;
    }
    
    void setXNodes(const array_1d &xn)
    {
        xnodes = xn;
    }
    
    void setYNodes(const array_1d &yn)
    {
        ynodes = yn;
    }
    
    void setValues(const array_2d &vs)
    {
        values = vs;
    }
    
    // Вычисляет значение полинома в точке (x, y)
    double operator()(double x, double y) const;

    // Вычисляет i-й полином по переменной x
    double xFundPoly(double x, size_t i)
    {
        return fundPoly(xnodes, x, i);
    }
    
    // Вычисляет j-й полином по переменной y
    double yFundPoly(double y, size_t j)
    {
        return fundPoly(ynodes, y, j);
    }
    
private:

    // Вычисляет k-й фундаментальный полином в точке x
    double fundPoly(const array_1d &g, double x, size_t k) const;

    array_1d xnodes; // Узлы по переменной x
    array_1d ynodes; // Узлы по переменной y
    array_2d values; // Значения в узлах
};

/*
 * Класс, реализующий двумерный сплайн,
 * состоящий из полиномов Лагранжа
 */
class Spline2D
{
public:
    
    typedef valarray_2d<LagrangePoly2D *> array_2p;
    
    // Конструктор
    Spline2D(
        const array_1d &_xgrid, // Разбиение области по переменной x
        const array_1d &_ygrid, // Разбиение области по переменной y
        const array_2p &_polys) // Массив полиномов Лагранжа
    : xgrid(_xgrid), ygrid(_ygrid), polys(_polys) {}
    
    // Конструктор копирования
    Spline2D(const Spline2D &spl)
    {
        *this = spl;
    }
    
    // Деструктор
    ~Spline2D();
    
    // Оператор присваивания
    const Spline2D &operator=(const Spline2D &spl);
    
    array_1d &getXGrid()
    {
        return xgrid;
    }
    
    array_1d &getYGrid()
    {
        return ygrid;
    }
    
    array_2p &getPolynoms()
    {
        return polys;
    }
    
    // Вычисляет значение сплайна в точке (x, y)
    double operator()(double x, double y) const;
    
private:

    // Сетка узлов по переменной x
    array_1d xgrid;
    
    // Сетка узлов по переменной y
    array_1d ygrid;
    
    // Массив интерполяционных полиномов
    array_2p polys;
};

// Создает неициализированный экземпляр класса Spline2D с сеткой Гаусса
Spline2D createGaussSpline2D(
    const double xval[2],  // Границы области по x
    const double yval[2],  // Границы области по y
    size_t xnum,           // Число разбиений по переменной x
    size_t ynum,           // Число разбиений по переменной y
    GaussGdType gtype = G7 // Число узлов Гаусса в полиномах Лагранжа
);

#endif