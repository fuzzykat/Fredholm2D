#ifndef FE2D_SOLVER_H_
#define FE2D_SOLVER_H_
#include "Array.h"
#include "Spline2D.h"
#include "LASSolver.h"

#include <functional>

// Типы указателей на функции
typedef double (*function_1d)(double);
typedef double (*function_2d)(double, double);
typedef double (*function_4d)(double, double, double, double);

/*
 * В пространстве имен FE2D определены методы решения
 * ИУФ 2-го рода
 */
namespace FE2D
{
    // Входные данные уравнения
    struct EqData
    {
        function_4d K; // Ядро уравнения
        function_2d f; // Правая часть
        double lm;     // Коэфф. \lambda
        double a1;     // Лев. граница отрезка
        double a2;     // Прав. граница отрезка
        double b1;     // Лев. граница отрезка
        double b2;     // Прав. граница отрезка
    };
    
    /*
     * Абстрактный класс, определяющий интерфейс для методов,
     * представляющих решение в виде сплайна
     */
    class FESolver
    {
    public:
        
        // Конструктор
        // _solver -- метод решения СЛАУ
        FESolver(LA::LASSolver *_solver)
            : laSolver(_solver) {}
        
        // Виртуальный деструктор
        virtual ~FESolver() {}
        
        // Решает уравнение Фредгольма
        virtual Spline2D solve(const EqData &data) = 0;
        
        // Задает метод решения СЛАУ
        void setLASSolver(LA::LASSolver *solver)
        {
            laSolver = solver;
        }
        
    protected:
    
        // Указатель на объект класса, реализующего решение СЛАУ
        LA::LASSolver *laSolver;
    };
    
    /*
     * Метод сплайн-коллокаций с внутренними узлами Гаусса
     */
    class SCSolver : public FESolver
    {
    public:
    
        // Параметры сетки сплайна
        struct SplineOpts
        {
            // Параметры внешней сетки
            size_t xnum; // разбиения по x
            size_t ynum; // разбиения по y
            
            // Параметры внутренней сетки
            GaussGdType gtype;
        };
        
        // Конструктор
        SCSolver(const SplineOpts &_opts, LA::LASSolver *laslv)
            : FESolver(laslv), opts(_opts) {}
        
        // Реализация метода
        Spline2D solve(const EqData &data);
        
    protected:
        
        // Параметры сплайна, аппроксимирующего решение
        SplineOpts opts;
    };
    
} // namespace FE2D

#endif