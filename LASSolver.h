#ifndef LAS_SOLVER_H_
#define LAS_SOLVER_H_
#include "Array.h"

/*
 * Пространство имен LA (Linear Algebra)
 * Здесь определены классы, реализующие 
 * численные методы решения СЛАУ
 */
namespace LA
{
    /*
     * Абстрактный класс, определяющий интерфейс для
     * методов решения СЛАУ
     */
    class LASSolver
    {
    public:
        
        // Виртуальный деструктор
        virtual ~LASSolver() {}
        
        // Решает СЛАУ
        virtual array_1d solve(array_2d &M, array_1d &G) = 0;
    };
    
    class LASIterativeSolver : public LASSolver
    {
    public:
    
        // Конструктор
        // _eps -- Требуемая точность
        // _inum -- Максимальное число итераций
        LASIterativeSolver(double _eps = 1e-14, bool _prep = false)
            : eps(_eps), prepForIter(_prep), maxIterNum(2048) {}
        
        // Задает максимальное число итераций
        void setMaxIterNum(size_t inum)
        {
            maxIterNum = inum;
        }
        
        // Включает/отключает подготовку к итерациям
        void setPrepForIter(bool pf)
        {
            prepForIter = pf;
        }
        
        // Возвращает число выполненных итераций
        size_t getPerfIterNum() const
        {
            return iterIndex;
        }
        
    protected:
        
        // Требуемая точность
        double eps;
        
        // Если задано true, то производится подготовка
        // СЛАУ для обеспечения условий сходимости итераций
        bool prepForIter;
        
        // Максимальное число итераций
        size_t maxIterNum;
        
        // Номер текущей итерации
        size_t iterIndex;
    };

    /*
     * Метод Гаусса
     */
    class GaussSolver : public LASSolver
    {
    public:
        
        array_1d solve(array_2d &M, array_1d &G);
    };
    
    /*
     * Метод последовательных приближений
     */
    class MPPSolver : public LASIterativeSolver
    {
    public:
        
        // Конструктор по умолчанию
        MPPSolver() {}
        
        // Конструктор:
        // Если prep == true -- выполняется подготока к итерациям
        MPPSolver(bool prep)
        {
            setPrepForIter(prep);
        }
        
        array_1d solve(array_2d &M, array_1d &G);
    };
    
    /*
     * Метод наискорейшего спуска
     */
    class FGDSolver : public LASIterativeSolver
    {
    public:
        
        // Конструктор по умолчанию
        FGDSolver() {}
        
        // Конструктор:
        // Если prep == true -- выполняется подготока к итерациям
        FGDSolver(bool prep)
        {
            setPrepForIter(prep);
        }
        
        array_1d solve(array_2d &M, array_1d &G);
    };
    
} // namespace LA

#endif