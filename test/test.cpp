#include "FE2DSolver.h"

#include <cstdio>
#include <cmath>
#include <omp.h>

// Точное решение
double phi(double x1, double x2)
{
    return (2.0 + cos(4.0*x1))/(2.0 - sin(4.0*x2));
}

// Ядро уравнения
double ker(double x1, double x2, double s1, double s2)
{
    return ( sin(s1)*cos(4.0*x1)*(2.0 - sin(4.0*s2)) )/( 2.0 - sin(4.0*x2) );
}

// Правая часть
double f(double x1, double x2)
{
    return 2.0/(2.0 - sin(4.0*x2));
}

// Вычисляет максимальную погрешность
double maxError(const Spline2D &fnca, function_2d fncb, double xr[2], double yr[2])
{
    double maxerr = 0.0;
    
    size_t n = 20;
    
    double xs = (xr[1] - xr[0])/n;
    double ys = (yr[1] - yr[0])/n;
    
    double x = xr[0];
    for(size_t i = 0; i < n; ++i, x += xs)
    {
        double y = yr[0];
        for(size_t j = 0; j < n; ++j, y += ys)
        {
            double err = std::fabs(fnca(x, y) - fncb(x, y));
            
            if(maxerr < err)
            {
                maxerr = err;
            }
        }
    }
    
    return maxerr;
}

int main(int argc, char **argv)
{
    // Исходные данные задачи
    FE2D::EqData data;
    // Параметры сплайна, аппроксимирующего решения
    FE2D::SCSolver::SplineOpts sopts;
    
    // Пределы интегрирования
    double xr[] = {0.0, M_PI};
    double yr[] = {0.0, M_PI};
    
    // Задаем входные данные
    data.K  = ker;
    data.f  = f;
    data.lm = 15.0/(58.0*M_PI);
    data.a1 = xr[0];
    data.a2 = xr[1];
    data.b1 = yr[0];
    data.b2 = yr[1];
    
    // Параметры сплайна
    sopts.xnum = 16; // Число полиномов по x1
    sopts.ynum = 16; // Число полиномов по x2
    sopts.gtype = G8;
    
    // Метод решения СЛАУ
    //LA::LASSolver *laSolver = new LA::GaussSolver();
    //LA::LASIterativeSolver *laSolver = new LA::MPPSolver(false);
    LA::LASIterativeSolver *laSolver = new LA::FGDSolver(false);
    
    // Максимальное число итераций
    laSolver->setMaxIterNum(32);
    
    // Метод решения ИУФ 2-го рода
    FE2D::FESolver *feSolver = new FE2D::SCSolver(sopts, laSolver);
    
    double time;
    
    // Засекаем время
    time = omp_get_wtime();
    // Запускаем метод, реализующий приближенное решение ИУФ
    Spline2D spl = feSolver->solve(data);
    time = omp_get_wtime() - time;
    
    // Вычисляем максимальную погрешность
    double maxErr = maxError(spl, phi, xr, yr);
    
    // Вывод результатов
    printf("Погрешность : %e\n", maxErr);
    printf("Выполнено итераций: %lu\n", laSolver->getPerfIterNum());
    printf("Времени затрачено : %lf\n", time);
    
    delete feSolver;
    delete laSolver;
    
    return 0;
}
