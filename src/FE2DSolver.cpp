#include "FE2DSolver.h"

namespace FE2D
{
    Spline2D SCSolver::solve(const EqData &data)
    {
        size_t nx = opts.xnum;
        size_t ny = opts.ynum;
        
        GaussNodes gnodes(opts.gtype);
        const double *gw = gnodes.getWeights();
    
        double xr[] = {data.a1, data.a2};
        double yr[] = {data.b1, data.b2};
    
        Spline2D spline = createGaussSpline2D(xr, yr, nx + 1, ny + 1, opts.gtype);
        Spline2D::array_2p &polynoms = spline.getPolynoms();
        array_1d &splXGrid = spline.getXGrid();
        array_1d &splYGrid = spline.getYGrid();
    
        size_t r = polynoms[0][0]->getXNodes().size();
        size_t sys_sz = nx*ny*r*r; // Размер системы
        array_1d sys_fvector(sys_sz); // Вектор правой части СЛАУ
        array_2d sys_matrix(sys_sz, sys_sz); // Матрица СЛАУ
        
        // Заполняем вектор правой части
        for(size_t i = 0; i < nx; ++i)
        {
        for(size_t j = 0; j < ny; ++j)
        {
            LagrangePoly2D *cpoly = polynoms[i][j];
            const array_1d &xgrid = cpoly->getXNodes();
            const array_1d &ygrid = cpoly->getYNodes();
            
            for(size_t k = 0; k < r; ++k)
            {
            for(size_t p = 0; p < r; ++p)
            {
                size_t index = i*ny*r*r + j*r*r + k*r + p;
                sys_fvector[index] = data.f(xgrid[k], ygrid[p]);
            }
            }
        }
        }
    
        // Заполняем матрицу системы
        #pragma omp parallel shared(sys_matrix, polynoms, data)
        {
        
        #pragma omp for collapse(2) // Распараллеливаем внешние циклы for по индексам i, j
        for(size_t i = 0; i < nx; ++i){ for(size_t j = 0; j < ny; ++j) // i j
        {
            LagrangePoly2D *cpolyij = polynoms[i][j];
            const array_1d &xgrid = cpolyij->getXNodes();
            const array_1d &ygrid = cpolyij->getYNodes();
            
            for(size_t q = 0; q < nx; ++q){ for(size_t s = 0; s < ny; ++s) // q s
            {
                LagrangePoly2D *cpolyqs = polynoms[q][s];
                const array_1d &xgridqs = cpolyqs->getXNodes();
                const array_1d &ygridqs = cpolyqs->getYNodes();
                
                double c2 = 0.5*(splXGrid[q+1] - splXGrid[q]);
                double c4 = 0.5*(splYGrid[s+1] - splYGrid[s]);
                
                for(size_t k = 0; k < r; ++k){ for(size_t p = 0; p < r; ++p) // k p
                {
                    double xik1 = xgrid[k];
                    double xjp2 = ygrid[p];
                    
                    for(size_t t = 0; t < r; ++t){ for(size_t u = 0; u < r; ++u) // t u
                    {
                        size_t index1 = i*ny*r*r + j*r*r + k*r + p;
                        size_t index2 = q*ny*r*r + s*r*r + t*r + u;
                        
                        double xqt1 = xgridqs[t];
                        double xsu2 = ygridqs[u];
                    
                        // Вычисление текущего элемента матрицы СЛАУ
                        double qint = data.lm*c2*c4*gw[t]*gw[u]*data.K(xik1, xjp2, xqt1, xsu2);
                        sys_matrix[index1][index2] = ( (index1 == index2) ? (1.0 - qint) : (-qint) );
                    }}
                }}
            }}
        }}
    
        }
    
        // Решаем систему
        array_1d sys_sln = laSolver->solve(sys_matrix, sys_fvector);
        
        // Достаем коэффициенты
        for(size_t i = 0; i < nx; ++i)
        {
        for(size_t j = 0; j < ny; ++j)
        {
            LagrangePoly2D *cpoly = polynoms[i][j];
            array_2d cpolyCoeffs(r, r);
            
            for(size_t k = 0; k < r; ++k)
            {
            for(size_t p = 0; p < r; ++p)
            {
                size_t index = i*ny*r*r + j*r*r + k*r + p;
                cpolyCoeffs[k][p] = sys_sln[index];
            }
            }
            
            // Записываем коэффициенты в текущий полином
            cpoly->setValues(cpolyCoeffs);
        }
        }
    
        return spline;
    }
    
} // namespace FE2D