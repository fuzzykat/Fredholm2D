#include "LA.h"

/*
 * Пространство имен LA (Linear Algebra)
 * Здесь реализованы операции с матрицами и векторами
 */
namespace LA
{
    void mvMult(const Matrix &A, Vector &V)
    {
        Vector W(V.size());
        
        #pragma omp parallel shared(W, A, V)
        {
        // Распараллеливаем внешний цикл for по индексу i
        #pragma omp for
        for(size_t i = 0; i < A.rows(); ++i)
        {
            W[i] = 0.0;
            for(size_t k = 0; k < A.cols(); k++)
            {
                W[i] += A[i][k]*V[k];
            }
        }
        }
        
        V = W;
    }
    
    void mmMult(const Matrix &A, Matrix &B)
    {
        Matrix C(A.rows(), B.cols());
        
        #pragma omp parallel shared(A, B, C)
        {
        // Распараллеливаем внешние циклы for по индексам i, j
        #pragma omp for collapse(2)
        for(size_t i = 0; i < A.rows(); ++i)
        {
            for(size_t j = 0; j < B.cols(); ++j)
            {
                C[i][j] = 0.0;
                for(size_t k = 0; k < A.cols(); k++)
                {
                    C[i][j] += A[i][k]*B[k][j];
                }
            }
        }
        }
        
        B = C;
    }
    
    void cmAdd(Matrix &A, double alpha)
    {
        #pragma omp parallel shared(A)
        {
        // Распараллеливаем внешний цикл for по индексу i
        #pragma omp for
        for(size_t i = 0; i < A.rows(); ++i)
        {
            A[i][i] += alpha;
        }
        }
    }
    
    void cvMult(Vector &A, double alpha)
    {
        #pragma omp parallel shared(A)
        {
        // Распараллеливаем внешний цикл for по индексу i
        #pragma omp for
        for(size_t i = 0; i < A.size(); ++i)
        {
            A[i] *= alpha;
        }
        }
    }
    
    void cmMult(Matrix &A, double alpha)
    {
        #pragma omp parallel shared(A)
        {
        // Распараллеливаем внешний цикл for по индексу i
        #pragma omp for collapse(2)
        for(size_t i = 0; i < A.rows(); ++i)
        {
            for(size_t j = 0; j < A.cols(); ++j)
            {
                A[i][j] *= alpha;
            }
        }
        }
    }
    
    void vvAdd(const Vector &A, Vector &V)
    {
        #pragma omp parallel shared(A, V)
        {
        // Распараллеливаем внешний цикл for по индексу i
        #pragma omp for 
        for(size_t i = 0; i < V.size(); ++i)
        {
            V[i] += A[i];
        }
        }
    }
    
    void vvSub(const Vector &A, Vector &V)
    {
        #pragma omp parallel shared(A, V)
        {
        // Распараллеливаем внешний цикл for по индексу i
        #pragma omp for 
        for(size_t i = 0; i < V.size(); ++i)
        {
            V[i] -= A[i];
        }
        }
    }
    
    void mTranspose(Matrix &A)
    {
        for(size_t i = 0; i < A.rows(); ++i)
        {
            for(size_t j = 0; j < A.cols(); ++j)
            {
                A[i][j] = A[j][i];
            }
        }
    }
    
    double mNorm(const Matrix &A)
    {
        double norm = 0.0;
        
        for(size_t i = 0; i < A.rows(); ++i)
        {
            for(size_t j = 0; j < A.cols(); ++j)
            {
                norm += A[i][j]*A[i][j];
            }
        }
        
        return sqrt(norm);
    }
    
    double vNorm(const Vector &A)
    {
        double norm = 0.0;
        
        for(size_t i = 0; i < A.size(); ++i)
        {
            norm += A[i]*A[i];
        }
        
        return sqrt(norm);
    }
    
    double vvProd(const Vector &A, const Vector &B)
    {
        double sclProd = 0.0;
        
        for(size_t i = 0; i < A.size(); ++i)
        {
            sclProd += A[i]*B[i];
        }
        
        return sclProd;
    }
    
} // namespace LA