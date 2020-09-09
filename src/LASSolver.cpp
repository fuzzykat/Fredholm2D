#include "LA.h"
#include "LASSolver.h"

namespace LA
{
    // Метод Гаусса. Решаем систему M*F = G
    array_1d GaussSolver::solve(array_2d &M, array_1d &G)
    {
        array_2d A(M.rows(), M.cols() + 1); // Расширенная матрица системы
        array_1d F(G.size());  // Вектор решения
        
        size_t i, j, k;
        size_t n = G.size();
        double pivot, m;
        double max;
        int max_index;
        
        {	// Формируем расширенную матрицу
            for(i = 0; i < n; i++)
                for(j = 0; j < n; j++)
                    A[i][j] = M[i][j];
            for(i = 0; i < n; i++)
                A[i][n] = G[i];
        }
	
        for(i = 0; i < n; i++) 
        {
            // Поиск максимального элемента строки
            max = fabs(A[i][i]); // Макс. Элемент
            max_index = i; // Индекс макс. элем.
            for(j = i + 1; j<n; j++)
            {
                double tmp = fabs(A[i][j]);
                if(tmp > max) 
                {
                    max = A[i][j];
                    max_index = j;
                }
            }

            //Меняем строки
            for(k = i; k < n + 1; k++)
            {
                pivot = A[max_index][k];
                A[max_index][k] = A[i][k];
                A[i][k] = pivot;
            }
    
            for(k = i+1; k < n; k++) 
            {
                m =- A[k][i]/A[i][i];
                for(j = i; j < n + 1; j++) 
                {
                    A[k][j] = A[k][j] + m*A[i][j];
                }
            }
        }
	
        for(int i = n - 1; i >= 0; i--) 
        {
            F[i] = A[i][n]/A[i][i];
            for(int k = i - 1; k >= 0; k--) 
            {
                A[k][n] -= A[k][i]*F[i];
            }
        }
        
        return F;
    }
    
    array_1d MPPSolver::solve(array_2d &B, array_1d &G)
    {
        Vector X = G;
        Vector X_p = X;
        
        if(prepForIter)
        {
            Matrix B_T = B;
            
            // Транспонируем
            mTranspose(B_T);
            // B = B^T * B
            mmMult(B_T, B);
            // G = B^T * G;
            mvMult(B_T, G);
        }
        
        // Вычисляем норму B
        double mu = mNorm(B);
        // Поправочный множитель
        double gamma = 2.0/mu;
        // B = B * (-gamma)
        cmMult(B, -gamma);
        // B = B + E
        cmAdd(B, 1.0);
        // G = G * gamma
        cvMult(G, gamma);
        
        for(size_t i = 0; i < maxIterNum; ++i)
        {
            X_p = X;
            iterIndex = i + 1;
            
            // X = B * X
            mvMult(B, X);
            // X = X + G
            vvAdd(G, X);
            
            if(vvMaxAbsDiff(X_p, X) < eps)
            {
                break;
            }
        }
        
        return X;
    }
    
    array_1d FGDSolver::solve(array_2d &B, array_1d &G)
    {
        
        if(prepForIter)
        {
            Matrix B_T = B;
            
            // Транспонируем
            mTranspose(B_T);
            // B = B^T * B
            mmMult(B_T, B);
            // G = B^T * G;
            mvMult(B_T, G);
        }
        
        // Начальное приближение
        Vector X = G;
        
        for(size_t i = 0; i < maxIterNum; ++i)
        {
            iterIndex = i + 1;
            
            // Вектор невязки
            Vector R_p = X;
            // R_p = B*X;
            mvMult(B, R_p);
            // R_p = B*X - G
            vvSub(G, R_p);
            
            if(vNorm(R_p) < eps)
            {
                break;
            }
            
            // BR_p = R_p
            Vector BR_p = R_p;
            // BR_p = B*R_p
            mvMult(B, BR_p);
            
            // Вычисляем коэффициент
            double alpha = vvProd(R_p, R_p)/vvProd(BR_p, R_p);
            
            // R_p = R_p * alpha
            cvMult(R_p, alpha);
            // X = X - R_p
            vvSub(R_p, X);
        }
        
        return X;
    }
    
} // namespace LA
