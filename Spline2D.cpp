#include "Spline2D.h"
#include "GaussNodes.h"

double LagrangePoly2D::operator()(double x, double y) const
{
    double res = 0.0;
    
    for(size_t i = 0; i < xnodes.size(); ++i)
    {
        for(size_t j = 0; j < ynodes.size(); ++j)
        {
            res += values[i][j]*fundPoly(xnodes, x, i)*fundPoly(ynodes, y, j);
        }
    }
    
    return res;
}

double LagrangePoly2D::fundPoly(const array_1d &g, double x, size_t k) const
{
    double R = 1.0;
    
    for(size_t i = 0; i < g.size(); i++)
    {
        if(i == k) continue;
        R *= (x - g[i])/(g[k] - g[i]);
    }
    
    return R;
}

Spline2D::~Spline2D() 
{
    for(size_t i = 0; i < polys.rows(); ++i)
    {
        for(size_t j = 0; j < polys.cols(); ++j)
        {
            delete polys[i][j];
        }
    }
}

const Spline2D &Spline2D::operator=(const Spline2D &spl)
{
    xgrid = spl.xgrid;
    ygrid = spl.ygrid;
    
    const array_2p &srcp = spl.polys;
    
    polys = array_2p(srcp.rows(), srcp.cols());
    
    for(size_t i = 0; i < polys.rows(); ++i)
    {
        for(size_t j = 0; j < polys.cols(); ++j)
        {
            polys[i][j] = new LagrangePoly2D(*srcp[i][j]);
        }
    }
    
    return *this;
}

double Spline2D::operator()(double x, double y) const
{
    size_t i, j;
    
    for(i = 0; i < xgrid.size()-1; ++i)
    {
        if((x >= xgrid[i]) && (x <= xgrid[i+1]))
        {
            break;
        }
    }
    
    for(j = 0; j < ygrid.size()-1; ++j)
    {
        if((y >= ygrid[j]) && (y <= ygrid[j+1]))
        {
            break;
        }
    }
    
    return (*polys[i][j])(x, y);
}

// Инициализация сетки Гаусса для полинома poly
void setupGaussGrid(LagrangePoly2D &poly, const double xval[2], const double yval[2], GaussGdType type)
{
    size_t gsize = type;
    array_1d xnodes(gsize);
    array_1d ynodes(gsize);
    
    GaussNodes gNodes(type);
    const double *gnodes = gNodes.getNodes();
    
    for(size_t i = 0; i < gsize; ++i)
    {
        xnodes[i] = (xval[0] + xval[1])/2 + (xval[1] - xval[0])/2*gnodes[i];
        ynodes[i] = (yval[0] + yval[1])/2 + (yval[1] - yval[0])/2*gnodes[i];
    }
    
    poly.setXNodes(xnodes);
    poly.setYNodes(ynodes);
}

// Создает неициализированный экземпляр класса Spline2D с сеткой Гаусса
Spline2D createGaussSpline2D(const double xval[2], const double yval[2], size_t xnum, size_t ynum, GaussGdType gtype)
{
    array_1d xnodes(xnum);
    array_1d ynodes(ynum);
    
    double xstep = (xval[1] - xval[0])/(xnum - 1);
    double ystep = (yval[1] - yval[0])/(ynum - 1);
    
    for(size_t i = 0; i < xnum; ++i)
    {
        xnodes[i] = xval[0] + i*xstep;
    }
    
    for(size_t j = 0; j < ynum; ++j)
    {
        ynodes[j] = yval[0] + j*ystep;
    }
    
    Spline2D::array_2p polyArray(xnum-1, ynum-1);
    
    for(size_t i = 0; i < xnum-1; ++i)
    {
        double xr[] = {xnodes[i], xnodes[i+1]};
        for(size_t j = 0; j < ynum-1; ++j)
        {
            double yr[] = {ynodes[j], ynodes[j+1]};
            
            LagrangePoly2D *cpoly = new LagrangePoly2D();
        
            setupGaussGrid(*cpoly, xr, yr, gtype);
            polyArray[i][j] =  cpoly;
        }
    }
    
    return Spline2D(xnodes, ynodes, polyArray);
}