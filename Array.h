#ifndef ARRAY_H_
#define ARRAY_H_
#include <cstdlib>
#include <valarray>

// Шаблон двумерного массива
template<class T> class valarray_2d : 
    public std::valarray< std::valarray<T> >
{
public:
    
    valarray_2d() {}
    
    valarray_2d(size_t rw, size_t cl)
        : std::valarray< std::valarray<T> >(rw)
    {
        for(size_t i = 0; i < rw; ++i)
            (*this)[i].resize(cl);
    }
    
    ~valarray_2d() {}
    
    size_t rows() const { return this->size(); }
    size_t cols() const { return (*this)[0].size(); }
};

// Одномерный массив
typedef std::valarray<double> array_1d;

// Двумерный массив
typedef valarray_2d<double> array_2d;

#endif