#ifndef GAUSS_NODES_H_
#define GAUSS_NODES_H_

// Тип сетки Гаусса GN, N -- число узлов
enum GaussGdType
{
    G2 = 2, G3, G4, G5, G6, G7, G8
};

/*
 * Класс предоставляет иформацию по весовым коэффициентам
 * и узлам квадратурных формул Гаусса-Лежандра
 */
class GaussNodes
{
public:
    
    // Конструктор
    GaussNodes(GaussGdType type);
    
    // Возвращает указатель на массив узлов
    const double *getNodes() const
    {
        return gnodes;
    }
    
    // Возвращает указатель на массив весов
    const double *getWeights() const
    {
        return gweights;
    }
    
protected:
    
    // Порядок формулы
    GaussGdType gt;
    
    // Указатель на массив узлов
    const double *gnodes;
    
    // Указатель на массив весов
    const double *gweights;
};

#endif