#ifndef SYMBLOCK_H
#define SYMBLOCK_H

#include <armadillo>
#include "channelindexpair.h"

class symblock
{
public:
    symblock(int row, int col):
        matrix(row, col){matrix.zeros();}

    arma::mat getMatrix(){return matrix;}
    void setMatrix(arma::mat X);
    void setZeros();
    void setElement(int, int, double);
    double getElement(int, int);

    // set mapping; ColMap - 2x2 - particles/holes: for ex. - (x + y) = (z + w)
    //              ColMap3x1 - for Q3 or Q4 diagrams: for ex. -(x) = (z + w - y)
    // first int - index to map, other ints - for creation of object channelindexpair. index is mapped to an object
    void setColMap(int, int, int);
    void setRowMap(int, int, int);
    void setRowMap3x1(int, int, int, int);
    void setColMap1x3(int, int);
    void setRowMap1x3(int, int);
    void setColMap3x1(int, int, int, int);

    // get reverse mapping; returns index corresponding to a cetain channelindexpair.
    int getRowMapInverse(int, int);
    int getColMapInverse(int, int);
    int getRowMapInverse3x1(int, int, int);
    int getColMapInverse3x1(int, int, int);
    int getRowMapInverse1x3(int);
    int getColMapInverse1x3(int);

    // get direct mapping
    channelindexpair getColMap(int);
    channelindexpair getRowMap(int);
    channelindexpair getColMap1x3(int);
    channelindexpair getRowMap1x3(int);
    channelindexpair getColMap3x1(int);
    channelindexpair getRowMap3x1(int);

private:
    arma::mat matrix;
    std::vector<channelindexpair> RowMap;  // 2 index default
    std::vector<channelindexpair> ColMap;  // 2 index default
    std::vector<channelindexpair> RowMap1x3; // one particle state maps to row index
    std::vector<channelindexpair> ColMap3x1; // three particle states map to col index
    std::vector<channelindexpair> RowMap3x1; // three particle states map to row index
    std::vector<channelindexpair> ColMap1x3; // one particle state maps to col index
};

#endif // SYMBLOCK_H
