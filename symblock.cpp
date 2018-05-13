#include "symblock.h"

void symblock::setElement(int row, int col, double value){
    matrix(row, col) = value;
}

double symblock::getElement(int row, int col){
    return matrix(row, col);
}

void symblock::setColMap(int index, int firstInd, int secondInd){
    bool found = false;
    for (unsigned i = 0; i < ColMap.size(); i++){
        channelindexpair IndexPair = ColMap.at(i);
        if (IndexPair.first() == firstInd && IndexPair.second() == secondInd && IndexPair.index() == index){
            found = true;
        }
    }
    if (found != true) {
        ColMap.emplace_back(channelindexpair());
        ColMap.back().setMap(index, firstInd, secondInd);
    }
}

void symblock::setRowMap(int index, int firstInd, int secondInd){
    bool found = false;
    for (unsigned i = 0; i < RowMap.size(); i++){
        channelindexpair IndexPair = RowMap.at(i);
        if (IndexPair.first() == firstInd && IndexPair.second() == secondInd && IndexPair.index() == index){
            found = true;
        }
    }
    if (found != true) {
        RowMap.emplace_back(channelindexpair());
        RowMap.back().setMap(index, firstInd, secondInd);
    }
}

void symblock::setColMap1x3(int index, int firstInd){
    ColMap1x3.emplace_back(channelindexpair());
    ColMap1x3.back().setMapOne(index, firstInd);
}

void symblock::setRowMap1x3(int index, int firstInd){
    RowMap1x3.emplace_back(channelindexpair());
    RowMap1x3.back().setMapOne(index, firstInd);
}

void symblock::setColMap3x1(int index, int firstInd, int secondInd, int thirdInd){
    ColMap3x1.emplace_back(channelindexpair());
    ColMap3x1.back().setMapThree(index, firstInd, secondInd, thirdInd);
}

void symblock::setRowMap3x1(int index, int firstInd, int secondInd, int thirdInd){
    RowMap3x1.emplace_back(channelindexpair());
    RowMap3x1.back().setMapThree(index, firstInd, secondInd, thirdInd);
}

channelindexpair symblock::getRowMap(int index){
    for (unsigned i = 0; i < RowMap.size(); i++){
        channelindexpair IndexPair = RowMap.at(i);
        if (IndexPair.index() == index){
            return IndexPair;
        }
    }
}

channelindexpair symblock::getRowMap1x3(int index){
    for (unsigned i = 0; i < RowMap1x3.size(); i++){
        channelindexpair IndexPair = RowMap1x3.at(i);
        if (IndexPair.index() == index){
            return IndexPair;
        }
    }
}

channelindexpair symblock::getRowMap3x1(int index){
    for (unsigned i = 0; i < RowMap3x1.size(); i++){
        channelindexpair IndexPair = RowMap3x1.at(i);
        if (IndexPair.index() == index){
            return IndexPair;
        }
    }
}

channelindexpair symblock::getColMap(int index){
    for (unsigned i = 0; i < ColMap.size(); i++){
        channelindexpair IndexPair = ColMap.at(i);
        if (IndexPair.index() == index){
            return IndexPair;
        }
    }
}

channelindexpair symblock::getColMap1x3(int index){
    for (unsigned i = 0; i < ColMap1x3.size(); i++){
        channelindexpair IndexPair = ColMap1x3.at(i);
        if (IndexPair.index() == index){
            return IndexPair;
        }
    }
}

channelindexpair symblock::getColMap3x1(int index){
    for (unsigned i = 0; i < ColMap3x1.size(); i++){
        channelindexpair IndexPair = ColMap3x1.at(i);
        if (IndexPair.index() == index){
            return IndexPair;
        }
    }
}
 // inverse mapping can be implemented in a smart way
int symblock::getRowMapInverse(int indexOne, int indexTwo){
    for (int i = 0; i < RowMap.size(); i++){
        channelindexpair IndexPair = RowMap.at(i);
        if (IndexPair.first() == indexOne && IndexPair.second() == indexTwo){
            return IndexPair.index();
        }
    }
}

int symblock::getColMapInverse(int indexOne, int indexTwo){
    for (int i = 0; i < ColMap.size(); i++){
        channelindexpair IndexPair = ColMap.at(i);
        if (IndexPair.first() == indexOne && IndexPair.second() == indexTwo){
            return IndexPair.index();
        }
    }
}

int symblock::getRowMapInverse1x3(int indexOne){
    for (int i = 0; i < RowMap1x3.size(); i++){
        channelindexpair IndexPair = RowMap1x3.at(i);
        if (IndexPair.first() == indexOne){
            return IndexPair.index();
        }
    }
}

int symblock::getColMapInverse3x1(int indexOne, int indexTwo, int indexThree){
    for (int i = 0; i < ColMap3x1.size(); i++){
        channelindexpair IndexPair = ColMap3x1.at(i);
        if (IndexPair.first() == indexOne && IndexPair.second() == indexTwo && IndexPair.third() == indexThree){
            return IndexPair.index();
        }
    }
}

int symblock::getRowMapInverse3x1(int indexOne, int indexTwo, int indexThree){
    for (int i = 0; i < RowMap3x1.size(); i++){
        channelindexpair IndexPair = RowMap3x1.at(i);
        if (IndexPair.first() == indexOne && IndexPair.second() == indexTwo && IndexPair.third() == indexThree){
            return IndexPair.index();
        }
    }
}

int symblock::getColMapInverse1x3(int indexOne){
    for (int i = 0; i < ColMap1x3.size(); i++){
        channelindexpair IndexPair = ColMap1x3.at(i);
        if (IndexPair.first() == indexOne){
            return IndexPair.index();
        }
    }
}

void symblock::setZeros(){
    matrix.zeros();
}

void symblock::setMatrix(arma::mat X){
    matrix = X;
}
