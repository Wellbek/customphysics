#ifndef CPMATRIX_H
#define CPMATRIX_H

#include <vector>
#include <iostream>
#include <string>

#include <LinearMath/btMatrix3x3.h>

using std::vector;
using std::string;

class cpMatrix{
    private:
        int n, m;
        vector<double> values = {};
    public:
        cpMatrix(int,int);

        int getN(){return n;};
        int getM(){return m;};

        int initializeWithBtMatrixArray(btMatrix3x3[],int);
        int initializeBlock(btMatrix3x3,int,int);
        cpMatrix transpose();

        btMatrix3x3 toBtMatrix3x3();
        btVector3 toBtVector3();

        btVector3 getBtVector3(int, int);
        int setWithBtVector3(btVector3, int, int);

        double& operator()(int i, int j);
        double operator()(int i, int j) const;

        cpMatrix& operator=(const btVector3&);
        cpMatrix operator+(const cpMatrix&);
        cpMatrix operator-(const cpMatrix&);
        cpMatrix operator*(const cpMatrix&);
        cpMatrix operator*(const double&);
        cpMatrix operator+=(const cpMatrix&);
        cpMatrix operator-=(const cpMatrix&);
        cpMatrix operator*=(const cpMatrix&);

        void print(string);

};

#endif