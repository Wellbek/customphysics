#include <vector>
#include <iostream>
#include <string>

#include <LinearMath/btMatrix3x3.h>

class cpMatrix{
    private:
        int n, m;
        std::vector<double> values = {};
    public:
        cpMatrix(int,int);

        int getN(){return n;};
        int getM(){return m;};

        int initializeWithBtMatrixArray(btMatrix3x3[],int);
        int initializeBlock(btMatrix3x3,int,int);

        btMatrix3x3 toBtMatrix3x3();

        double& operator()(int i, int j);
        double operator()(int i, int j) const;

        //cpMatrix& operator=(const cpMatrix&);
        cpMatrix operator+(const cpMatrix&);
        cpMatrix operator-(const cpMatrix&);
        cpMatrix operator*(const cpMatrix&);
        cpMatrix operator*(const double&);
        cpMatrix operator+=(const cpMatrix&);
        cpMatrix operator-=(const cpMatrix&);
        cpMatrix operator*=(const cpMatrix&);

        void print(std::string);

};

void cpMatrix::print(std::string name){
    std::cout << "Matrix " << name << ":" << std::endl;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            std::cout << std::dec << (*this)(i,j) << " , ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

cpMatrix::cpMatrix(int _n, int _m){
    n = _n;
    m = _m;
    for(int i = 0; i < n+m; i++){
        values.push_back(0);
    }
}

/*turns 3x3 cpMatrix to btMatrix3x3
if the format is not 3x3 it will return the 3x3 identity matrix*/
btMatrix3x3 cpMatrix::toBtMatrix3x3(){
    btMatrix3x3 res = btMatrix3x3::getIdentity();
    if(n != 3 || m != 3) return res;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            res[i][j] = (*this)(i,j);
        }
    }
    return res;
}

/*initializes values of the matrix with an array of bullet matrices,
for this to work, blocks need to be big enough to cover the entire Matrix and size needs to be the length of the array
The dimension of the Matrix needs to be 3xm or nx3 with n or m divisible by 3.
Returns -1 if failed and 1 if succeded*/
int cpMatrix::initializeWithBtMatrixArray(btMatrix3x3 blocks[], int size){
    if(n == 3){
        if(3*size != m) return -1;

        for(int k = 0; k < size; k++){
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < 3; j++){
                    (*this)(i,j+k*3) = blocks[k][i][j];
                }
            }
        }

    }else if (m == 3){
        if(3*size != n) return -1;

        for(int k = 0; k < size; k++){
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < 3; j++){
                    (*this)(i+k*3,j) = blocks[k][i][j];
                }
            }
        }

    }else{
        return -1;
    }
    
    return 0;
}

/*initializes a block inside a matrix, can start anywhere not just in steps of 3*/
int cpMatrix::initializeBlock(btMatrix3x3 block, int start_row, int start_col){
    if(start_row+3 <= n && start_col+3 <= m){
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                (*this)(start_row+i,start_col+j) = block[i][j];
            }
        }
    }else{
        return -1;
    }
    
    return 0;
}


double& cpMatrix::operator()(int i, int j)
{
    if(i >= n || j >= m){ 
        std::cout << "Matrix operation failed; Out of bounds exception: i,n: " << 
                    std::dec << i << " , " << std::dec << n << ";  j,m:" << std::dec << j << " , " << std::dec << m << std::endl; 
    }
    return values[i * m + j];
}

double cpMatrix::operator()(int i, int j) const
{
    if(i > n || j > m){ 
        std::cout << "Matrix operation failed; Out of bounds exception: i,n: " << 
                    std::dec << i << " , " << std::dec << n << ";  j,m:" << std::dec << j << " , " << std::dec << m << std::endl; 
        return 0;
    }
    return values[i * m + j];
}

/*cpMatrix& cpMatrix::operator=(const cpMatrix& other){
    n = other.n;
    m = other.m;
    values = other.values;
    return *this;
}*/

cpMatrix cpMatrix::operator+(const cpMatrix& other){
    if(n != other.n || m != other.m) return cpMatrix(0,0);
    cpMatrix res(n,m);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            res(i,j) = (*this)(i,j) + other(i,j);
        }
    }
    return res;
}

cpMatrix cpMatrix::operator+=(const cpMatrix& other){
    (*this) = (*this) + other;
    return *this;
}

cpMatrix cpMatrix::operator*=(const cpMatrix& other){
    (*this) = (*this) * other;
    return *this;
}

cpMatrix cpMatrix::operator-=(const cpMatrix& other){
    (*this) = (*this) - other;
    return *this;
}

cpMatrix cpMatrix::operator-(const cpMatrix& other){
    if(n != other.n || m != other.m) return cpMatrix(0,0);
    cpMatrix res(n,m);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            res(i,j) = (*this)(i,j) - other(i,j);
        }
    }
    return res;
}

cpMatrix cpMatrix::operator*(const cpMatrix& other){
    //std::cout << std::dec << m << " " << other.n;
    if(m != other.n) return cpMatrix(0,0);
    cpMatrix res(n,other.m);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < other.m; j++){
            for(int k = 0; k < m; k++){
                res(i,j) += (*this)(i,k) * other(k,j);
            }
        }
    }
    return res;
}

cpMatrix cpMatrix::operator*(const double& scalar){
    cpMatrix res(n,m);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            res(i,j) *= scalar;
        }
    }
    return res;
}

/*if you want to test this class, go in your terminal to /modules and compile it with:
g++ -o cpMatrix customphysics/cpMatrix.cpp -I ../thirdparty/bullet*/
int main(){
    btMatrix3x3 I = btMatrix3x3::getIdentity();

    cpMatrix m1(3,3);
    m1(0,0) = 1;
    m1(1,2) = 12;
    m1(2,0) = 3;
    m1.initializeBlock(I, 0, 0);
    m1.print("m1");

    cpMatrix m2(3,12);

    btMatrix3x3 arr[4] = {I,I*2,I*3,I*4};

    m2.initializeWithBtMatrixArray(arr, 4);
    m2.print("m2");

    //cpMatrix m3 = m1*m2;
    (m1*m2).print("m3");

    return 0;
}