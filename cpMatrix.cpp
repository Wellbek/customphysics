#include "cpMatrix.h"

using std::cout;
using std::vector;
using std::string;
using std::endl;
using std::dec;

void cpMatrix::print(string name){
    cout << "Matrix " << name << ":" << endl;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            cout << dec << (*this)(i,j) << " , ";
        }
        cout << endl;
    }
    cout << endl;
}

cpMatrix::cpMatrix(int _n, int _m){
    n = _n;
    m = _m;
    for(int i = 0; i < n*m; i++){
        values.push_back(0);
    }
}

/*Returns the transposed Matrix*/
cpMatrix cpMatrix::transpose(){
    cpMatrix res(m,n);
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            res(i,j) = (*this)(j,i);
        }
    }
    return res;
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

/*turns 3x1 cpMatrix to btVector3
if the format is not 3x1 it will return (0,0,0)*/
btVector3 cpMatrix::toBtVector3(){
    btVector3 res = {0,0,0};
    if(n != 3 || m != 1) return res;
    for(int i = 0; i < 3; i++){
        res[i] = (*this)(i,0);
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

/*initializes a block inside a matrix, can start anywhere not just in steps of 3
returns -1 if failed*/
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

btVector3 cpMatrix::getBtVector3(int start_row, int col){
    btVector3 res = {0,0,0};
    if(start_row+3 > n || col >= m) return res;
    res.setX((*this)(start_row,col));
    res.setY((*this)(start_row+1,col));
    res.setZ((*this)(start_row+2,col));
    return res;
}

/*sets a 3x1 block inside a matrix given a btVector3 and starting row and column
returns -1 if failed*/
int cpMatrix::setWithBtVector3(btVector3 vec, int start_row, int col){
    if(start_row+3 > n || col >= m) return -1;
    (*this)(start_row, col) = vec.getX();
    (*this)(start_row+1, col) = vec.getY();
    (*this)(start_row+2, col) = vec.getZ();
    return 0;
}

double& cpMatrix::operator()(int i, int j)
{
    if(i >= n || j >= m){ 
        cout << "Matrix operation failed; Out of bounds exception: i,n: " << 
                    dec << i << " , " << dec << n << ";  j,m:" << dec << j << " , " << dec << m << endl; 
    }
    return values[i * m + j];
}

double cpMatrix::operator()(int i, int j) const
{
    if(i > n || j > m){ 
        cout << "Matrix operation failed; Out of bounds exception: i,n: " << 
                    dec << i << " , " << dec << n << ";  j,m:" << dec << j << " , " << dec << m << endl; 
        return 0;
    }
    return values[i * m + j];
}

cpMatrix& cpMatrix::operator=(const btVector3& other){
    n = 3;
    m = 1;
    values = {other.getX(), other.getY(), other.getZ()};
    return *this;
}

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
    //cout << dec << m << " " << other.n;
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
            res(i,j) = (*this)(i,j) * scalar;
        }
    }
    return res;
}

/*if you want to test this class, go in your terminal to /modules and compile it with:
g++ -o cpMatrix customphysics/cpMatrix.cpp -I ../thirdparty/bullet*/
/*int main(){
    btMatrix3x3 I = btMatrix3x3::getIdentity();

    cpMatrix m1(3,3);
    m1(0,0) = 1;
    m1(1,2) = 12;
    m1(2,0) = 3;
    m1.print("m1");

    cpMatrix m2(3,12);

    btMatrix3x3 arr[4] = {I,I*2,I*3,I*4};

    m2.initializeWithBtMatrixArray(arr, 4);
    m2.print("m2");

    cpMatrix m3 = m1*m2;
    m3.print("m3");

    cpMatrix m4 = m2+m3;
    m4.print("m4");

    m4 -= m2;
    m4.print("m5");

    return 0;
}*/