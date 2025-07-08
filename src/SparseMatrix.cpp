#include <iostream>
#include <cassert> // for the assert function
#include "SparseMatrix.hpp"


//1. Constructor of zero sparse matrix
SparseMatrix::SparseMatrix(int numRows, int numCols){
    this->m = numRows;
    this->n = numCols;
    this->M = 0;
    this->nnz = 0;
    this->rowidx = NULL;
    this->nzval = NULL;
}

//2. Constructor of deep copy
SparseMatrix::SparseMatrix(int m, int n, int nnz, int M, int* rowidx, double* nzval){
    this->m = m;
    this->n = n;
    this->nnz = nnz;
    this->M = M;
    this->rowidx = new int[M*n];
    this->nzval = new double[M*n]; 

    for (int i =0; i<M*n; i++){
        this->nzval[i] = nzval[i];
        this->rowidx[i] = rowidx[i];
    }
}

//3. Constructor sparse matrix

SparseMatrix::SparseMatrix(int nnz, int const *ridx, int const *cidx, double const *nzval, int m, int n){
    this->nnz = nnz;
    this->m = m;
    this->n = n;

    // Calcul de M :
    int M = 0;
    int temp = 0;
    for (int i =0; i<this->n;i ++){
        for (int j=0; j<this->nnz;j++){
            if (cidx[j] == i){
                temp++;
            }
        }
        if (temp > M){
            M = temp;
        }
        temp = 0;
    }
    this->M = M;


    int *true_ridx = new int[M*n];
    double *true_nzval = new double[M*n];

    int index;
    int nnz_per_col;
    for (int i=0; i<this->n; i++){
        index = 0;
        nnz_per_col = 0;
        int *temp_ridx = new int[M];
        double *temp_val = new double[M];  

        for (int j=0; j<this->nnz;j++){
            if (cidx[j]== i){
                temp_ridx[index] = ridx[j];
                temp_val[index] = nzval[j];
                index++;
                nnz_per_col++;
            }
        }

        //Trier le tableau en fonction du row

        //Rajouter dans le tableau final par ordre croissant de temp_ridx
        for (int k = 0; k<nnz_per_col;k++){
            true_nzval[M*i +k] = temp_val[k];
            true_ridx [M*i +k] = temp_ridx[k];
        }
        //Ajouter les -1 restants dans rowidx : 
        for (int k = nnz_per_col; k<M; k++){
            true_ridx[M*i +k] = -1;
        }
      
        delete[] temp_ridx;
        delete[] temp_val;
    }

    this->rowidx = true_ridx;
    this->nzval = true_nzval;
}

// Destructor
SparseMatrix ::~SparseMatrix(){
    delete[] rowidx;
    delete[] nzval;
}

//4. Accessor Getsize
int SparseMatrix::GetSize(int val)const{
    assert(val == 1 || val == 2);
    if (val == 1){
        return this->m;
    }
    else if (val == 2){
        return this->n;
    }
    return 0;
}

//5. Operator = with another sparse matrix
SparseMatrix& SparseMatrix::operator=(const SparseMatrix& otherSparseMatrix){

    assert(this->m == otherSparseMatrix.m && this->n == otherSparseMatrix.n);

    int *true_ridx = new int[otherSparseMatrix.M*otherSparseMatrix.n];
    double *true_nzval = new double[otherSparseMatrix.M*otherSparseMatrix.n];

    this->M = otherSparseMatrix.M;
    this->nnz = otherSparseMatrix.nnz;

    for (int i=0; i< otherSparseMatrix.M*otherSparseMatrix.n;i++){
        true_ridx[i] = otherSparseMatrix.rowidx[i];
        true_nzval[i] = otherSparseMatrix.nzval[i];
    }

    delete[] this->rowidx;
    delete[] this->nzval;

    this->rowidx = true_ridx;
    this->nzval = true_nzval;

    return * this;
}

//6. +, -, and scalar multiplication

SparseMatrix SparseMatrix::operator+() const{
    SparseMatrix newMatrix(this->m, this->n, this->nnz, this->M, this->rowidx,this->nzval);
    return newMatrix;
}

SparseMatrix SparseMatrix::operator-() const {
    
    double  *newval = new double[this->M*this->n];
    int *newid = new int[this->M*this->n];
    for (int i=0;i<this->M*this->n;i++){
        newval[i] = - this->nzval[i];
        newid[i] = this->rowidx[i];
    }
    SparseMatrix newMatrix(this->m, this->n, this->nnz, this->M, newid,newval);
    delete[] newval;
    delete[] newid;
    return newMatrix;
}

SparseMatrix SparseMatrix::operator*(double a)const{
    double  *newval = new double[this->M*this->n];
    int *newid = new int[this->M*this->n];
   for (int i=0; i< this->M*this->n; i++){
           newval[i] = a *this->nzval[i];
           newid[i] = this->rowidx[i];
       
   }
    SparseMatrix newMatrix(this->m, this->n, this->nnz, this->M, newid,newval);
    delete[] newval;
    delete[] newid;
    return newMatrix;
}


//7. Matrix-vector and matrix vector multiplication

Vector operator *(const SparseMatrix& m,const Vector& v){
    assert(m.n == v.GetSize());

    Vector final_vect(m.m);

    for (int i=0; i<m.n; i++){

        for (int j=0; j<m.M; j++){
            int j_bis = j+m.M*i;
            int c = m.rowidx[j_bis];
            if ((c>=0) && (c<m.n)){
                final_vect(c) += m.nzval[j_bis]*v.Read(i);
            }
            else{
                break;
            } 
        }
    }
    return final_vect;
}

