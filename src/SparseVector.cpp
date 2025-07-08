#include <iostream>
#include <cassert> // for the assert function
#include "SparseVector.hpp"
#include <cmath>
#include <stdio.h>

//1. Create empty constructor
SparseVector::SparseVector() : AbstractVector(1){
    this->nnz = 0;
    this->rowidx = NULL;
    this->nzval = NULL;
}

//2. Create destructor
 SparseVector::~ SparseVector(){
     delete[] rowidx;
     delete[] nzval;
 }

// size = vraie taille du vecteur
// nnz = nombre de valeur non nulles : longueur de rowidx et de nzval
// rowidx = index des lignes non nulles
// nzval = valeurs aux indexs non nuls

//3. Create constructor
SparseVector::SparseVector(int nnz, int const *rowidx, double const *nzval, int size) : AbstractVector(size)
{
    this->nnz = nnz;

    this->rowidx = new int[nnz];
    this->nzval = new double[nnz]; 

    for (int i=0; i<nnz; i++){
        
        this->rowidx[i] = rowidx[i];
        this->nzval[i] = nzval[i]; 
    }
    
}


//4. Assignation operator - copy vector into the vector
SparseVector& SparseVector::operator=(const SparseVector& otherVector){
    // On change size si les deux vecteurs ne sont pas de la même taille
    if (GetSize() != otherVector.GetSize()){
        mSize = otherVector.GetSize();
    }

    if (this->nnz == otherVector.nnz){
        for (int i=0; i<this->nnz ; i++){
            this->rowidx[i] = otherVector.rowidx[i];
            this->nzval[i] = otherVector.nzval[i];
        }
    }

    else{
        this->nnz = otherVector.nnz;
        delete[] this->rowidx;
        delete[] this->nzval;

        int *rowidx2 = new int[this->nnz];
        double *nzval2 = new double[this->nnz];
        this->rowidx = rowidx2;
        this->nzval = nzval2;

        for (int i=0; i<this->nnz;i++){
            this->rowidx[i] = otherVector.rowidx[i];
            this->nzval[i] = otherVector.nzval[i];
        }
    }

    return * this; // On demande de retourner l'adresse
}



//5. Operator+ to perform elementwise addition
SparseVector SparseVector::operator+(const SparseVector& v1) const{
    assert(GetSize() == v1.GetSize());
    
    int *temp_idx = new int[nnz + v1.nnz];
    double *temp_value = new double[nnz + v1.nnz];

    int final_nnz = 0;
    int counter_this = 0;
    int counter_v1 = 0;

    for (int i=0; i< (nnz + v1.nnz); i++){

        if (this->rowidx[counter_this] < v1.rowidx[counter_v1]){
            temp_idx[final_nnz] = this->rowidx[counter_this];
            temp_value[final_nnz] = this->nzval[counter_this];
            counter_this++;
            final_nnz++;
        }

        else if (this->rowidx[counter_this] > v1.rowidx[counter_v1]){
            temp_idx[final_nnz] = v1.rowidx[counter_v1];
            temp_value[final_nnz] = v1.nzval[counter_v1];

            counter_v1++;
            final_nnz++;
        }

        else if (this->rowidx[counter_this] == v1.rowidx[counter_v1]){

            temp_idx[final_nnz] = this->rowidx[counter_this];
            temp_value[final_nnz] = this->nzval[counter_this] + v1.nzval[counter_v1];
            counter_this++;
            counter_v1++;
            final_nnz++;
        }

        if (counter_this >= this->nnz){
            break;
        }

        else if (counter_v1 >= v1.nnz){
            break;
        }
    }

    if (counter_this >= this->nnz){
        for (int i=counter_v1; i< v1.nnz;i++){
            temp_idx[final_nnz] = v1.rowidx[i];
            temp_value[final_nnz] = v1.nzval[i];
            final_nnz++;
        }
    }

    else if (counter_v1 >= v1.nnz){
        for (int i=counter_this; i< this->nnz;i++){
            temp_idx[final_nnz] = this->rowidx[i];
            temp_value[final_nnz] = this->nzval[i];
            final_nnz++;
        }
    }
    

    // Creer la structure du vecteur à retourner :
    SparseVector FinalSparse(final_nnz, temp_idx,temp_value, GetSize());
    delete[] temp_idx;
    delete[] temp_value;
    return FinalSparse;
}


// Read zero based indexing 


double SparseVector::Read(int i) const{
    assert(i > -1);
    assert(i < GetSize());
    for (int j = 0; j<nnz; j++){
        if (rowidx[j] == i){
            return nzval[j];
        }
    }
    return 0;
}
