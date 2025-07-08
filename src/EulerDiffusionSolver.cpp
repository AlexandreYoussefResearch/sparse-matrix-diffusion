#include <cmath>
#include <cstdio>
#include <cassert>
#include <vector>
#include <iostream>
#include "globals.hpp"
#include "EulerDiffusionSolver.hpp"
#include "SparseMatrix.hpp"
#include "DiffusionProblem.hpp"

//useful macro not mandatory to use it
#define I(a, b) (b)*(nx+1)+(a)


// Constructor :
EulerDiffusionSolver :: EulerDiffusionSolver(int nt, int nx, int ny, const DiffusionProblem &p) : A((nx+1)*(ny+1),(nx+1)*(ny+1)),b((nx+1)*(ny+1)),p(p){
    this->nt = nt;
    this->nx = nx;
    this->ny = ny;

    this->dt = (p.getLt())/nt;
    this->dx = (p.getLx())/nx;
    this->dy = (p.getLy())/ny;

    this->alphax = dt/(dx*dx);
    this->alphay = dt/(dy*dy);

    this->A = computeAb(this->b);
}

    // Matrix vector product :
SparseMatrix EulerDiffusionSolver::computeAb(Vector &b){
    int m = (nx+1)*(ny+1);
    int n = (nx+1)*(ny+1);

    int *ridx = new int[m*n];
    int *cidx = new int[m*n];
    double *nzval = new double[m*n];
    //fill in the values of the matrix A:
    int nnz = 0;
    for (int i=0; i<=nx;i++){
        for (int j=0; j<=ny;j++){
            // Element x,y
            ridx[nnz] = I(i,j);
            cidx[nnz] = I(i,j);
            nzval[nnz] = -2*alphax*p.getD(i*dx,j*dy)-2*alphay*p.getD(i*dx,j*dy) +1;
            nnz++;

            // Element x+1,y
            if (i != nx){
                ridx[nnz] = I(i,j);
                cidx[nnz] = I(i+1,j);
                nzval[nnz] = alphax*p.getD(i*dx,j*dy);
                nnz++;
            }

            // Element x-1,y
            if (i != 0){
                ridx[nnz] = I(i,j);
                cidx[nnz] = I(i-1,j);
                nzval[nnz] = alphax*p.getD(i*dx,j*dy);
                nnz++;
            }


            // Element x y+1
            if (j != ny){
                ridx[nnz] = I(i,j);
                cidx[nnz] = I(i,j+1);
                nzval[nnz] = alphay*p.getD(i*dx,j*dy);
                nnz++;
            }

            // Element x y-1
            if (j != 0){
                ridx[nnz] = I(i,j);
                cidx[nnz] = I(i,j-1);
                nzval[nnz] = alphay*p.getD(i*dx,j*dy);
                nnz++;
            }

            b(I(i,j)) = dt*p.getf(i*dx,j*dy);
        }
    }


    int *ridx2 = new int[nnz];
    int *cidx2 = new int[nnz];
    double *nzval2 = new double[nnz];

    for (int i = 0; i<nnz ; i++){
        ridx2[i] = ridx[i];
        cidx2[i] = cidx[i];
        nzval2[i] = nzval[i];
    }

    delete[] ridx;
    delete[] cidx;
    delete[] nzval;

    //Create matrix A:
    SparseMatrix A(nnz,ridx2,cidx2,nzval2,m,n);
    
    delete[] ridx2;
    delete[] cidx2;
    delete[] nzval2;
    
    return A;
}

void EulerDiffusionSolver::solve(Vector **sol){
    for (int r=0; r<=this->nx; r++){
        for (int c=0; c<=this->ny;c++){
            (*sol[0])[I(r,c)] = this->p.getInitial(r*dx,c*dy);
        }
    }

    // Iterate over time until Lt :
    for (int i = 1;i<= this->nt;i++){
        *sol[i] = (this->A* (*sol[i-1])) + this->b;

        //Boundary conditions:
        for (int j=0; j<=nx;j++){
            (*sol[i])[I(j,0)] = this->p.getBoundaryy0(j*dx);
            (*sol[i])[I(j,ny)] = this->p.getBoundaryy1(j*dx);
        }

        for (int j=0; j<=ny;j++){
            (*sol[i])[I(0,j)] = this->p.getBoundaryx0(j*dy);

            (*sol[i])[I(-1,j+1)] = this->p.getBoundaryx1(j*dy);
        }
    
    }
}



