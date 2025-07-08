#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "EulerDiffusionSolver.hpp"
#include "DiffusionProblem.hpp"
#include "globals.hpp"

using namespace std;

 #define I(a, b) (b)*(nx+1)+(a)


void write_tensor(const char fname[], Vector** A, int nt, int nx, int ny) {
    FILE* f = fopen(fname, "w");
    if (f == NULL) {
        perror("fopen");
    }
    int err = fprintf(f, "%d %d %d\n", nt+1, nx+1, ny+1);
    if (err < 0) fprintf(stderr, "Error when writing. Do you have enough space left ?");
    for (int k = 1; k <= nt+1; k++)
        for (int i = 1; i <= nx+1; i++)
            for (int j = 1; j <= ny+1; j++) {
                    err = fprintf(f, "%.16lf\n", (*A[k-1]).Read((j-1)*(nx+1)+i-1));
                    if (err < 0) fprintf(stderr, "Error when writing. Do you have enough space left ?");
                }
        
    err = fclose(f);
    if (err != 0) perror("fclose");
}

int main() {
    int nt = 2000, nx =80, ny = 80;
    int N = (nx+1) * (ny+1);
    Vector *sol[nt+1];
    for (int i = 0; i <= nt; i++) {
        sol[i] = new Vector(N);
    }

    // Your problem here (we use lambda function but not mandatory)
    std::function<double(double)> constanty0 = [](double) {return 283.15;};
    std::function<double(double)> constanty1 = [](double) {return 293.15;};
    std::function<double(double)> constantx1 = [](double x) {return 0.523956*exp(x)+282.626;};
    //DiffusionProblem(double Lt, double Lx, double Ly, getD, getInitial, getf, getBoundaryx0, getBoundaryx1, getBoundaryy0, getBoundaryy1);
    DiffusionProblem p(.05, 5., 3., [](double x,double y) {return 4.47*pow(10,-3)*x + 20.36*pow(10,-2);}, [](double x,double y) {return 283.15 + (10./3.)*y;}, [](double x, double y) {return -x +5;}, constantx1, constantx1, constanty0, constanty1);
    EulerDiffusionSolver solver(nt, nx, ny, p);
    solver.solve(sol);
    

    
     //FOR THE RATE OF CONVERGENCE
    int nx2 = 80;
    int ny2 = 80;
    int nt2 = 10000;
    int N2 = (nx2+1) * (ny2+1);
    Vector *sol_true[nt2+1];
    for (int i = 0; i <= nt2; i++) {
        sol_true[i] = new Vector(N2);
    }

    EulerDiffusionSolver solver2(nt2, nx2, ny2, p);
    solver2.solve(sol_true);
    #define I2(a, b) (b)*(nx2+1)+(a)

    double error = 0;
    for (int i=0;i<=nx;i++){
        for (int j=0;j<=ny;j++){
        error += ((abs((*sol_true[nt2])[I2(i,j)] - (*sol[nt])[I(i,j)])))/((*sol_true[nt2])[I2(i,j)]);
        }
    }
    error = 100*(error/((nx+1)*(ny+1)));
    printf("Error = %.16lf \n", error);
    
    for (int i = 0; i <= nt2; i++) {
        delete sol_true[i];
    }
    
   
    //write_tensor("sol.txt", sol, nt, nx, ny);
    for (int i = 0; i <= nt; i++) {
        delete sol[i];
    }
    return EXIT_SUCCESS;
}
