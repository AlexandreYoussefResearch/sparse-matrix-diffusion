#include <cassert>
#include "DiffusionProblem.hpp"




// Constructor
DiffusionProblem::DiffusionProblem(double Lt, double Lx,double Ly,std::function<double(double,double)> getD, std::function<double(double,double)> getInitial, std::function<double(double,double)> getf, std::function<double(double)> getBoundaryx0, std::function<double(double)> getBoundaryx1, std::function<double(double)> getBoundaryy0, std::function<double(double)> getBoundaryy1){
    this->Lt = Lt;
    this->Lx = Lx;
    this->Ly = Ly;
    this->getD = getD;
    this->getInitial = getInitial;
    this->getf = getf;
    this->getBoundaryx0 = getBoundaryx0;
    this->getBoundaryx1 = getBoundaryx1;
    this->getBoundaryy0 = getBoundaryy0;
    this->getBoundaryy1 = getBoundaryy1;
}



// Create the accessors
double DiffusionProblem::getLt()const{
    return this->Lt;
}

double DiffusionProblem::getLx()const{
    return this->Lx;
}

double DiffusionProblem::getLy()const{
    return this->Ly;
}

