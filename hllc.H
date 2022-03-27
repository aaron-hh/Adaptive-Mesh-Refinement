#include <array>
#include <iostream>
#include <fstream>
#include <string>

#include "AMREX_Array4.H"

typedef std::array<double,4> Arrayofdouble;


class numerical_method
{
private:
    Arrayofdouble Vanleer(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar);

public:
    numerical_method();
    
    void wavespeedestimate(Arrayofdouble const& p, Arrayofdouble const& q, double gamma, Arrayofdouble& wavespeed, int D, int nVar);

    void ywavespeedestimate(Arrayofdouble const& p, Arrayofdouble const& q, double gamma, Arrayofdouble& wavespeed, int D, int nVar);

    void compute_nplushalf_variables(Arrayofdouble u_n, Arrayofdouble u, Arrayofdouble u_p, double dx, double dt, double gamma, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, int nVar, int D);

    void compute_HLLCflux(amrex::Array4<amrex::Real> const& arr, Arrayofdouble& Fhllc, int nVar, double gamma, int D, int i, int j, int k, double dx, double dt);
       
    void ycompute_nplushalf_variables(Arrayofdouble u_n, Arrayofdouble u, Arrayofdouble u_p, double dx, double dt, double gamma, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, int nVar, int D);

    void ycompute_HLLCflux(amrex::Array4<amrex::Real> const& arr, Arrayofdouble& Fhllc, int nVar, double gamma, int D, int i, int j, int k, double dx, double dt);
};



