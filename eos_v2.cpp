// 2020H - File added which contains functions to compute flux and convert variables
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "eos_v2.H"

typedef std::array<double,4> Arrayofdouble;
//one-D: rho, vx, p, 0
//two-D: rho, vx, vy, p

//constructor 
eos::eos()
{}

//setter for gamma 
void eos::set_gamma(double GAMMA)
{
    m_gamma = GAMMA;
}

//defining compute pressure function
double eos::compute_pressure(double density, double v_x, double v_y, double energy, double m_gamma)
{
    double pressure = (m_gamma-1) * (energy - 0.5*density*((v_x*v_x) + (v_y*v_y)));

    return pressure;
}

//defining compute energy function
double eos::compute_energy(double pressure, double density, double v_x, double v_y, double m_gamma)
{   
    double energy = pressure / (m_gamma-1) + 0.5*density*((v_x*v_x) + (v_y*v_y));

    return energy;
}

//defining primitive to conservative function
Arrayofdouble eos::prim_to_u(Arrayofdouble prim, double m_gamma, int nVar, int D) 
{
    Arrayofdouble u;
    eos* E = new eos();

    u[0] = prim[0];
    u[1] = prim[1] * prim[0];
    
    if(D==2)
    {
    u[D] = prim[D] * prim[0];
    u[D+1] = E->compute_energy(prim[D+1], prim[0], prim[1], prim[D], m_gamma);
    }
    else if(D==1)
    {
    u[D+1] = E->compute_energy(prim[D+1], prim[0], prim[1], 0, m_gamma);
    }

    return u;
}

//defining conservative to primitive function
Arrayofdouble eos::u_to_prim(Arrayofdouble u, double m_gamma, int nVar, int D)
{  
    Arrayofdouble prim;
    eos* E = new eos();

    prim[0] = u[0];
    prim[1] = u[1]/u[0];

    if(D==2)
    {
    prim[D] = u[D]/u[0];
    prim[D+1] = E->compute_pressure(prim[0], prim[1], prim[D], u[D+1], m_gamma);
    }
    else if(D==1)
    {
    prim[D+1] = E->compute_pressure(prim[0], prim[1], 0, u[D+1], m_gamma);
    }

    return prim;
}

//defining flux function
Arrayofdouble eos::fluxf(Arrayofdouble u, Arrayofdouble prim, int nVar, int D)
{
    Arrayofdouble x_fluxfunction;

    x_fluxfunction[0] = u[1];
    x_fluxfunction[1] = u[1] * prim[1] + prim[D+1];
    x_fluxfunction[D+1] = (u[D+1] + prim[D+1]) * prim[1]; 

    if(D==2)
    {
    x_fluxfunction[2] = prim[0] * prim[1] * prim[2];
    }

    return x_fluxfunction;
}

//defining flux function
Arrayofdouble eos::y_fluxf(Arrayofdouble u, Arrayofdouble prim, int nVar, int D)
{
  Arrayofdouble y_fluxfunction;

  y_fluxfunction[0] = u[2];
  y_fluxfunction[1] = prim[0] * prim[1] * prim[2];
  y_fluxfunction[2] = u[2] * prim[2] + prim[3]; 
  y_fluxfunction[3] = (u[3] + prim[3]) * prim[2]; 

  return y_fluxfunction;
}

