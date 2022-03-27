#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

typedef std::array<double,4> Arrayofdouble;

class eos
{
private:
    double m_gamma;
    double compute_pressure(double density, double v_x, double v_y, double energy, double m_gamma);
    double compute_energy(double pressure, double density, double v_x, double v_y, double m_gamma);

public:
    eos();
    //~eos();

    void set_gamma(double GAMMA);
    void compute_wavespeed(double cs, double pressure, double density, double m_gamma);
    virtual Arrayofdouble prim_to_u(Arrayofdouble prim, double gamma, int nVar, int D);
    virtual Arrayofdouble u_to_prim(Arrayofdouble u, double gamma, int nVar, int D);
    Arrayofdouble fluxf(Arrayofdouble u, Arrayofdouble prim, int nVar, int D);
    Arrayofdouble y_fluxf(Arrayofdouble u, Arrayofdouble prim, int nVar, int D);
};
