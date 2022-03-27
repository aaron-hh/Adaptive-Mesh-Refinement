// 2020H - File added which contains functions for HLLC solver
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>

#include "eos_v2.H"
#include "hllc.H"
#include "AMREX_Array4.H"

eos* E1 = new eos();

using namespace amrex;

//constructor
numerical_method::numerical_method()
{}

//function to find Vanleer constant
Arrayofdouble numerical_method::Vanleer(Arrayofdouble ui, Arrayofdouble uiMinus1, Arrayofdouble uiPlus1, int nVar)
{
  Arrayofdouble r;
  Arrayofdouble xi;
  Arrayofdouble Vl;

  for(int i=0; i<nVar; i++)
  {
    if (fabs(uiPlus1[i] - ui[i])>0)
    {
      r[i] = (ui[i] - uiMinus1[i]) / (uiPlus1[i] - ui[i]);
    }
    else 
    {
      r[i] = 0;
    }
    xi[i] = 2/(1+r[i]);
    //r[i] = 0;
  }

  for(int i=0; i<nVar; i++)
  {
    if(r[i]<=0)
    {
      Vl[i] = 0;
    }
    else if(r[i]>0)
    {
      if(xi[i]>(2*r[i]/(1+r[i])))
      {
        Vl[i] = 2*r[i]/(1+r[i]);
      }
      else 
      {
        Vl[i] = xi[i];
      }
    }
  }
  return Vl;
}

//computing nplushalf left and right
void numerical_method::compute_nplushalf_variables(Arrayofdouble u_n, Arrayofdouble u, Arrayofdouble u_p, double dx, double dt, double gamma, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, int nVar, int D)
{
  Arrayofdouble ui, uiMinus1, uiPlus1, Vl, Mb, delta;
  Arrayofdouble uiL, uiR, primiL, primiR, fuiL, fuiR;

  //defining ui, uiMinus1, uiPlus1
  for(int l=0; l<nVar; l++)
  {
  uiMinus1[l] = u_n[l];
  ui[l] = u[l];
  uiPlus1[l] = u_p[l];
  }

  //finding Vl
  numerical_method *Vanleer_ui = new numerical_method();
  Vl = Vanleer_ui->Vanleer(ui, uiMinus1, uiPlus1, nVar);
  double Vl_min;
  Vl_min = *std::min_element(Vl.begin(), Vl.end());

  //define delta
  for(int l=0; l<nVar; l++)
  {
    delta[l] = 0.5*(uiPlus1[l]-uiMinus1[l]);
  }

  //define uiL, uiR
  for(int l=0; l<nVar; l++)
  {
    uiL[l] = ui[l] - 0.5*Vl_min*delta[l];
    uiR[l] = ui[l] + 0.5*Vl_min*delta[l];
  }

  //convert uiL, uiR to primiL, primiR
  E1->set_gamma(gamma);
  primiL = E1->u_to_prim(uiL, gamma, nVar, D);
  primiR = E1->u_to_prim(uiR, gamma, nVar, D);

  fuiL = E1->fluxf(uiL, primiL, nVar, D);
  fuiR = E1->fluxf(uiR, primiL, nVar, D);

  //calculate uiL_nplushalf, uiR_nplushalf
  for(int l=0; l<nVar; l++)
  {
    uiL_nplushalf[l] = uiL[l] - 0.5*dt/dx*(fuiR[l]-fuiL[l]);
    uiR_nplushalf[l] = uiR[l] - 0.5*dt/dx*(fuiR[l]-fuiL[l]);
  }
}

//computing nplushalf
void numerical_method::ycompute_nplushalf_variables(Arrayofdouble u_n, Arrayofdouble u, Arrayofdouble u_p, double dx, double dt, double gamma, Arrayofdouble& uiL_nplushalf, Arrayofdouble& uiR_nplushalf, int nVar, int D)
{
  Arrayofdouble ui, uiMinus1, uiPlus1, Vl, Mb, delta;
  Arrayofdouble uiL, uiR, primiL, primiR, fuiL, fuiR;

  //defining ui, uiMinus1, uiPlus1
  for(int l=0; l<nVar; l++)
  {
  uiMinus1[l] = u_n[l];
  ui[l] = u[l];
  uiPlus1[l] = u_p[l];
  }

  //finding Vl
  numerical_method *Vanleer_ui = new numerical_method();
  Vl = Vanleer_ui->Vanleer(ui, uiMinus1, uiPlus1, nVar);
  double Vl_min;
  Vl_min = *std::min_element(Vl.begin(), Vl.end());

  //define delta
  for(int l=0; l<nVar; l++)
  {
    delta[l] = 0.5*(uiPlus1[l]-uiMinus1[l]);
  }

  //define uiL, uiR
  for(int l=0; l<nVar; l++)
  {
    uiL[l] = ui[l] - 0.5*Vl_min*delta[l];
    uiR[l] = ui[l] + 0.5*Vl_min*delta[l];
  }

  //convert uiL, uiR to primiL, primiR
  E1->set_gamma(gamma);
  primiL = E1->u_to_prim(uiL, gamma, nVar, D);
  primiR = E1->u_to_prim(uiR, gamma, nVar, D);

  fuiL = E1->y_fluxf(uiL, primiL, nVar, D);
  fuiR = E1->y_fluxf(uiR, primiL, nVar, D);

  //calculate uiL_nplushalf, uiR_nplushalf
  for(int l=0; l<nVar; l++)
  {
    uiL_nplushalf[l] = uiL[l] - 0.5*dt/dx*(fuiR[l]-fuiL[l]);
    uiR_nplushalf[l] = uiR[l] - 0.5*dt/dx*(fuiR[l]-fuiL[l]);
  }
}

void numerical_method::wavespeedestimate(Arrayofdouble const& a, Arrayofdouble const&  b, double gamma, Arrayofdouble& wavespeed, int D, int nVar)
{
  Arrayofdouble P, Q;
  P = E1->u_to_prim(a, gamma, nVar, D);
  Q = E1->u_to_prim(b, gamma, nVar, D);

  double rho_l = P[0];
  double rho_r = Q[0];

  double vx_l = P[1];
  double vx_r = Q[1];

  double p_l = P[D+1];
  double p_r = Q[D+1];

  double a_l = sqrt(gamma*p_l/rho_l);
  double a_r = sqrt(gamma*p_r/rho_r);

  double S_lr = std::max(fabs(vx_l)+a_l, fabs(vx_r)+a_r);
  double Sl = - S_lr;
  double Sr = + S_lr;

  double S_star = (p_r-p_l+((rho_l*vx_l)*(Sl-vx_l))-((rho_r*vx_r)*(Sr-vx_r)))/((rho_l*(Sl-vx_l))-(rho_r*(Sr-vx_r)));

  wavespeed[0] = Sl;
  wavespeed[1] = Sr;
  wavespeed[2] = S_star;
}

void numerical_method::ywavespeedestimate(Arrayofdouble const& a, Arrayofdouble const&  b, double gamma, Arrayofdouble& wavespeed, int D, int nVar)
{
  Arrayofdouble P, Q;
  P = E1->u_to_prim(a, gamma, nVar, D);
  Q = E1->u_to_prim(b, gamma, nVar, D);

  double rho_l = P[0];
  double rho_r = Q[0];

  double vx_l = P[1];
  double vx_r = Q[1];

  double p_l = P[D+1];
  double p_r = Q[D+1];

  double a_l = sqrt(gamma*p_l/rho_l);
  double a_r = sqrt(gamma*p_r/rho_r);

  double vy_l;
  double vy_r;

  if(D == 1)
  {
  vy_l = 0;
  vy_r = 0;
  }
  else if(D == 2)
  {
  vy_l = P[D];
  vy_r = Q[D];
  }

  double S_lr = std::max(fabs(vy_l)+a_l, fabs(vy_r)+a_r);
  double Sl = - S_lr;
  double Sr = + S_lr;

  double S_star = (p_r-p_l + ((rho_l*vy_l)*(Sl-vy_l)) - ((rho_r*vy_r)*(Sr-vy_r))) / ((rho_l*(Sl-vy_l))-(rho_r*(Sr-vy_r)));

  wavespeed[0] = Sl;
  wavespeed[1] = Sr;
  wavespeed[2] = S_star;
}

void numerical_method::compute_HLLCflux(amrex::Array4<amrex::Real> const& arr, Arrayofdouble& Fhllc, int nVar, double gamma, int D, int i, int j, int k, double dx, double dt)
{
  //defining u_n, u, u_p for input into compute_nplushalf_variables left
  Arrayofdouble ul_n, ul, ul_p;
  for(int l=0; l<nVar; l++)
  {
  //ui-1
  ul_n[l] = arr(i-2,j,k,l);

  //ui
  ul[l] = arr(i-1,j,k,l);

  //ui+1
  ul_p[l] = arr(i,j,k,l);
  }

  //std::cout<<ul[0]<<std::endl;

  //calling compute_nplushalf function
  Arrayofdouble uiL_nplushalf_l, uiR_nplushalf_l;
  compute_nplushalf_variables(ul_n, ul, ul_p, dx, dt, gamma, uiL_nplushalf_l, uiR_nplushalf_l, nVar, D);

  //defining u_n, u, u_p for input into compute_nplushalf_variables right
  Arrayofdouble ur_n, ur, ur_p;
  for(int l=0; l<nVar; l++)
  {
  //ui-1
  ur_n[l] = arr(i-1,j,k,l);

  //ui
  ur[l] = arr(i,j,k,l);

  //ui+1
  ur_p[l] = arr(i+1,j,k,l);
  }

  //std::cout<<ur[3]<<std::endl;

  //calling compute_nplushalf function
  Arrayofdouble uiL_nplushalf_r, uiR_nplushalf_r;
  compute_nplushalf_variables(ur_n, ur, ur_p, dx, dt, gamma, uiL_nplushalf_r, uiR_nplushalf_r, nVar, D);

  //defining a,b for HLLC solver
  Arrayofdouble a, b;
  for(int l=0; l<nVar; l++)
  {
    a[l] = uiR_nplushalf_l[l];
    b[l] = uiL_nplushalf_r[l];
  }

  //std::cout<<a[3]<<std::endl;

  //calling the compute wavespeed function
  Arrayofdouble wavespeed;
  wavespeedestimate(a, b, gamma, wavespeed, D, nVar);

  //std::cout<<wavespeed[0]<<std::endl;

  //define input variables 
  double Sl = wavespeed[0];
  double Sr = wavespeed[1];
  double S_star = wavespeed[2];

  //convert conservative to primitive variables
  Arrayofdouble P, Q;
  P = E1->u_to_prim(a, gamma, nVar, D);
  Q = E1->u_to_prim(b, gamma, nVar, D);

  double rho_l = P[0];
  double vx_l = P[1];
  double U_l = a[D+1];
  double p_l = P[D+1];

  double rho_r = Q[0];
  double vx_r = Q[1];
  double U_r = b[D+1];
  double p_r = Q[D+1];

  Arrayofdouble Fl,Fr;
  E1->set_gamma(gamma);
  Fl = E1->fluxf(a, P, nVar, D);
  Fr = E1->fluxf(b, Q, nVar, D);

  //compute Ql* & Fl*, Qr* & Fr*
  Arrayofdouble Ql_star,Qr_star,Q_hllc,prim_hllc;
  Arrayofdouble Fl_star,Fr_star;

  Ql_star[0] = rho_l*((Sl-vx_l)/(Sl-S_star));
  Ql_star[1] = rho_l*((Sl-vx_l)/(Sl-S_star))*(S_star);
  Ql_star[D+1] = rho_l*((Sl-vx_l)/(Sl-S_star))*((U_l/rho_l)+(S_star-vx_l)*(S_star+(p_l/(rho_l*(Sl-vx_l)))));

  Qr_star[0] = rho_r*((Sr-vx_r)/(Sr-S_star));
  Qr_star[1] = rho_r*((Sr-vx_r)/(Sr-S_star))*(S_star);
  Qr_star[D+1] = rho_r*((Sr-vx_r)/(Sr-S_star))*((U_r/rho_r)+(S_star-vx_r)*(S_star+(p_r/(rho_r*(Sr-vx_r)))));

  if(D == 2)
  {
  double vy_l = P[2];
  double vy_r = Q[2];

  Ql_star[2] = rho_l*((Sl-vx_l)/(Sl-S_star))*(vy_l);
  Qr_star[2] = rho_r*((Sr-vx_r)/(Sr-S_star))*(vy_r);
  }

  for(int l=0; l<nVar; l++)
  {
  Fl_star[l] = Fl[l]  + Sl*(Ql_star[l] - a[l]);
  Fr_star[l] = Fr[l]  + Sr*(Qr_star[l] - b[l]);
  }
	
	if(Sl >= 0)
	{	
		Fhllc = Fl;
		return; 
	}

	if(Sl < 0 && S_star >= 0)
	{
		Fhllc = Fl_star;
		return;
	}
	
	if(S_star < 0 && Sr > 0)
	{
		Fhllc = Fr_star;
		return;
	}

	if(Sr <= 0)
	{	Fhllc = Fr;
		return; 
	}		
}

void numerical_method::ycompute_HLLCflux(amrex::Array4<amrex::Real> const& arr, Arrayofdouble& Fhllc, int nVar, double gamma, int D, int i, int j, int k, double dx, double dt)
{
  //defining u_n, u, u_p for input into compute_nplushalf_variables left
  Arrayofdouble ul_n, ul, ul_p;
  for(int l=0; l<nVar; l++)
  {
  //ui-1
  ul_n[l] = arr(i,j-2,k,l);

  //ui
  ul[l] = arr(i,j-1,k,l);

  //ui+1
  ul_p[l] = arr(i,j,k,l);
  }

  //calling compute_nplushalf function
  Arrayofdouble uiL_nplushalf_l, uiR_nplushalf_l;
  ycompute_nplushalf_variables(ul_n, ul, ul_p, dx, dt, gamma, uiL_nplushalf_l, uiR_nplushalf_l, nVar, D);

  //defining u_n, u, u_p for input into compute_nplushalf_variables right
  Arrayofdouble ur_n, ur, ur_p;
  for(int l=0; l<nVar; l++)
  {
  //ui-1
  ur_n[l] = arr(i,j-1,k,l);

  //ui
  ur[l] = arr(i,j,k,l);

  //ui+1
  ur_p[l] = arr(i,j+1,k,l);
  }

  //calling compute_nplushalf function
  Arrayofdouble uiL_nplushalf_r, uiR_nplushalf_r;
  ycompute_nplushalf_variables(ur_n, ur, ur_p, dx, dt, gamma, uiL_nplushalf_r, uiR_nplushalf_r, nVar, D);

  //defining a,b for HLLC solver
  Arrayofdouble a, b;
  for(int l=0; l<nVar; l++)
  {
    a[l] = uiR_nplushalf_l[l];
    b[l] = uiL_nplushalf_r[l];
  }

  //calling the compute wavespeed function
  Arrayofdouble wavespeed;
  ywavespeedestimate(a, b, gamma, wavespeed, D, nVar);

  //define input variables 
  double Sl = wavespeed[0];
  double Sr = wavespeed[1];
  double S_star = wavespeed[2];

  //convert conservative to primitive variables
  Arrayofdouble P, Q;
  P = E1->u_to_prim(a, gamma, nVar, D);
  Q = E1->u_to_prim(b, gamma, nVar, D);

  double rho_l = P[0];
  double vx_l = P[1];
  double U_l = a[D+1];
  double p_l = P[D+1];

  double rho_r = Q[0];
  double vx_r = Q[1];
  double U_r = b[D+1];
  double p_r = Q[D+1];

  Arrayofdouble Fl,Fr;
  E1->set_gamma(gamma);
  Fl = E1->y_fluxf(a, P, nVar, D);
  Fr = E1->y_fluxf(b, Q, nVar, D);

  //compute Ql* & Fl*, Qr* & Fr*
  Arrayofdouble Ql_star,Qr_star,Q_hllc,prim_hllc;
  Arrayofdouble Fl_star,Fr_star;

  double vy_l = P[2];
  double vy_r = Q[2];

  Ql_star[2] = rho_l*((Sl-vy_l)/(Sl-S_star))*(S_star);
  Qr_star[2] = rho_r*((Sr-vy_r)/(Sr-S_star))*(S_star);

  Ql_star[0] = rho_l*((Sl-vy_l)/(Sl-S_star));
  Ql_star[1] = rho_l*((Sl-vy_l)/(Sl-S_star))*(vx_l);
  Ql_star[D+1] = rho_l*((Sl-vy_l)/(Sl-S_star))*((U_l/rho_l)+(S_star-vy_l)*(S_star+(p_l/(rho_l*(Sl-vy_l)))));

  Qr_star[0] = rho_r*((Sr-vy_r)/(Sr-S_star));
  Qr_star[1] = rho_r*((Sr-vy_r)/(Sr-S_star))*(vx_r);
  Qr_star[D+1] = rho_r*((Sr-vy_r)/(Sr-S_star))*((U_r/rho_r)+(S_star-vy_r)*(S_star+(p_r/(rho_r*(Sr-vy_r)))));

  for(int l=0; l<nVar; l++)
  {
  Fl_star[l] = Fl[l]  + Sl*(Ql_star[l] - a[l]);
  Fr_star[l] = Fr[l]  + Sr*(Qr_star[l] - b[l]);
  }

  if(Sl >= 0)
  {	
    Fhllc = Fl;
    return; 
  }

  if(Sl <0 && S_star >= 0)
  {
    Fhllc = Fl_star;
    return;
  }
  
  if(S_star < 0 && Sr > 0)
  {
    Fhllc = Fr_star;
    return;
  }

  if(Sr <= 0)
  {	Fhllc = Fr;
    return; 
  }		
}



