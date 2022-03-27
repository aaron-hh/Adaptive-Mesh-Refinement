// 2020H - File added to compute exact solution for Euler system

//Solving conservative form Eulers equation using Riemann exact solver
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <array>

/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Code for solving Euler using Riemann exact                */
/*                                                                    */
/* -------------------------------------------------------------------*/

//function for pressure guess
void pressure_guess(double gamma, double& p_old, std::array<double,3> priml, std::array<double,3> primr)
{
  //defining variables
  double rhol = priml[0];
  double vl = priml[1];
  double pl = priml[2]; 

  double rhor = primr[0]; 
  double vr = primr[1]; 
  double pr = primr[2]; 

  //compute Csl, Csr
  double Csl = sqrt(gamma*pl/rhol);
  double Csr = sqrt(gamma*pr/rhor);

  //logic for choosing pressure approximation
  int Quser = 2;
  double rho_a = 0.25*(rhol+rhor)*(Csl+Csr);
  double ppv = std::max(0.0, 0.5*(pl+pr) + 0.5*(vl-vr)*rho_a);
  double p_min = std::min(pl,pr);
  double p_max = std::max(pl,pr);
  double qmax = p_max/p_min;

  bool Cond1 = (qmax <= Quser) && ((p_min <= ppv) && (p_max >= ppv));
        
  if(Cond1)
    {
      //Primitive variable pressure approximation
      p_old = ppv; 
    }
  else
    {		
      if(ppv < p_min)
		  {	
        //Two rarefaction pressure approximation
        double top = Csl + Csr - 0.5*(gamma-1)*(vr-vl);
        double bottom_l = Csl / pow(pl,((gamma-1)/(2*gamma))); 
        double bottom_r = Csr / pow(pr,((gamma-1)/(2*gamma))); 
        p_old = pow((top/(bottom_l + bottom_r)), (2*gamma/(gamma-1)));
		  }
      else
      {		
        //Two shock pressure approximation
        double Al = 2/((gamma+1)*rhol);
        double Bl = (gamma-1)*pl/(gamma+1);
        double p0 = std::max(10e-6, ppv);
        double gl = sqrt(Al/(p0+Bl));
        
        double Ar = 2/((gamma+1)*rhor);
        double Br = (gamma-1)*pr/(gamma+1);
        double gr = sqrt(Ar/(p0+Br));

        p_old = (gl*pl + gr*pr - (vr-vl))/(gl+gr);
      }		
	  }
}

//sampling solution
std::array<double,3> sample(double gamma, double& p_old, std::array<double,3> priml, std::array<double,3> primr, double v_star, double S, int wave)
{

  //defining variables
  double rhol = priml[0];
  double vl = priml[1];
  double pl = priml[2]; 

  double rhor = primr[0]; 
  double vr = primr[1]; 
  double pr = primr[2]; 

  double p_new = p_old;

  //compute Csl, Csr
  double Csl = sqrt(gamma*pl/rhol);
  double Csr = sqrt(gamma*pr/rhor);

  //compute rhol*, vl*
  double Sl, rhol_star, vl_star;
  if (p_new>pl)//shock
  {
    Sl = vl - Csl*sqrt(((gamma+1)/(2*gamma))*(p_new/pl)+((gamma-1)/(2*gamma)));
    rhol_star = rhol*((p_new/pl)+((gamma-1)/(gamma+1)))/(((gamma-1)/(gamma+1))*(p_new/pl)+1);
  }
  else if (p_new<=pl)//rarefaction
  {
    rhol_star = rhol*pow((p_new/pl),(1/gamma));
  }

  //compute rhor*, vr* 
  double Sr, rhor_star, vr_star;
  if (p_new>pr)//shock
  {
    Sr = vr + Csr*sqrt(((gamma+1)/(2*gamma))*(p_new/pr)+((gamma-1)/(2*gamma)));
    rhor_star = rhor*((p_new/pr)+((gamma-1)/(gamma+1)))/(((gamma-1)/(gamma+1))*(p_new/pr)+1);
  }
  else if (p_new<=pr)//rarefaction
  {
    rhor_star = rhor*pow((p_new/pr),(1/gamma));
  }

  double Shl = vl - Csl;
  double Csl_star = Csl*pow((p_new/pl),((gamma-1)/(2*gamma)));
  double Stl = v_star - Csl_star;
  double Shr = vr + Csr;
  double Csr_star = Csr*pow((p_new/pr),((gamma-1)/(2*gamma)));
  double Str = v_star + Csr_star;

  //compute rho, v, pressure for intermediate state across a rarefaction fan
  double rhol_raref = rhol*pow(((2/(gamma+1))+((gamma-1)/((gamma+1)*Csl))*(vl-S)),(2/(gamma-1)));
  double vl_raref = (2/(gamma+1))*(Csl+((gamma-1)*vl/2)+S);
  double pl_raref = pl*pow(((2/(gamma+1))+((gamma-1)/((gamma+1)*Csl))*(vl-S)),((2*gamma)/(gamma-1)));

  double rhor_raref = rhor*pow(((2/(gamma+1))-((gamma-1)/((gamma+1)*Csr))*(vr-S)),(2/(gamma-1)));
  double vr_raref = (2/(gamma+1))*(-Csr+((gamma-1)*vr/2)+S);
  double pr_raref = pr*pow(((2/(gamma+1))-((gamma-1)/((gamma+1)*Csr))*(vr-S)),((2*gamma)/(gamma-1)));

  //compute u for left contact
  std::array<double,3> ul;
  ul[0] = priml[0];
  ul[1] = priml[1];
  ul[2] = priml[2];

  std::array<double,3> ul_fan;
  ul_fan[0] = rhol_raref;
  ul_fan[1] = vl_raref;
  ul_fan[2] = pl_raref;

  std::array<double,3> ul_star_shock;
  ul_star_shock[0] = rhol_star;
  ul_star_shock[1] = v_star;
  ul_star_shock[2] = p_new;

  std::array<double,3> ul_star_fan;
  ul_star_fan[0] = rhol*pow((p_new/pl),(1/gamma));
  ul_star_fan[1] = vl_star;
  ul_star_fan[2] = p_new;

  //compute u for right contact
  std::array<double,3> ur;
  ur[0] = primr[0];
  ur[1] = primr[1];
  ur[2] = primr[2];

  std::array<double,3>ur_fan;
  ur_fan[0] = rhor_raref;
  ur_fan[1] = vr_raref;
  ur_fan[2] = pr_raref;

  std::array<double,3> ur_star_shock;
  ur_star_shock[0] = rhor_star;
  ur_star_shock[1] = v_star;
  ur_star_shock[2] = p_new;

  std::array<double,3> ur_star_fan;
  ur_star_fan[0] = rhor*pow((p_new/pr),(1/gamma));
  ur_star_fan[1] = vr_star;
  ur_star_fan[2] = p_new;

  //logic for computing exact solution
  std::array<double,3> u;

  switch(wave)
  {
    case 1: //rarefaction-contact-rarefaction
    {
      if(S < Shl)
      {
        u = ul;
      }
      else if(S >= Shl && S < Stl)
      {
        u = ul_fan;
      }
      else if(S >= Stl && S < v_star)
      {
        u = ul_star_shock;
      }
      else if(S >= v_star && S < Str)
      {
        u = ur_star_shock;
      }
      else if(S >= Str && S < Shr)
      {
        u = ur_fan;
      }
      else if(S >= Shr)
      {
        u = ur;
      }

      break;
    }

    case 2: //rarefaction-contact-shock
    {
      if (S < Shl)
			{ 
        u = ul;				
      }
      else if(S >= Shl && S < Stl)
      { 
        u = ul_fan;
      }	
      else if(S >= Stl && S < v_star)
      { 				
        u = ul_star_shock;
      }		    
      else if(S >= v_star && S < Sr)
      {				
        u = ur_star_shock;				
      }		    		    
      else if(S >= Sr)
      {			
        u = ur;	
      }	

      break;
    }

    case 3: //shock-contact-rarefaction
    {
      if(S < Sl)
			{ 
        u = ul;					
      }	   		    
      else if(S >= Sl && S < v_star)
      {		
        u = ul_star_shock;
      }
      else if(S >= v_star && S < Str)
      { 				
        u = ur_star_shock;
      }		    
      else if(S >= Str && S < Shr)
      {	
        u = ur_fan;
      }		    
      else if(S >= Shr)
      {			
        u = ur;	
      }			
			
			break;			
    }

    case 4: //shock-contact-shock
    {
			if(S < Sl)
			{
        u = ul;			
      }	    
            
      else if(S >= Sl && S < v_star)
      { 			
        u = ul_star_shock;				
      }		    
      else if(S >= v_star && S < Sr )
      {	
        u = ur_star_shock;
      }		    
      else if(S >= Sr)
      {		
        u = ur;			
      }
    
			break;	
    }
  }
  return u;
}

int main()
{
double nCells = 128;
double x0 = 0;
double x1 = 1;
double t0 = 0.0;
double t1 = 0.15;
double C = 0.9;
double gamma = 1.4;
double dx = (x1-x0)/nCells;

std::vector<std::array<double, 3> >prim(nCells+2);//an array of primitive elements
std::vector<std::array<double, 3> >u(nCells+2);//an array of conservative elements

//setting initial data in primitive form
for(int i=0; i<nCells+2; i++)
{
  double x = x0 + (i-1)*dx;

  if(x<=0.5)
  {
    prim[i][0] = 1.0;
    prim[i][1] = -2.0;
    prim[i][2] = 0.4;
  }
  else
  {
    prim[i][0] = 1.0;
    prim[i][1] = 2.0;
    prim[i][2] = 0.4;
  }

  // if(x<=0.5)
  // {
  //   prim[i][0] = 5.99924;
  //   prim[i][1] = 19.5975;
  //   prim[i][2] = 460.894;
  // }
  // else
  // {
  //   prim[i][0] = 5.99242;
  //   prim[i][1] = -6.19633;
  //   prim[i][2] = 46.0950;
  // }
}

//compute p* using root finding Newton Raphson
double rhol = prim[0][0];
double vl = prim[0][1];
double pl = prim[0][2]; 

double rhor = prim[nCells+1][0]; 
double vr = prim[nCells+1][1]; 
double pr = prim[nCells+1][2]; 

double p_old = 0.0;//guess pressure
double p_new = 0.0;
double plfluxfunc;
double prfluxfunc;
double pfluxfunc;
double plfluxfunc_d;
double prfluxfunc_d;
double pfluxfunc_d;

//compute Csl, Csr
double Csl = sqrt(gamma*pl/rhol);
double Csr = sqrt(gamma*pr/rhor);

//call the compute pressure function here
pressure_guess(gamma, p_old, prim[0], prim[nCells+1]);

//Newton Raphson
do
{
//compute flux function left
if(p_old>pl)//shock
{
  double Al = 2/((gamma+1)*rhol);
  double Bl = (gamma-1)*pl/(gamma+1);
  plfluxfunc = (p_old-pl)*pow(Al/(p_old+Bl),0.5);
  plfluxfunc_d = pow((Al/(Bl+p_old)),0.5)*(1-0.5*((p_old-pl)/((Bl+p_old))));
}
else if (p_old<=pl)//rarefaction
{
  plfluxfunc = ((2*Csl)/(gamma-1))*(pow((p_old/pl),((gamma-1)/(2*gamma)))-1); 
  plfluxfunc_d = (1/(rhol*Csl))*pow((p_old/pl),(-(gamma+1)/(2*gamma)));
}

//compute flux function right
if(p_old>pr)
{
  double Ar = 2/((gamma+1)*rhor);
  double Br = (gamma-1)*pr/(gamma+1);
  prfluxfunc = (p_old-pr)*pow(Ar/(p_old+Br),0.5);
  prfluxfunc_d = pow((Ar/(Br+p_old)),0.5)*(1-0.5*((p_old-pr)/((Br+p_old))));
}
else if (p_old<=pr)
{
  prfluxfunc = ((2*Csr)/(gamma-1))*(pow((p_old/pr),((gamma-1)/(2*gamma)))-1);
  prfluxfunc_d = (1/(rhor*Csr))*pow((p_old/pr),(-(gamma+1)/(2*gamma)));  
}

//compute flux function
pfluxfunc = plfluxfunc + prfluxfunc + vr - vl;
pfluxfunc_d = plfluxfunc_d + prfluxfunc_d;

//compute new pressure
p_new = p_old - (pfluxfunc/pfluxfunc_d);
p_old = std::max(p_new, 1e-16);

} while (2.0*std::fabs((p_new-p_old)/(p_new+p_old))>1e-9);

//compute pstar and vstar
p_new = p_old;//p_new is p_star
//p_new = 1691.64;
double v_star = 0.5*(vl+vr) + 0.5*(prfluxfunc - plfluxfunc);

//std::cout<<p_new<<std::endl;

//Riemann exact solver
std::array<double,512> xd; 
std::array<double,512> S;

for(int i=1; i<nCells; i++)
{
  xd[0] = dx;
  xd[i] = xd[i-1] + dx;
}

//deciding wave pattern based on presence of shock / rarefaction in solution
int wave; 	

if(p_new < pl && p_new < pr)
{
  // The pattern is rarefaction-contact-rarefaction		
  wave = 1; 		
} 
else if(p_new < pl && p_new >= pr)
{	
  // The pattern is rarefaction-contact-shock
  wave = 2; 
}	
else if(p_new >= pl && p_new < pr)
{	
  // The pattern is shock-contact-rarefaction
  wave = 3 ; 
} 
else if(p_new >= pl && p_new >= pr)
{	
  // The pattern is shock-contact-shock
  wave = 4 ; 
} 

for(int i=0; i<nCells+2; i++)
{
double xd0 = 0.5;
double t = t1-t0;
S[i] = (xd[i]-xd0)/t;

u[i] = sample(gamma, p_new, prim[0], prim[nCells+1], v_star, S[i], wave);
}

//outputing data
std::ofstream output ("Riemann.txt");
for(int i=0; i<nCells; i++)
{
  double x = x0 + (i+0.5)*dx;
  output<<x<<" "<<u[i][0]<<std::endl;
}
}
