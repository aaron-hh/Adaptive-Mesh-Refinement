// 2020H - File added to perform convergence analysis

//Process text files being output from Riemann solver and AMREX to make sure text files for each refinement levels only contain data uncovered by grid above
//Read in files containing numerical solution and exact solutions simultaneously and compute the L1, L2, Linf error
//Change the number of cells at line 78. 86 according to the final cells number in the text files
//Repeat steps from no refinement to level 2 refinement and sum up the error of all levels

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

typedef std::vector<double> Vectorofdouble;
typedef std::vector<std::vector<double> > Vector_vectorofdouble;

//defining some global variables
int num_L = 3;

void compute_error(Vectorofdouble& Error, Vectorofdouble num_sol, Vectorofdouble exact_sol, double dx, int N)
{
    double sum1, sum2;
	
	double LInfnorm = 0;

	for(int i=0; i<N; i++)
	{
	sum1 += fabs(num_sol[i]-exact_sol[i]);
	sum2 += pow((num_sol[i]-exact_sol[i]),2);
	LInfnorm = std::max(LInfnorm,fabs(num_sol[i]-exact_sol[i]));	
	}

	double L1norm = 0;
	L1norm = sum1/N;
	
	double L2norm = 0;
	L2norm += sqrt(sum2/N);

    //used division by number of cells as stated in Nandan's paper instead of multiplication by dx in Toro

	Error[0] = L1norm;
	Error[1] = L2norm;
	Error[2] = LInfnorm;

}

int main()
{
    //read in files containing numerical solution for 128 cells
    std::ifstream ns_l("test1_128_l2.txt");

    Vectorofdouble v_l;

    double data1_l, data2_l;

    while(ns_l >> data1_l && ns_l >> data2_l)
    {
        v_l.push_back(data2_l);
    }

    ns_l.close();

    //read in files containing exact solution 
    std::ifstream es("test1_exact_l2.txt");

    Vectorofdouble v_exact;

    double data1_exact, data2_exact;

    while(es >> data1_exact && es >> data2_exact)
    {
        v_exact.push_back(data2_exact);
    }

    es.close();

    //calculate erros for all resolutions
    Vectorofdouble mesh(1);
    mesh[0] = 132;

    int num_mesh = mesh.size();

    Vector_vectorofdouble error_all;
    error_all.resize(num_mesh, Vectorofdouble(num_L));

    // 200 cells
    compute_error(error_all[0], v_l, v_exact, (1.0/double(132)), 132);

    //output
    std::ofstream output_convergence ("CR.dat");
    output_convergence << "# Mesh - L1Error - L2Error - LInfError" << std::endl;
            

    for(int i=0; i<num_mesh; i++)
    {
        output_convergence<<mesh[i]<<" "<<error_all[i][0]<<" "<<error_all[i][1]<<" "<<error_all[i][2]<<std::endl;
    }
}
