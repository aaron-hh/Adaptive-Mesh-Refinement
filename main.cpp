
#include <new>
#include <iostream>
#include <iomanip>

#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>

using namespace amrex;

//calling getLevelBld function, refers to LevelBldAdv.cpp for details
amrex::LevelBld* getLevelBld ();

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);//initialising AMREX

    Real dRunTime1 = amrex::second();

    int  max_step;
    Real strt_time;
    Real stop_time;

    {
        ParmParse pp;

        max_step  = 1000;
        strt_time =  0.0;
        stop_time = 1.0;

        //reading from input file
        pp.query("max_step",max_step);//maximum number of time step to take
        pp.query("strt_time",strt_time);
        pp.query("stop_time",stop_time);//maximum time to reach
    }

    if (strt_time < 0.0) {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0.0) {
	amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    {
        Amr amr(getLevelBld());//getlevelbld function is defined in amrlevel.cpp
        //Amr is a class defined under #include <AMReX_Amr.H>

	amr.init(strt_time,stop_time);

	//if the left-hand side expression is true, the combined result is true (the right-hand side expression is never evaluated).
	while ( amr.okToContinue() &&
  	       (amr.levelSteps(0) < max_step || max_step < 0) &&
	       (amr.cumTime() < stop_time || stop_time < 0.0) )//Cumulative measurement of execution time.

	{
	    //
	    // Do a coarse timestep.  Recursively calls timeStep()
	    //
	    amr.coarseTimeStep(stop_time);
	}

	// Write final checkpoint and plotfile
	if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
	    amr.checkPoint();
	}

    //overwrite this virtual function for format compatible with gnuplot
	if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
	    amr.writePlotFile();
	}

    }

    Real dRunTime2 = amrex::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run time = " << dRunTime2 << std::endl;

    amrex::Finalize();//ending AMREX

    return 0;
}
