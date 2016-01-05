#include<fstream>
#include <cassert>
#include<iostream>
#include <stdlib.h>

#include <mpi.h>

#include <mpi.h>
#include"Vector.h"  //Expr_CG function
//#include "cg_naive.h" //CG class

#include "Timer.h"

#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}

#endif
using namespace std;


int main(int argc, char *argv[]){ 
    
    (void) argc; //to suppress Warnings about unused argc
    int size(0); 
    int rank(0);
	
	MPI_Init( &argc, &argv );
	
	MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	
    assert(argc>3);
    int nx   = atoi(argv[1]);
    int ny   = atoi(argv[2]);
    int c    = atoi(argv[3]);	
    double eps  = stod(argv[4]);
    
#ifdef USE_LIKWID
   likwid_markerInit();
   likwid_markerStartRegion( "CG" );
#endif

    
   
    siwir::Timer timer;    
	
	//using template expretions
    double r = Expr_CG(nx,ny,c,eps,rank,size);

    
    double time = timer.elapsed();

#ifdef USE_LIKWID
   likwid_markerStopRegion( "CG" );
   likwid_markerClose();
#endif
    
   //fast.print_gnuplot();
    if(rank==0){
       cout<<"time:"<<time<<'\n';
       cout<<"R:"<<r<<'\n';
    }
 MPI_Finalize();
}
