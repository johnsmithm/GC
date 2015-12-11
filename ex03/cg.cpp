#include<fstream>
#include <cassert>
#include<iostream>
#include <stdlib.h>


#include"Vector.h"  //Expr_CG function
#include "cg_naive.h" //CG class

#include "CG.h"
#include "Timer.h"

#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}

#endif
using namespace std;


int main(int argc, char *argv[]){
    (void) argc; //to suppress Warnings about unused argc
    assert(argc>3);
    int nx   = atoi(argv[1]);
    int ny   = atoi(argv[2]);
    int c    = atoi(argv[3]);	
    int eps  = stod(argv[4]);
    
#ifdef USE_LIKWID
   likwid_markerInit();
   likwid_markerStartRegion( "CG" );
#endif

   
    siwir::Timer timer;
    
    //using simple function operation on vectors
    CG fast(nx,ny,c,eps);
    double  r = fast.solve_naive();
	
	//using template expretions
    //double r = Expr_CG(nx,ny,c,eps);
    
    double time = timer.elapsed();

#ifdef USE_LIKWID
   likwid_markerStopRegion( "CG" );
   likwid_markerClose();
#endif

   fast.print_gnuplot();
   cout<<"time:"<<time<<'\n';
   cout<<"R:"<<r<<'\n';


}
