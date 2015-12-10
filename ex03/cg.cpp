#include<fstream>
#include <cassert>
#include<iostream>
#include <stdlib.h>

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
    int eps  = atoi(argv[4]);
    
#ifdef USE_LIKWID
   likwid_markerInit();
   likwid_markerStartRegion( "RBGS" );
#endif

    RBGS fast(nx,ny,c);
    siwir::Timer timer;
    
    double  r = fast.solve();
    
    double time = timer.elapsed();

#ifdef USE_LIKWID
   likwid_markerStopRegion( "RBGS" );
   likwid_markerClose();
#endif

   fast.print_gnuplot();
   cout<<"time:"<<time<<'\n';
   cout<<"R:"<<r<<'\n';


}
