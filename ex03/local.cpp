#include<fstream>
#include <cassert>
#include<iostream>
#include <time.h>
#include <stdlib.h> 
#include <string>
#include <mpi.h>

#include"Vector.h"  //Expr_CG function
#include "cg_naive.h" //CG class

using namespace std;




int main(int argc, char *argv[]){
	int size(0); 
    int rank(0);
	
	MPI_Init( &argc, &argv );
	
	MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	
    assert(argc>3);
    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int c  = atoi(argv[3]);	
    double eps  = stod(argv[4]);
    //cout<<eps<<" "<<ny<<" "<<c<<endl;

    clock_t t1,t2;
    t1=clock();
    
	//using simple function operation on vectors
    //CG fast(nx,ny,c,eps);
    //double  r = fast.solve_naive();
	
	//using template expretions
    double r = Expr_CG(nx,ny,c,eps);
	
	
    t2=clock();
    float diff ((float)t2-(float)t1);
    float time = diff / CLOCKS_PER_SEC;

   //fast.print_gnuplot_naive();

    if(rank==0){
    cout<<"time:"<<time<<'\n';
    cout<<"R:"<<r<<'\n';
	}

    MPI_Finalize();
}