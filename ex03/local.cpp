#include<fstream>
#include <cassert>
#include<iostream>
#include <time.h>
#include <stdlib.h> 

#include"Vector.h"  //Expr_CG function
#include "cg_naive.h" //CG class

using namespace std;




int main(int argc, char *argv[]){
	
    assert(argc>3);
    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int c  = atoi(argv[3]);	
    double eps  = stod(argv[4]);
    //cout<<eps<<" "<<ny<<" "<<c<<endl;

    clock_t t1,t2;
    t1=clock();
    
    CG fast(nx,ny,c,eps);
    double  r = fast.solve_naive();
	
    //double r = Expr_CG(nx,ny,c,eps);
	
    t2=clock();
    float diff ((float)t2-(float)t1);
    float time = diff / CLOCKS_PER_SEC;

   //fast.print_gnuplot_naive();

  
    cout<<"time:"<<time<<'\n';
    cout<<"R:"<<r<<'\n';


}