#include<fstream>
#include <cassert>
#include<iostream>
#include <time.h>
//#include "cg_naive.h"
#include <string> 
#include <stdlib.h> 
#include"Vector.h"

using namespace std;

double Expr_CG(int nx,int ny,int c,double eps){
	//initialization	
	int pg = (1+nx)*(ny+1);	
	Vector r(pg), d(pg), z(pg); 
	       double pi = 3.141592653589793;
	       double hx_ = 2.0/nx;
           double hy_ = 1.0/ny;
		   double C = 4*pi*pi;   
		   double freqx = 2*pi*hx_;   
		   double freqy = 2*pi*hy_;     
		   //4π^2 sin(2πx) sinh(2πy)
	Vector f(nx+1,ny+1,[C,freqx,freqy](int x,int y)->double{return C*sin(freqx*x)*sinh(freqy*y);});	
			int last_row = (ny)*(nx+1);	
			double SINH = sinh(2*pi); 
		   //sin(2πx) sinh(2πy)	
	Vector u(pg,last_row,last_row+nx+1,[SINH,freqx](int x)->double{return sin(x*freqx) * SINH;});	
	Stencil A(nx,ny);
	double delta0 = 0, delta1 = 0, beta = 0,alfa=0;
	//initialization
	//CG
	r = f - A*u;
	std::cout<<u<<(A*u);
	delta0 = r.LNorm();
	if(sqrt(delta0)<eps)return sqrt(delta0);
	
	d = r;return 0;
	for(int i=0;i<c;++i){
		z = A*d;
		
		alfa = delta0 / (d^z);
		u = u + d*alfa;
		r = r - z*alfa;
		delta1 = r.LNorm();
		if(sqrt(delta1)<eps)return sqrt(delta1);
		beta = delta1/delta0;
		d = r + d*beta;
		delta0=delta1;
	}
	return sqrt(delta0);
	
}


int main(int argc, char *argv[]){
	/*int size_ = 5;
	Vector a(size_,1.), b(size_,2.),c(size_);
	//double d = 3.0;
	c = a - b  ;
	std::cout<<c<<"\n";
	return 0;*/
    assert(argc>3);
    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int c  = atoi(argv[3]);	
    double eps  = stod(argv[4]);
    //cout<<eps<<" "<<ny<<" "<<c<<endl;

    clock_t t1,t2;
    t1=clock();
    
    //CG fast(nx,ny,c,eps);
    //double  r1 = fast.solve_naive();
    //double  r = fast.solve_naive();
    double r = Expr_CG(nx,ny,c,eps);
    t2=clock();
    float diff ((float)t2-(float)t1);
    float time = diff / CLOCKS_PER_SEC;

   //fast.print_gnuplot_naive();

  
    cout<<"time:"<<time<<'\n';
    cout<<"R:"<<r<<'\n';


}