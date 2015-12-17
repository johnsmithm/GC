#include<math.h>
#include<iostream>
#include<fstream>

class CG{
    public:
        CG(int nx,int ny,int c,double eps):eps_(eps),nx_(nx+1),ny_(ny+1),c_(c),pi(3.141592653589793){
            //distance between grid points    
            hx_ = 2.0/nx;
            hy_ = 1.0/ny;
            //stencil    
            xst =  1.0/(hx_*hx_);
            yst =  1.0/(hy_*hy_);
            mst = 2.0/(hx_*hx_)+2.0/(hy_*hy_)+4*pi*pi;          
            pg = nx_*ny_;           

        }

        ~CG(){}   
    private:

	void matrix_multiplication(double r[],double u[]){
	  for(int i=1;i<ny_-1;++i)
                for(int j=1;j<nx_-1;++j)
                    r[i*nx_+j] = mst*u[i*nx_+j] - xst*(u[i*nx_+j+1]+u[i*nx_+j-1])-yst*(u[j+(i+1)*nx_]+u[j+(i-1)*nx_]);   
	}
	

	void vector_addition(double r[],double a[],double b[],double scalar){
	    for(int i=1;i<ny_-1;++i)
                for(int j=1;j<nx_-1;++j)
			           r[i*nx_+j]=a[i*nx_+j]+b[i*nx_+j]*scalar;
	}
	
	void vector_subtraction(double r[],double a[],double b[],double scalar){
	    for(int i=1;i<ny_-1;++i)
                for(int j=1;j<nx_-1;++j)
			          r[i*nx_+j]=a[i*nx_+j]-b[i*nx_+j]*scalar;
	}
	
	double LNorm(double a[],double b[]){
		double rs=0;
	    for(int i=1;i<ny_-1;++i)
                for(int j=1;j<nx_-1;++j)
			         rs+=a[i*nx_+j]*b[i*nx_+j];
		return rs;
	}
	
    public:  
    double solve_naive(){
       initialization_naive();
      
     
		matrix_multiplication(z,u);
		// view(f);
		vector_subtraction(d,f,z,1);		
		vector_subtraction(r,f,z,1);
		delta0 = LNorm(d,d);
		if(sqrt(delta0)>eps_)
		   for(int it=0;it<c_;++it){//nr iterations
			   matrix_multiplication(z,d);
			   alfa= delta0/LNorm(d,z);
			   vector_addition(u,u,d,alfa);
			   vector_subtraction(r,r,z,alfa);
			   delta1 = LNorm(r,r);			   
			   beta = delta1/delta0;
			   delta0=delta1; 
			   if(sqrt(delta1)<=eps_)break;
			   vector_addition(d,r,d,beta);		   
			   //if(it%20==0||1)std::cout<<delta1<<" -\n";
		   }
        //view(u);
       return  sqrt(delta0);
    }
    private:
    //grid points for right to left, bottom to top
    void initialization_naive(){    
       //using one vector for grid points    
       f = new double[pg];    
       u = new double[pg];    
       r = new double[pg];    
       d = new double[pg]; 		   
       z = new double[pg]; 	    
        //TODO OpenMP.
       double C = 4*pi*pi;   
       double freqx = 2*pi*hx_;   
       double freqy = 2*pi*hy_;     
        
       //TODO  threads
       for(int i=0;i<ny_;++i)
          for(int j=0;j<nx_;++j){
           f[j+i*nx_]=C*sin(freqx*j)*sinh(freqy*i);//4π^2 sin(2πx) sinh(2πy)
           u[j+i*nx_]=0;
		   r[j+i*nx_]=0;
		   d[j+i*nx_]=0;
		   z[j+i*nx_]=0;
        }       
        
       int last_row = (ny_-1)*nx_;  
       double SINH = sinh(2*pi);  
       for(int i=0;i<nx_;++i){           
           u[last_row+i] = sin(i*freqx) * SINH;//sin(2πx) sinh(2πy)
       }
    }
    
    void view(double *a){
        for(int i=0;i<ny_;i++){
            for(int j=0;j<nx_;j++)std::cout<<a[i*nx_+j]<<" ";
            std::cout<<"\n";
        }
    }
 
    
    
    
 
    
    public: void print_gnuplot_naive(){
    
        std::ofstream out("solution.txt");
        
        for(int i=0;i<ny_;++i){//every red grid point
                for(int j=0;j<nx_;j+=1){                   
                    
                    out<<(j*hy_)<<" "<<(i*hx_)<<" "<<u[i*nx_+j]<<"\n";
                }
        }       
    }
   
  

    private:
            //we make nx_ even
            bool even;
            //optimize implemantation: red point, black points
            double * u ,*f , *r, *z, *d, beta,alfa, delta0, delta1,eps_;
            //delta x, delta y, left and right points stencil, top and bottom point stencil, middle point stencil
            double hx_,hy_, xst, yst, mst;
            //nr of iteration, nr of grid points x ax, nr of grid points y ax, #points grid, # threads
            int nx_,ny_, c_, pg;
            //pi from maths
            const double pi;
};
