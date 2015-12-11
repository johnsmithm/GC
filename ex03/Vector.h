
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <math.h>
#include <functional>


template<class A>
	struct  Expr {
    operator const A& () const{
	return *static_cast<const A*>(this);
	}
};

class Vector : public Expr<Vector > {
private:
	double *data;
	int nx,ny;
	int n;
public :
	Vector(int n_,double w = 0):n(n_){
	   nx=ny=0;	
       data = new double[n];
	   for(int i=0;i<n;++i)data[i]=w;
	}
	
	Vector(int n_,int nx_,int ny_, double w = 0):n(n_){
	   nx=nx_;
	   ny=ny_;	
       data = new double[n];
	   for(int i=0;i<n;++i)data[i]=w;
	}
	
	 Vector(size_t nx_,size_t ny_,  std::function<double(size_t,size_t)> f)  :n(nx_*ny_)   
	{
		    nx=nx_;
			ny=ny_;	
	      data = new double[n];	 
	  for(size_t i=0;i<ny_;++i)
		      for(size_t j=0;j<nx_;++j)
			    data[j+nx_*i] = f(j,i);
	}
	
	    Vector( int n_,size_t beg, size_t end, std::function<double(size_t)> f) :n(n_)
	{   
	      nx = end - beg;
	      ny = n/nx;	
	      data = new double[n];	
	      for(int i=0;i<n;++i)data[i]=0;
	      for(size_t i=beg;i<end;++i)data[i] = f(i);
	}//end Vector Constructors
    Vector(Vector&& o) noexcept : data(std::move(o.data)) {std::cout<<'|';}
	~Vector(){delete [] data;}
	
	//vector operators
	double operator[] (int i)const{
	return data[i];
	}
	
	double& operator[] (int i){
	return data[i];
	}
	
	template<class A>
	void operator = (const Expr<A>& a_){
	const A& a(a_);
	        //nx  = a.nx_();
	        //ny = a.ny_();
		for(int i=1;i<ny-1;++i)
			for(int j=1;j<nx-1;++j)
		{
		data[j+nx*i] = a[j+nx*i];
		}
	}
	
	void operator = (const Vector& a){
	        nx  = a.nx_();
	        ny = a.ny_();
			n = a.size();
			for(int i=0;i<n;++i){
			data[i] = a[i];
			}
	}
	void operator = ( Vector&& a){//move constructor
		   // std::cout<<"-";
	        nx  = a.nx_();
	        ny = a.ny_();
			n = a.size();
			data = a.data;
			a.data = NULL;
	}
	
	double operator^(const Vector & a)const{
		double l = 0;		
		for(int i=0;i<n;++i)l+=data[i]*a.data[i];
		return l;
	}
	//end vector operation
	//brgin vector members
	double LNorm(){
		double l = 0;		
		for(int i=0;i<n;++i)l+=data[i]*data[i];
		return l;
	}	
    size_t size()const{return n;}
	size_t nx_()const{return nx;}
	size_t ny_()const{return ny;}
	//end vector members
};

//vector print operators
std::ostream & operator<<( std::ostream & os, const Vector & v )
{
	if(v.nx_()!=0){
	 for( size_t i=0; i < v.ny_(); ++i )
		{
			for( size_t j=0; j < v.nx_(); ++j )
			{
				os << v[i*v.nx_()+j] << " ";        
			}  
		    os<<"\n";
		}
	}else
    for( size_t i=0; i < v.size(); ++i )
    {
        os << v[i] << " ";        
    }
	os<<"\n";
    return os;
}

std::ostream & operator&( std::ostream & os, const Vector & v )
{
	if(v.nx_()!=0){
	 for( size_t i=0; i < v.ny_(); ++i )
		{
			for( size_t j=0; j < v.nx_(); ++j )
			{
				os<<i<<' '<<j<<' ' << v[i*v.nx_()+j] << "\n";        
			}  
		    //os<<"\n";
		}
	}else
    for( size_t i=0; i < v.size(); ++i )
    {
        os << v[i] << " ";        
    }
	os<<"\n";
    return os;
}//vector end

class Stencil : public Expr<Stencil > {


public :
	Stencil(int nx,int ny):nx_(nx+1),ny_(ny+1),pi(3.141592653589793){
	    double hx_ = 2.0/nx;
            double hy_ = 1.0/ny;
            xst =  1.0/(hx_*hx_);
            yst =  1.0/(hy_*hy_);
            mst = 2.0/(hx_*hx_)+2.0/(hy_*hy_)+4*pi*pi;          
            pg = (nx+1)*(ny+1);  
	}
	~Stencil(){}	
	
	Vector operator - (Vector& u){
		//std::cout<<"--";
		 Vector r(pg,u.nx_(),u.ny_());
		 for(int i=1;i<ny_-1;++i)
					for(int j=1;j<nx_-1;++j)
                    	r[i*nx_+j] = mst*u[i*nx_+j] - xst*(u[i*nx_+j+1]+u[i*nx_+j-1])-yst*(u[j+(i+1)*nx_]+u[j+(i-1)*nx_]);
		 return r;
	}
	int nx() const{return nx_;}
	int ny() const{return ny_;}
	double yst_()const{return yst;}
	double xst_()const{return xst;}
	double mst_()const{return mst;}
	private :
		int nx_,ny_,pg;
		double yst, xst, mst;
		const double pi;
		
};//stencil end

//operation begin
template <class A, class B>
	class Add : public Expr<Add<A,B>> {
	const A& a_;
	const B& b_;
	public :
		Add(const A& a,const B& b): a_(a),b_(b){}
		double operator[](int i)const{
		return (a_[i] + b_[i]);
		}
	};

template <class A, class B>
	class Minus : public Expr<Minus<A,B>> {
	const A& a_;
	const B& b_;
	public :
		Minus(const A& a,const B& b): a_(a),b_(b){}
		double operator[](int i)const{
		return (a_[i] - b_[i]);
		}
	};

template <class A>
	class Add<A,double> : public Expr<Add<A,double>> {
	const A& a_;
	const double& b_;
	public :
		Add(const A& a,const double& b): a_(a),b_(b){}
		double operator[](int i)const{
		return (a_[i] * b_);
		}
	};

template <class A>
	class Add<Stencil,A> : public Expr<Add<Stencil,A>> {
	const A& a_;
	const Stencil& b_;
	public :
		Add(const A& a,const Stencil& b): a_(a),b_(b){}
		
		double operator[](int i)const{
		return (b_.mst_()*a_[i] - b_.xst_()*(a_[i+1]+a_[i-1])-b_.yst_()*(a_[i+b_.nx()]+a_[i-b_.nx()]));
		}
	};



template <class A, class B>
	inline Add<A,B> operator+ (const Expr<A>& a, const Expr<B>& b){
	return Add<A,B>(a,b);
	}

template <class A, class B>
	inline Minus<A,B> operator- (const Expr<A>& a, const Expr<B>& b){
	return Minus<A,B>(a,b);
	}

template <class A>
	inline Add<A,double> operator* (const Expr<A>& a, const double& b){
	return Add<A,double>(a,b);
	}

//that can be slower cause of checking
template <class A>
	inline Add<Stencil,A> operator* ( const Stencil& b, const Expr<A>& a){
	return Add<Stencil,A>(a,b);
	}
//operation end


double Expr_CG(int nx,int ny,int c,double eps){
	//initialization	
	int pg = (1+nx)*(ny+1);	
	Vector r(pg,1+nx,1+ny), d(pg,1+nx,1+ny), z(pg,1+nx,1+ny); 
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
	//std::cout<<(A*u);//u<<(A*u);
	delta0 = r.LNorm();
	d = r;
	if(sqrt(delta0)>eps)
	  for(int i=0;i<c;++i){
		  z = A*d;		  
		  alfa = delta0 / (d^z);
		  u = u + d*alfa;
		  r = r - z*alfa;
		  delta1 = r.LNorm();
		  beta = delta1/delta0;
		  delta0=delta1;
		  if(sqrt(delta1)<eps)break;		  
		  d = r + d*beta;		  
	  }
	//std::cout<<u;
	std::ofstream out("solution.txt");
	out&u;
	return sqrt(delta0);	
}
