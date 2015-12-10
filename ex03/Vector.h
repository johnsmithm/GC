
#include <algorithm>
#include <cassert>
#include <iostream>
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
	 Vector(size_t nx,size_t ny,  std::function<double(size_t,size_t)> f)  :n(nx*ny)   
    {
		 nx=ny=0;	
	   data = new double[n];	 
       for(size_t i=0;i<ny;++i)
		   for(size_t j=0;j<nx;++j)
		        data[j+nx*i] = f(i,j);
    }
	Vector( int n_,size_t beg, size_t end, std::function<double(size_t)> f) :n(n_)
    {   
	   nx = end - beg;
	   ny = n/nx;	
	   data = new double[n];	
	   for(int i=0;i<n;++i)data[i]=0;
       for(size_t i=beg;i<end;++i)data[i] = f(i);
    }
	~Vector(){delete [] data;}
	double operator[] (int i)const{
	return data[i];
	}
	double& operator[] (int i){
	return data[i];
	}
	template<class A>
	void operator = (const Expr<A>& a_){
	const A& a(a_);
		for(int i=0;i<n;++i){
		data[i] = a[i];
		}
	}
	double LNorm(){
		double l = 0;		
		for(int i=0;i<n;++i)l+=data[i]*data[i];
		return l;
	}
	double operator^(const Vector & a)const{
		double l = 0;		
		for(int i=0;i<n;++i)l+=data[i]*a[i];
		return l;
	}
    size_t size()const{return n;}
	size_t nx_()const{return nx;}
	size_t ny_()const{return ny;}
};


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

class Stencil : public Expr<Vector > {


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
	Vector operator * (Vector& u){
		 Vector r(pg,u.nx_(),u.ny_());
		 for(int i=1;i<ny_-1;++i)
					for(int j=1;j<nx_-1;++j)
                    	r[i*nx_+j] = mst*u[i*nx_+j] - xst*(u[i*nx_+j+1]+u[i*nx_+j-1])-yst*(u[j+(i+1)*nx_]+u[j+(i-1)*nx_]);
		 return r;
	}
	private :
		int nx_,ny_,pg;
		double yst, xst, mst;
		const double pi;
		
};


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

