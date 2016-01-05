
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <math.h>
#include <functional>

const int UP    = 0;
const int DOWN  = 1;
const int LEFT  = 2;
const int RIGHT = 3;
MPI_Comm cartcomm( MPI_COMM_NULL );


template<class A>
	struct  Expr {
    operator const A& () const{
	return *static_cast<const A*>(this);
	}
};




class Vector : public Expr<Vector > {

public :
	void commit_datatypes(){
	   MPI_Type_vector( nyi, 1, nxs, MPI_DOUBLE, &columntype );
	   MPI_Type_commit( &columntype );
	
	   MPI_Type_vector( ny, nx, nxs, MPI_DOUBLE, &printtype );
	   MPI_Type_commit( &printtype );	 
	}
	
	Vector(int n_,int nx_,int ny_,int nxi_,int nyi_,int nxs_,int nys_,int x_,int y_, double w = 0)//d,r,z
		:n(n_),x(x_),y(y_),nx(nx_),ny(ny_),nxs(nxs_),nys(nys_),nxi(nxi_),nyi(nyi_){	  	   		
       data = new double[n];
	   for(size_t i=0;i<n;++i)data[i]=w;
	   commit_datatypes();		
	}	
	
	 Vector(int n_,int nx_,int ny_,int nxi_,int nyi_,int nxs_,int nys_,int x_,int y_,
			size_t offsetx1,size_t offsety1,std::function<double(size_t,size_t)> f) //f
		 :n(n_),x(x_),y(y_),nx(nx_),ny(ny_),nxs(nxs_),nys(nys_),nxi(nxi_),nyi(nyi_)
	{	
	  data = new double[n];	 
	  for(size_t i=0;i<n;++i)data[i]=0;		 
	  for(size_t i=0;i<ny;++i)
		      for(size_t j=0;j<nx;++j)
			    data[j+x+nxs*(i+y)] = f(j+offsetx1,i+offsety1);
	  commit_datatypes();	  
	}
	
	    Vector(int n_,int nx_,int ny_,int nxi_,int nyi_,int nxs_,int nys_,int x_,int y_,size_t offsetx1, std::function<double(size_t)> f)//u
			:n(n_),x(x_),y(y_),nx(nx_),ny(ny_),nxs(nxs_),nys(nys_),nxi(nxi_),nyi(nyi_)
	{   	
	      data = new double[n];	
	      for(size_t i=0;i<n;++i)data[i]=0;	
				
	      for(size_t i=0;i<nx;++i)
			  data[nxs*(nys-1)+x+i] = f(i+offsetx1);
		  commit_datatypes();
	}//end Vector Constructors
    //Vector(Vector&& o) noexcept : data(std::move(o.data)) {std::cout<<'|';}
	~Vector(){
		delete [] data;
		MPI_Type_free( &columntype );
		MPI_Type_free( &printtype );
	}
	
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
		
		for(size_t i=0;i<nyi;++i)
			for(size_t j=0;j<nxi;++j)
		{
		data[j+1+nxs*(1+i)] = a[(1+i)*nxs+j+1];
		}
	}
	
	void operator = (const Vector& a)
	//:n(a.n),x(a.x),y(a.y),nx(a.nx),ny(a.ny),nxs(a.nxs),nys(a.nys)
	{
		n=a.n;x=a.x;y=a.y;nx=a.nx;ny=a.ny;nxs=a.nxs;nys=a.nys;
			for(size_t i=0;i<n;++i){
			data[i] = a[i];
			}
	}
	
	void operator = ( Vector&& a)
	//:n(a.n),x(a.x),y(a.y),nx(a.nx),ny(a.ny),nxs(a.nxs),nys(a.nys)
	{//move constructor
		   // std::cout<<"-";
		    n=a.n;x=a.x;y=a.y;nx=a.nx;ny=a.ny;nxs=a.nxs;nys=a.nys;
		    nxi=a.nxi;nyi=a.nyi;
		    offsetx=a.offsetx;offsety=a.offsety;
			data = a.data;
			a.data = NULL;
	}
	
	double operator^(const Vector & a)const{
		double l = 0;		
		for(size_t i=0;i<nyi;++i)
			for(size_t j=0;j<nxi;++j)			
			      l+=data[(1+i)*nxs+j+1]*a.data[(1+i)*nxs+j+1];
		double l1;		
	    MPI_Allreduce( &l, &l1 ,1 , MPI_DOUBLE ,MPI_SUM, cartcomm );
		return l1;
	}
	//end vector operation
	//brgin vector members
	double LNorm(){
		double l = 0;		
		for(size_t i=0;i<nyi;++i)
			for(size_t j=0;j<nxi;++j)
			    l+=data[(1+i)*nxs+j+1]*data[j+1+nxs*(i+1)];
		double l1;		
	    MPI_Allreduce( &l, &l1 ,1 , MPI_DOUBLE ,MPI_SUM, cartcomm );	
		return l1;
	}	
	
	void get_info(int rank){
	               MPI_Status status;
				   MPI_Recv( &data[nxs+x], 1, printtype, rank, 10, cartcomm, &status );
	}
	void sent_info(){
        
	               MPI_Send( &data[nxs+x],    1, printtype, 0, 10, cartcomm );
	}
	
	//template <>
	//friend Add<Stencil,Vector>::Add(Stencil,Vector);
	void  set_distace(double const x1,double const y1){hx=x1;hy=y1; }
	void  set_offset(int const x1,int const y1)       {offsetx=x1;offsety=y1; }
	
    size_t size()const{return n;}
	size_t offsetx__()const{return offsetx;}
	size_t offsety__()const{return offsety;}
	double hx__() const{return hx;}
	double hy__() const{return hy;}
	size_t nx__() const{return nx;}
	size_t ny__() const{return ny;}
	size_t nxs__()const{return nxs;}
	size_t nys__()const{return nys;}
	size_t x__()  const{return x;}
	size_t y__()  const{return y;}
	//end vector members
    public:
    //:n(n_),x(x_),y(y_),nx(nx_),ny(ny_),nxs(nxs_),nys(nys_),nxi(nxi_),nyi(nyi_)
	
	size_t n;
    size_t x;//first  point for f initialization
	size_t y;
	size_t nx;
    size_t ny;//inner points for f
	size_t nxs;
    size_t nys;//total points
	size_t nxi;
    size_t nyi;//inner points 
    double *data;
	double hx;
    double hy;//for print
	size_t offsetx;
    size_t offsety;
	MPI_Datatype columntype; 
	MPI_Datatype printtype;
};

//vector print operators
std::ostream & operator<<( std::ostream & os, const Vector & v )
{

	 for( size_t i=0; i < v.nys__(); ++i )
		{
			for( size_t j=0; j < v.nxs__(); ++j )
			{
				os << v[i*v.nxs__()+j] << " ";        
			}  
		    os<<"\n";
		}
	
	os<<"\n";
    return os;
}

std::ostream & operator&( std::ostream & os, const Vector & v )
{
	
	 for( size_t i=0; i < v.ny__(); ++i )
		{
			for( size_t j=0; j < v.nx__(); ++j )
			{  // os<<j<<" "<<i<<" ";
				os<<(j+v.offsetx__())*v.hx__()<<' '<<(i+v.offsety__())*v.hy__()<<' ' << v[(v.y__()+i)*v.nxs__()+j+v.x__()] << "\n";        
			}  
		    //MPI_Barrier( cartcomm );//??
		    //os<<"\n";
		}
	
	os<<"\n";
    return os;
}//vector end



class Stencil : public Expr<Stencil > {


public :
	Stencil(int nx1,int ny1, int * nbrs_):nx_(nx1+1),ny_(ny1+1),pi(3.141592653589793),nbrs(nbrs_){
	        double hx_ = 2.0/nx1;
            double hy_ = 1.0/ny1;
            xst =  1.0/(hx_*hx_);
            yst =  1.0/(hy_*hy_);
            mst = 2.0/(hx_*hx_)+2.0/(hy_*hy_)+4*pi*pi;          
            pg = (nx1+1)*(ny1+1);  
	}
	~Stencil(){}	
	/*
	Vector operator - (Vector& u){
		//std::cout<<"--";
		 Vector r(pg,u.nx__(),u.ny__());
		 for(int i=1;i<ny_-1;++i)
					for(int j=1;j<nx_-1;++j)
                    	r[i*nx_+j] = mst*u[i*nx_+j] - xst*(u[i*nx_+j+1]+u[i*nx_+j-1])-yst*(u[j+(i+1)*nx_]+u[j+(i-1)*nx_]);
		 return r;
	}*/
	
	
	
	int nx() const{return nx_;}
	int ny() const{return ny_;}
	
	int left() const{return nbrs[LEFT];}
	int right()const{return nbrs[RIGHT];}
	int up()   const{return nbrs[UP];}
	int down() const{return nbrs[DOWN];}
	
	double yst_()const{return yst;}
	double xst_()const{return xst;}
	double mst_()const{return mst;}
	private :
		int nx_,ny_,pg;
		double yst, xst, mst;
		const double pi;
		int * nbrs;
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
		int x__()const{return b_.x__();}
		int y__()const{return a_.y__();}
		int nxs__()const{return b_.nxs__();}
	};

template <class A, class B>
	class Minus : public Expr<Minus<A,B>> {
	const A& a_;
	const B& b_;
	public :
		Minus(const A& a,const B& b): a_(a),b_(b){}
		
		int x__()const{return b_.x__();}
		int y__()const{return a_.y__();}
		int nxs__()const{return b_.nxs__();}
		
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
		int x__()const{return a_.x__();}
		int y__()const{return a_.y__();}
		int nxs__()const{return a_.nxs__();}
	};

template <>
	class Add<Stencil,Vector> : public Expr<Add<Stencil,Vector>> {
	const Vector& a_;
	const Stencil& b_;
    MPI_Request reqs[8];
    MPI_Status stats[8];	
	public :
		Add(const Vector& a,const Stencil& b): a_(a),b_(b){
			//TODO check index and add if to nx ny making
			MPI_Isend( &a.data[1+a.nxs],                a.nxi, MPI_DOUBLE,   b.up()  ,  0, cartcomm, &reqs[0]   ); // first row
			MPI_Isend( &a.data[1+a.nxs*(a.nys-2)],      a.nxi, MPI_DOUBLE,   b.down(),  1, cartcomm, &reqs[1]   ); //last row
			MPI_Isend( &a.data[1+a.nxs],                1,    a.columntype, b.left() , 2, cartcomm, &reqs[2]   ); //first column
			MPI_Isend( &a.data[1+a.nx+a.nxs],           1,    a.columntype, b.right(), 3, cartcomm, &reqs[3]   );//last column

			MPI_Irecv( &a.data[1],                     a.nxi, MPI_DOUBLE,   b.up(),    1, cartcomm, &reqs[4]   );//first row
			MPI_Irecv( &a.data[1+a.nxs*(a.nys-1)],     a.nxi, MPI_DOUBLE,   b.down(),  0, cartcomm, &reqs[5]   );//last row
			MPI_Irecv( &a.data[a.nxs*2-1],             1,    a.columntype, b.right(), 2, cartcomm, &reqs[6]   );//last column
			MPI_Irecv( &a.data[a.nxs],                 1,    a.columntype, b.left(),  3, cartcomm, &reqs[7]   );//first column


            MPI_Waitall( 8, reqs, stats );
		}
		double operator[](int i)const{
		return (b_.mst_()*a_[i] - b_.xst_()*(a_[i+1]+a_[i-1])-b_.yst_()*(a_[i+a_.nxs__()]+a_[i-a_.nxs__()]));
		}
		int x__()const{return a_.x__();}
		int y__()const{return a_.y__();}
		int nxs__()const{return a_.nxs__();}
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

Vector * get_vector(int rank,int nx,int ny,int px = 2,int py = 2){
	    
	    int coords[2] = {0,0};
        MPI_Cart_coords(cartcomm, rank, 2, coords);
	    std::swap(coords[0],coords[1]);
	
	    int nx_ = (nx+1)/px , ny_ = (ny+1)/py; //inner point in each block for f
		int nxi = nx_, nyi = ny_;//inner point for stencil calculation
		int nxs = nx_ , nys = ny_;//total size of a block
		int x=0,y=0; //first inner point coordonates
		int offsetx=coords[1]*nx_, offsety =coords[0]*ny_;
	

		if(coords[0]==0 || coords[0]==py-1){//if we are last or first row
			if(coords[0]==py-1)
				ny_ += (ny+1)%py; // add the remaining rows
			nys = ny_+1;
		}else nys +=2; //add the layers

		if(coords[1]==0 || coords[1]==px-1){//if we are last or first column
			if(coords[1]==px-1)
				nx_ += (nx+1)%px; // add the remaining columns
			nxs = nx_+1;
		}else nxs +=2; //add the layers


		if(coords[0]==0 && coords[1]==0){//first point
			x=0;y=0;		
		}else if(coords[0]==0){//first column
		   x=1;y=0;
		}else if(coords[1]==0)//first row
		{
		x=0; y=1;
		}else{
			x=1;y=1;
		}//others

		nxi = nx_;nyi = ny_;
		if(coords[0]==0 || coords[0]==py-1){//first row or last row
			--nyi;		
		}
		if(coords[1]==0 || coords[1]==px-1){//first column
			--nxi;
		}
        int pg = nxs*nys;
		Vector * r = new Vector(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y);	
		r->set_offset(offsetx,offsety);
		return r;
}

int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89};

int getHeight(int size){
    int height = 1;
    bool use = true;
    for(int i = 20;i >= 0;--i){
        while(size%primes[i] == 0 && size != 1){
            size /= primes[i];
            use = !use;
            if(use) {
                height *= primes[i];
            }
        }
    }
    return height;
}

double Expr_CG(int nx,int ny,int c,double eps, int rank, int nrpr){
	//MPI topology      
	  int px = getHeight(nrpr);
      int py = nrpr/px;
      if(rank == 0){
        std::cout<<"Height="<<px<<", Width="<<py<<";\n";
      }
	  int dims[2] = {px,py};
	  int periods[2] = {0,0}; // not periodic!
      const int reorder = 1;  // allow reordering of process ranks      
      int cartrank(0);
      int coords[2] = {0,0};
      int nbrs[4] = {0,0,0,0};
      MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm );
      MPI_Comm_rank( cartcomm, &cartrank );
      MPI_Cart_coords( cartcomm, cartrank, 2, coords );
      MPI_Cart_shift( cartcomm, 0, 1, &nbrs[LEFT], &nbrs[RIGHT] );
      MPI_Cart_shift( cartcomm, 1, 1, &nbrs[UP], &nbrs[DOWN] );
	  std::swap(coords[0],coords[1]);
	
 
	
	
	//initialization	
	int nx_ = (nx+1)/px , ny_ = (ny+1)/py; //inner point in each block for f
	int nxi = nx_, nyi = ny_;//inner point for stencil calculation
	int nxs = nx_ , nys = ny_;//total size of a block
	int x=0;
    int y=0; //first inner point coordonates
	int offsetx=coords[1]*nx_, offsety =coords[0]*ny_;
	
	if(coords[0]==0 || coords[0]==py-1){//if we are last or first row
	    if(coords[0]==py-1)
			ny_ += (ny+1)%py; // add the remaining rows
		nys = ny_+1;
	}else nys +=2; //add the layers
	
	if(coords[1]==0 || coords[1]==px-1){//if we are last or first column
	    if(coords[1]==px-1)
			nx_ += (nx+1)%px; // add the remaining columns
		nxs = nx_+1;
	}else nxs +=2; //add the layers
	
	
	if(coords[0]==0 && coords[1]==0){//first point
	    x=0;y=0;		
	}else if(coords[0]==0){//first column
	   x=1;y=0;
	}else if(coords[1]==0)//first row
	{
	x=0; y=1;
	}else{
		x=1;y=1;
	}//others
	
	nxi = nx_;nyi = ny_;
	if(coords[0]==0 || coords[0]==py-1){//first row or last row
	    --nyi;		
	}
	if(coords[1]==0 || coords[1]==px-1){//first column
	    --nxi;
	}
	
	
	
	int pg = nxs*nys;	
	//int n_,int nx_,int ny_,int  nxi,int nyi,int nxs_,int nys_,int x_,int y_,size_t offsetx,size_t offsety
	Vector  r(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y),
			d(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y),
			z(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y); 
	double pi = 3.141592653589793;
	double hx_ = 2.0/nx;
    double hy_ = 1.0/ny;
	double C = 4*pi*pi;   
	double freqx = 2*pi*hx_;   
	double freqy = 2*pi*hy_;     
	//4π^2 sin(2πx) sinh(2πy)
	Vector f(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y,offsetx,offsety,
			 [C,freqx,freqy](int x1,int y1)->double{return C*sin(freqx*x1)*sinh(freqy*y1);});	//set freq
			double SINH = sinh(2*pi); 
	
	//sin(2πx) sinh(2πy)
	Vector *ut;
	if(coords[0]==py-1)//we are at border's block
	     ut= new Vector(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y,offsetx, [SINH,freqx](int x1)->double{return sin(x1*freqx) * SINH;});
	else ut = new Vector(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y);
	Vector u(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y);
	u = *ut;
	
	Stencil A(nx,ny,nbrs);
	double delta0 = 0, delta1 = 0, beta = 0,alfa=0 ,scalar =0;
	//initialization end
	
	

	
	//CG
	r = f - A*u;
	delta0 = r.LNorm();
	d = r;
	if(sqrt(delta0)>eps)
	  for(int i=0;i<c;++i){
		  z = A*d;		  
		  scalar = (d^z);
		  alfa = delta0 / scalar;
		  u = u + d*alfa;
		  r = r - z*alfa;
		  delta1 = r.LNorm();
		  beta = delta1/delta0;
		  delta0=delta1;
		  if(sqrt(delta1)<eps)break;		  
		  d = r + d*beta;		  
	  }
	
		/*for( int i = 0; i < px*py; ++i )
      {
         if( cartrank == i )
         {

            std::cout<<u << "------------------------------------------------------------\n"
			          << "rank (Cartesian topology):         " << cartrank << "\n"
                      << "Cartesian coordinates:             ( " << coords[0] << ", " << coords[1] << " )\n" 
		              << "neighbors (x-direction, expected): " << nbrs[LEFT] << " (left), " << nbrs[RIGHT] << " (right)\n"
                      << "neighbors (y-direction, expected): " << nbrs[DOWN] << " (down), " << nbrs[UP] << " (up)\n"   
				      << "nx="<<nx_<<" ny="<<ny_<<"\n"
				      << "nxs="<<nxs<<" nys="<<nys<<"\n"
				      << "nxi="<<nxi<<" nyi="<<nyi<<"\n"					      
				      << "x="<<x<<" y="<<y<<"\n"
                      << std::endl;
         }
         MPI_Barrier( MPI_COMM_WORLD );
      }*/
	
	//print
	std::ofstream out("solution.txt");
	for( int i = 0; i < px*py; ++i )
		  {
			 if( cartrank == 0 )
			 {
				if( i == 0 )
				{
				   u.set_distace(hx_,hy_);
				   u.set_offset(offsetx,offsety);	
				   out&u;
				}
				else
				{
				   //out<<i<<'\n';	
				   Vector * pr = get_vector(i,nx,ny,px,py);	
				   pr->get_info(i);	
				   pr->set_distace(hx_,hy_);		
				   out&(*pr);	
				  
					
				}
			 }
			 else if( cartrank == i )
			 {
				   u.sent_info();	
			 }
		  }
	
	return sqrt(delta0);	
}
