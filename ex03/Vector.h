
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
	//virtual int x_() const= 0;
	//virtual int nxs_()const = 0;
};

/*template <class A, class B>
	class Add : public Expr<Add<A,B>>;
template <>
	class Add<Stencil,Vector> : public Expr<Add<Stencil,Vector>> ;
*/


class Vector : public Expr<Vector > {
public:
	double *data;
	size_t nx,ny;//inner points for f
	size_t nxi,nyi;//inner points 
	size_t nxs,nys;//total points
	size_t n;
	size_t x,y;//first  point for f initialization
	MPI_Datatype columntype;   	   
public :
	void commit_datatypes(){
	   MPI_Type_vector( ny, 1, nys, MPI_DOUBLE, &columntype );
	   MPI_Type_commit( &columntype );
	}
	/*Vector(int n_,double w = 0):n(n_){
	   nx=ny=0;	
       data = new double[n];
	   for(int i=0;i<n;++i)data[i]=w;
	}*/
	
	Vector(int n_,int nx_,int ny_,int nxi_,int nyi_,int nxs_,int nys_,int x_,int y_, double w = 0)//d,r,z
		:n(n_),x(x_),y(y_),nx(nx_),ny(ny_),nxs(nxs_),nys(nys_),nxi(nxi_),nyi(nyi_){	  	   		
       data = new double[n];
	   for(int i=0;i<n;++i)data[i]=w;
	   commit_datatypes();		
	}	
	
	 Vector(int n_,int nx_,int ny_,int nxi_,int nyi_,int nxs_,int nys_,int x_,int y_,
			size_t offsetx,size_t offsety,std::function<double(size_t,size_t)> f) //f
		 :n(n_),x(x_),y(y_),nx(nx_),ny(ny_),nxs(nxs_),nys(nys_),nxi(nxi_),nyi(nyi_)
	{	
	  data = new double[n];	 
	  for(size_t i=0;i<ny;++i)
		      for(size_t j=0;j<nx;++j)
			    data[j+x+nxs*(i+y)] = f(j+offsetx,i+offsety);
	  commit_datatypes();	  
	}
	
	    Vector(int n_,int nx_,int ny_,int nxi_,int nyi_,int nxs_,int nys_,int x_,int y_,size_t offsetx, std::function<double(size_t)> f)//u
			:n(n_),x(x_),y(y_),nx(nx_),ny(ny_),nxs(nxs_),nys(nys_),nxi(nxi_),nyi(nyi_)
	{   	
	      data = new double[n];	
	      for(int i=0;i<n;++i)data[i]=0;	
				
	      for(size_t i=0;i<nx;++i)
			  data[nxs*(nys-1)+x+i] = f(i+offsetx);
		  commit_datatypes();
	}//end Vector Constructors
    //Vector(Vector&& o) noexcept : data(std::move(o.data)) {std::cout<<'|';}
	~Vector(){delete [] data; MPI_Type_free( &columntype );}
	
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
		
		for(int i=0;i<nyi;++i)
			for(int j=0;j<nxi;++j)
		{
		data[j+1+nxs*(1+i)] = a[(1+i)*nxs+j+1];
		}
	}
	
	void operator = (const Vector& a)
	//:n(a.n),x(a.x),y(a.y),nx(a.nx),ny(a.ny),nxs(a.nxs),nys(a.nys)
	{
		n=a.n;x=a.x;y=a.y;nx=a.nx;ny=a.ny;nxs=a.nxs;nys=a.nys;
			for(int i=0;i<n;++i){
			data[i] = a[i];
			}
	}
	
	void operator = ( Vector&& a)
	//:n(a.n),x(a.x),y(a.y),nx(a.nx),ny(a.ny),nxs(a.nxs),nys(a.nys)
	{//move constructor
		   // std::cout<<"-";
		    n=a.n;x=a.x;y=a.y;nx=a.nx;ny=a.ny;nxs=a.nxs;nys=a.nys;
			data = a.data;
			a.data = NULL;
	}
	
	double operator^(const Vector & a)const{
		double l = 0;		
		for(int i=0;i<nyi;++i)
			for(int j=0;j<nxi;++j)			
			      l+=data[(1+i)*nxs+j+1]*a.data[(1+i)*nxs+j+1];
		double l1;		
	    MPI_Allreduce( &l, &l1 ,1 , MPI_DOUBLE ,MPI_SUM, cartcomm );
		return l1;
	}
	//end vector operation
	//brgin vector members
	double LNorm(){
		double l = 0;		
		for(int i=0;i<nyi;++i)
			for(int j=0;j<nxi;++j)
			    l+=data[(1+i)*nxs+j+1]*data[j+1+nxs*(i+1)];
		double l1;		
	    MPI_Allreduce( &l, &l1 ,1 , MPI_DOUBLE ,MPI_SUM, cartcomm );	
		return l1;
	}	
	
	//template <>
	//friend Add<Stencil,Vector>::Add(Stencil,Vector);
	
    size_t size()const{return n;}
	size_t nx_()const{return nx;}
	size_t ny_()const{return ny;}
	size_t nxi_()const{return nxi;}
	size_t nyi_()const{return nyi;}
	size_t nxs_()const{return nxs;}
	size_t nys_()const{return nys;}
	size_t x_()const{return x;}
	size_t y_()const{return y;}
	//end vector members
};

//vector print operators
std::ostream & operator<<( std::ostream & os, const Vector & v )
{
	if(v.nxs_()!=0){
	 for( size_t i=0; i < v.nys_(); ++i )
		{
			for( size_t j=0; j < v.nxs_(); ++j )
			{
				os << v[i*v.nxs_()+j] << " ";        
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
	 for( size_t i=0; i < v.nyi_(); ++i )
		{
			for( size_t j=0; j < v.nxi_(); ++j )
			{
				os<<i<<' '<<j<<' ' << v[(1+i)*v.nxs_()+j+1] << "\n";        
			}  
		    //MPI_Barrier( cartcomm );//??
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
	Stencil(int nx,int ny, int * nbrs_):nx_(nx+1),ny_(ny+1),pi(3.141592653589793),nbrs(nbrs_){
	        double hx_ = 2.0/nx;
            double hy_ = 1.0/ny;
            xst =  1.0/(hx_*hx_);
            yst =  1.0/(hy_*hy_);
            mst = 2.0/(hx_*hx_)+2.0/(hy_*hy_)+4*pi*pi;          
            pg = (nx+1)*(ny+1);  
	}
	~Stencil(){}	
	/*
	Vector operator - (Vector& u){
		//std::cout<<"--";
		 Vector r(pg,u.nx_(),u.ny_());
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
		int x_()const{return b_.x_();}
		int y_()const{return a_.y_();}
		int nxs_()const{return b_.nxs_();}
	};

template <class A, class B>
	class Minus : public Expr<Minus<A,B>> {
	const A& a_;
	const B& b_;
	public :
		Minus(const A& a,const B& b): a_(a),b_(b){}
		
		int x_()const{return b_.x_();}
		int y_()const{return a_.y_();}
		int nxs_()const{return b_.nxs_();}
		
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
		int x_()const{return a_.x_();}
		int y_()const{return a_.y_();}
		int nxs_()const{return a_.nxs_();}
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
		return (b_.mst_()*a_[i] - b_.xst_()*(a_[i+1]+a_[i-1])-b_.yst_()*(a_[i+a_.nxs_()]+a_[i-a_.nxs_()]));
		}
		int x_()const{return a_.x_();}
		int y_()const{return a_.y_();}
		int nxs_()const{return a_.nxs_();}
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
	//MPI topology
	  int px = 2, py = 2;
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
	int x=0,y=0; //first inner point coordonates
	int offsetx=coords[0]*nx_, offsety =coords[1]*ny_;
	
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
			 [C,freqx,freqy](int x,int y)->double{return C*sin(freqx*x)*sinh(freqy*y);});	//set freq
			int last_row = nxs*(nys-1);	
			double SINH = sinh(2*pi); 
	//sin(2πx) sinh(2πy)		
	Vector u(pg,nx_,ny_,nxi,nyi,nxs,nys,x,y,offsetx, [SINH,freqx](int x)->double{return sin(x*freqx) * SINH;});
	
	Stencil A(nx,ny,nbrs);
	double delta0 = 0, delta1 = 0, beta = 0,alfa=0 ,scalar =0;
	double mydelta0 = 0, mydelta1 = 0,  myscalar = 0;
	//initialization
	
	
	for( int i = 0; i < 4; ++i )
      {
         if( cartrank == i )
         {

            std::cout<<f << "------------------------------------------------------------\n"
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
      }
	
	
	//CG
	r = f - A*u;
	//std::cout<<(A*u);//u<<(A*u);
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
	//std::cout<<u;
	//return sqrt(delta0);	
	
	std::ofstream out("solution.txt");
	//todo all sent to 0 proccess
	//out&u;//naive print
	
	//advanced print
	int a=1;
	MPI_Status status;	
	if(0)
	if(coords[0]==0 && coords[1]==0){//first block
		
	    out&u;
		a=nbrs[DOWN];
		MPI_Send( &a, 1, MPI_INT, nbrs[RIGHT], 10, cartcomm );
	}else if(coords[0]==py-1 && coords[1]==px-1){//last block
	   
		MPI_Recv( &a, 1, MPI_INT, nbrs[LEFT], 10, cartcomm, &status );
	   out&u;
	}else if(coords[1]==px-1){//last block in row
		
		MPI_Recv( &a, 1, MPI_INT, nbrs[LEFT], 10, cartcomm, &status );
	    out&u;
		MPI_Send( &a, 1, MPI_INT, a, 10, cartcomm );
	}else if(coords[1] == 0){//first block in row
		int last;
		int lastcoords[2] = {coords[0]-1,py-1};
		MPI_Cart_rank(cartcomm,lastcoords,&last);
	    MPI_Recv( &a, 1, MPI_INT, last, 10, cartcomm, &status );
	    out&u;
		a=nbrs[DOWN];
		MPI_Send( &a, 1, MPI_INT, nbrs[RIGHT], 10, cartcomm );
	}else{//inner block
		
	    MPI_Recv( &a, 1, MPI_INT, nbrs[LEFT], 10, cartcomm, &status );
	    out&u;
		MPI_Send( &a, 1, MPI_INT, nbrs[RIGHT], 10, cartcomm );
	}
	
	return sqrt(delta0);	
}
