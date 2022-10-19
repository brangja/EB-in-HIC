struct vec{ 
    double x; 
    double y;
    double z;
    double f;};
	
class Interpolation
{
private:
	
public:
	double S1(double x)
	{
		
		double res=0.0;
		if(fabs(x)<=1.0) res=1.0-fabs(x);
		else res=0.0;
		return res;

	};
	
	double d1( struct vec *array, int nx, double xp );
	double d2( struct vec *array, int nx, int ny, double xp, double yp );
	double d3( struct vec *array, int nx, int ny, int nz, double xp, double yp, double zp);
	
};

double Interpolation::d1( struct vec *array, int nx, double xp )
{
	//check the point in the range of data
	double x0, xf, dx;
	x0=array[0].x; xf=array[nx-1].x; dx=array[1].x-array[0].x;
	 
	if(xp<x0 || xp> xf) { cout<<"overflow the range of data" <<endl; return 0.0; }

	int i, ix, ixe, ixen;
	double res=0.0, xi, dxi;

	ix=floor( (xp-x0)/dx );
	//the range of integration 
	ixe=ix-1; ixen=ix+1;
	if(ixe<0) ixe=ix;
	if(ixen>nx-1) ixen=ix;

	//the integration 
	for(i=ixe; i<=ixen; i++)
	{
		xi=x0+i*dx; 
		dxi=(xp-xi)/dx;
		
		res=res+array[i].f*S1(dxi);
		
	}

	return res;
		
}

double Interpolation::d2( struct vec *array, int nx, int ny, double xp, double yp )
{
	//check the point in the range of data
	double x0, xf, dx, y0, yf, dy;
	x0=array[0].x; xf=array[(nx-1)*ny+ny-1].x; dx=array[1*ny].x-array[0].x;
	y0=array[0].y; yf=array[0*ny+ny-1].y; dy=array[0*ny+1].y-array[0].y;

	//cout<< "y0="<<y0<<"\tyf="<<yf<<"\tdy="<<dy<<endl; 
	 
	if(xp<x0 || xp> xf || yp<y0 || yp> yf) { cout<<"overflow the range of data" <<endl; return 0.0; }

	int i, ix, ixe, ixen;
	int j, iy, iye, iyen;
	double res=0.0, xi, dxi, yi, dyi;

	ix=floor( (xp-x0)/dx );
	iy=floor( (yp-y0)/dy );
	//the range of integration 
	ixe=ix-1; ixen=ix+1;
	if(ixe<0) ixe=ix;
	if(ixen>nx-1) ixen=ix;

	iye=iy-1; iyen=iy+1;
	if(iye<0) iye=iy;
	if(iyen>ny-1) iyen=iy;

	//the integration 
	for(i=ixe; i<=ixen; i++)
		for(j=iye; j<=iyen; j++)
			{
				xi=x0+i*dx; yi=y0+j*dy; 
				dxi=(xp-xi)/dx; dyi=(yp-yi)/dy;
				
				res=res+array[i*ny+j].f*S1(dxi)*S1(dyi);
				
			}

	return res;
		
}



double Interpolation::d3(  struct vec *array, int nx, int ny, int nz, double xp, double yp, double zp)
{
	//check the point in the range of data
	double x0, xf, dx, y0, yf, dy, z0, zf, dz;
	x0=array[0].x; xf=array[(nx-1)*ny*nz+(ny-1)*nz+nz-1].x; dx=array[ny*nz].x-array[0].x;
	y0=array[0].y; yf=array[(ny-1)*nz+nz-1].y; dy=array[nz].y-array[0].y;
	z0=array[0].z; zf=array[nz-1].z; dz=array[1].z-array[0].z;

	// cout<< "x0="<<x0<<"\txf="<<xf<<"\tdx="<<dx<<endl;
	// cout<< "y0="<<y0<<"\tyf="<<yf<<"\tdy="<<dy<<endl;
	// cout<< "z0="<<z0<<"\tzf="<<zf<<"\tdz="<<dz<<endl; 
	 
	if(xp<x0 || xp> xf || yp<y0 || yp> yf || zp<z0 || zp> zf) { cout<<"overflow the range of data" <<endl; return 0.0; }

	int i, ix, ixe, ixen;
	int j, iy, iye, iyen;
	int k, iz, ize, izen;
	double res=0.0, xi, dxi, yi, dyi, zi, dzi;

	ix=floor( (xp-x0)/dx );
	iy=floor( (yp-y0)/dy );
	iz=floor( (zp-z0)/dz );
	//the range of integration 
	ixe=ix-1; ixen=ix+1;
	if(ixe<0) ixe=ix;
	if(ixen>nx-1) ixen=ix;

	iye=iy-1; iyen=iy+1;
	if(iye<0) iye=iy;
	if(iyen>ny-1) iyen=iy;

	ize=iz-1; izen=iz+1;
	if(ize<0) ize=iz;
	if(izen>nz-1) izen=iz;

	//the integration 
	for(i=ixe; i<=ixen; i++)
		for(j=iye; j<=iyen; j++)
			for(k=ize; k<=izen; k++)
			{
				xi=x0+i*dx; yi=y0+j*dy; zi=z0+k*dz;
				dxi=(xp-xi)/dx; dyi=(yp-yi)/dy; dzi=(zp-zi)/dz;
				
				res=res+array[i*ny*nz+j*nz+k].f*S1(dxi)*S1(dyi)*S1(dzi);
				
			}

	return res;
		
}
