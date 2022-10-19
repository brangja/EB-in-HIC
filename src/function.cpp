//density of spactator
double Maxwell ::density(double x, double y, double z)
{
	double R=Radius;
	double gamma=1.0/sqrt(1.0-Beta*Beta) ;
	double r=sqrt( pow(x, 2)+pow(y,2)+pow(z,2) );
	
	double a=0.535;
	double w=0.0;
	double rho_0=gamma*Zp/1163.2987;
	
	double rho=rho_0*( 1+w*pow(r/R, 2) )/( 1+exp( (r-R)/a ) );
	
	
	
	return rho;
	

}

double Maxwell ::charge(double r)
{
	double R=Radius;
	double qr=-0.145794 - 0.0413548*exp(r) + 0.719714*r - 0.479308*pow(r,2) + 0.396885*pow(r,3);
	
	if( r>R ) return Zp;
	else return qr;
	

}



double Maxwell :: EBsingle_ext(double Rx, double Ry, double Rz, double xi, int path) //xi is denoted the direction of velocity along the z axis
{
	int sn=path;
	double v=xi*Beta;
	double gamma=sqrt(1.0/(1.0-pow(v,2)));
	double vcrossRx=-v*Ry;
	double vcrossRy=v*Rx;
	double vcrossRz=0.0;
	double Denominator=pow( pow(Rx,2)+pow(Ry,2)+gamma*gamma*pow(Rz,2),1.5); 
	
	
	double eb=0.0;
	switch(sn)
	{
	case 0: {eb=alpha*gamma*Rx/Denominator; break;}
	case 1: {eb=alpha*gamma*Ry/Denominator; break;}
	case 2: {eb=alpha*gamma*Rz/Denominator; break;}
	case 3: {eb=alpha*gamma*vcrossRx/Denominator; break;}
	case 4: {eb=alpha*gamma*vcrossRy/Denominator; break;}
	case 5: {eb=alpha*gamma*vcrossRz/Denominator; break;}
	default: {cout<<"something orrer !"<<endl; break;}
	
	}
	
	if(Denominator<1e-10) return 0.0; //just for avoid the singularity
	else return eb;
}


//EB for external electric field
double Maxwell :: EB_ext(double t, double x, double y, double z, int order)
{
		
	double b=Impact_parameter;
	int sn=order-1;
	
	double R=Radius;
	double eb=0.0;
	double v=Beta;
	double gamma=sqrt(1.0/(1.0-pow(v,2)));
	
	
	double xp1, xp2, yp1, yp2, Rx1, Ry1, Rx2, Ry2, Rz1, Rz2, Rsqrt1, Rsqrt2;
	xp1=-b/2.0; yp1=0.0; xp2=b/2.0; yp2=0.0;
    Rx1=x-xp1; Ry1=y-yp1; Rz1=z-v*t;  Rx2=x-xp2; Ry2=y-yp2; Rz2=z+v*t;
    Rsqrt1=sqrt( pow(Rx1,2)+pow(Ry1,2)+gamma*gamma*pow(Rz1,2) );  Rsqrt2=sqrt( pow(Rx2,2)+pow(Ry2,2)+gamma*gamma*pow(Rz2,2) );
    eb= charge(Rsqrt1)*EBsingle_ext(Rx1, Ry1, Rz1, 1.0, sn)+ charge(Rsqrt2)*EBsingle_ext(Rx2, Ry2, Rz2, -1.0, sn) ;
        
	return eb;

}



double Maxwell :: sigma_v(double it, double ix, double iy, double iz)
{
		
	double t=t0+it*dt; 
	double x=x0+ix*dx;
	double y=y0+iy*dy;
	double z=dz+iz*dz;
	
	double t_0=time_hydro;
	double Tc=135.0/197.0;
	double Tau_0=3*Tc;
	double sigma_t_Bjorken=sigmaTc*Tau_0/Tc*(t_0/(t+1e-10));
	double sigma_t_Hubble=sigmaTc*Tau_0/Tc*pow(t_0/(t+1e-10), 1.0/3.0);
	double sigma_tau0=sigmaTc*Tau_0/Tc;
	double sigma_tx=0.0;
	
	double R=Radius, b=Impact_parameter;
	double Rx=R-b/2.0, Ry=sqrt(R*R-b*b/4.0), Rz=0.1; 
	double r=sqrt(x*x+y*y+z*z);
	double r_range=3.0;//the range of radius of the QGP 
	double spacefactor=exp(-r*r/20.0);

	switch(Model)
	{
		case 000: if(t<time_hydro)sigma_tx=0.0; else sigma_tx=sigmaTc; break;
		case 001: if(t<time_HIC) sigma_tx=0.0; else sigma_tx=sigmaTc; break;
		case 111: if(t>=time_hydro) sigma_tx=sigma_t_Bjorken; else sigma_tx=0.0; break; //Bjorken-Homogeneous-vacuum
		// case 111: {
		// 		if(r>3.0) sigma_tx=0.0; 
		// 		else {if(t>t_0) sigma_tx=sigma_t_Bjorken; else sigma_tx=0.0;} 
		// 		break;	
		// 	  }//Bjorken-Sphere-vacuum
		case 121: if(t>=time_hydro) sigma_tx=sigma_t_Bjorken*spacefactor; else sigma_tx=0.0; break; //Bjorken-Gauss-vacuum
		case 211: if(t>=time_hydro) sigma_tx=sigma_t_Hubble; else sigma_tx=0.0; break; //Hubble-Homogeneous-vacuum
		// case 211: {
		// 		if(r>3.0) sigma_tx=0.0; 
		// 		else {if(t>t_0) sigma_tx=sigma_t_Hubble; else sigma_tx=0.0;} 
		// 		break;	
		// 	  }//Hubble-Sphere-vacuum
		case 221: if(t>=time_hydro) sigma_tx=sigma_t_Hubble*spacefactor; else sigma_tx=0.0; break; //Hubble-Gauss-vacuum
		case 112: if(t>=time_hydro) sigma_tx=sigma_t_Bjorken; else if(t<time_HIC) sigma_tx=0.0; else sigma_tx=sigma_tau0; break; //Bjorken-Homogeneous-turn-on
		// case 112: {
		// 		if(r>3.0) sigma_tx=0.0; 
		// 		else {if(t>t_0) sigma_tx=sigma_t_Bjorken; else if(t<0.0) sigma_tx=0.0; else sigma_tx=sigma_tau0; } 
		// 		break;	
		// 	  }//Bjorken-Sphere-turn-on
		case 122: if(t>=time_hydro) sigma_tx=sigma_t_Bjorken*spacefactor; else if(t<time_HIC) sigma_tx=0.0; else sigma_tx=sigma_tau0*spacefactor; break;//Bjorken-Gauss-turn-on
		case 212: if(t>=time_hydro) sigma_tx=sigma_t_Hubble; else if(t<time_HIC) sigma_tx=0.0; else sigma_tx=sigma_tau0; break; //Hubble-Homogeneous-turn-on
		// case 212: {
		// 		if(r>3.0) sigma_tx=0.0; 
		// 		else {if(t>t_0) sigma_tx=sigma_t_Hubble; else if(t<0.0) sigma_tx=0.0; else sigma_tx=sigma_tau0; } 
		// 		break;	
		// 	  }//Hubble-Sphere-turn-on
		case 222: if(t>time_hydro) sigma_tx=sigma_t_Hubble*spacefactor; else if(t<time_HIC) sigma_tx=0.0; else sigma_tx=sigma_tau0*spacefactor; break;//Hubble-Gauss-turn-on
		default: cout<<"something wrong !" << endl;
	}
	
	
	
	return sigma_tx;
	
	//return 10.0*sigmaTc;

}


double Maxwell :: sigma_chi(double it, double ix, double iy, double iz)
{
		
	double t=t0+it*dt; 
	double x=x0+ix*dx;
	double y=y0+iy*dy;
	double z=z0+iz*dz;
	double mu5=5.0; //unit of fm^-1
	double Nc=3.0; //just consider u,d,s quark
	double sigmachi=2.0/PI*alpha*Nc*6.0/9.0*mu5; 
	
	
	if(t>=time_hydro) return sigmachi;
	else return 0.0;
	  

}


double Maxwell :: jsx(double it, double ix, double iy, double iz )
{
		
	
	double t=t0+it*dt; 
	double x=x0+ix*dx;
	double y=y0+iy*dy;
	double z=z0+iz*dz;
	double jx=sigma_v(it, ix, iy, iz)*EB_ext(t, x, y, z, 1)+sigma_chi(it, ix, iy, iz)*EB_ext(t, x, y, z, 4);
	
	return   jx;  

}

double Maxwell :: jsy(double it, double ix, double iy, double iz )
{
		
	double t=t0+it*dt; 
	double x=x0+ix*dx;
	double y=y0+iy*dy;
	double z=z0+iz*dz;	
	double jy=sigma_v(it, ix, iy, iz)*EB_ext(t, x, y, z, 2)+sigma_chi(it, ix, iy, iz)*EB_ext(t, x, y, z, 5);
	
	return  jy;  

}

double Maxwell :: jsz(double it, double ix, double iy, double iz)
{
		
	double t=t0+it*dt; 
	double x=x0+ix*dx;
	double y=y0+iy*dy;
	double z=z0+iz*dz;	
	double jz=sigma_v(it, ix, iy, iz)*EB_ext(t, x, y, z, 3)+sigma_chi(it, ix, iy, iz)*EB_ext(t, x, y, z, 6);
	
	return jz ;  

}

double Maxwell :: COE(double it, double ix, double iy, double iz)
{
	
	
	double ce_p=2.0/(2.0+dt*sigma_v(it, ix, iy, iz) );
	double ce_m=2.0/(2.0-dt*sigma_v(it, ix, iy, iz) );
	double coe=ce_p/ce_m;
	
	
	return  coe;

}

double Maxwell :: CE(double it, double ix, double iy, double iz)
{
	
	double ce_p=2.0/(2.0+dt*sigma_v(it, ix, iy, iz) );
		
	
	return  ce_p;

}
