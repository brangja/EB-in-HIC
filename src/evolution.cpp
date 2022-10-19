void Maxwell::Compute_EB(int time_step )
{

	int it,i,j,k;
	it=time_step;	
	
	//int ix,iy,iz;
	
	
	//exchange the data of Ex0 and Ex1
	#pragma omp parallel for collapse(3) 
	for(int ix=0; ix<Nx; ix++){
	  for(int iy=0; iy<Ny; iy++){		
	    for(int iz=0; iz<Nz; iz++){     
	      
	      
	      int id0=ix*Ny*Nz+iy*Nz+iz;
	      
	      //E field
	      Ex0[id0]=Ex1[id0]; 
	      Ey0[id0]=Ey1[id0];
	      Ez0[id0]=Ez1[id0];
	      
	      
	      //B field
	      Bx0[id0]=Bx1[id0];
	      By0[id0]=By1[id0];
	      Bz0[id0]=Bz1[id0];
	     
	    /*Bx1[id]=0.0; 
	      By1[id]=0.0;
	      Bz1[id]=0.0;
	      Ex1[id]=0.0; 
	      Ey1[id]=0.0;
	      Ez1[id]=0.0;
	      */
	      
	     
	     
	    }
	  }
	}
	
	
	#pragma omp parallel for collapse(3) 
	for(int ix=0; ix<Nx; ix++){
	  for(int iy=0; iy<Ny; iy++){		
	    for(int iz=0; iz<Nz; iz++){     
	      
	      int id=0;
		  int ixf, iyf, izf;
		  int idxf, idyf, idzf;
		  int ixb, iyb, izb;
		  int idxb, idyb, idzb;

	      ixf=ix+1; iyf=iy+1; izf=iz+1;
		  if(ix==Nx-1)ixf=ix;
		  if(iy==Ny-1)iyf=iy;
		  if(iz==Nz-1)izf=iz;
		  id=ix*Ny*Nz+iy*Nz+iz;
		  idxf=ixf*Ny*Nz+iy*Nz+iz;
		  idyf=ix*Ny*Nz+iyf*Nz+iz;
		  idzf=ix*Ny*Nz+iy*Nz+izf;

		  ixb=ix-1; iyb=iy-1; izb=iz-1;
		  if(ix==0)ixb=ix;
		  if(iy==0)iyb=iy;
		  if(iz==0)izb=iz;
		  id=ix*Ny*Nz+iy*Nz+iz;
		  idxb=ixb*Ny*Nz+iy*Nz+iz;
		  idyb=ix*Ny*Nz+iyb*Nz+iz;
		  idzb=ix*Ny*Nz+iy*Nz+izb;


			Bx1[id]=Bx0[id]+ dt/dz*(Ey0[idzf]-Ey0[id]) \
	      				     -dt/dy*(Ez0[idyf]-Ez0[id]);
	      
	    	By1[id]=By0[id]+ dt/dx*(Ez0[idxf]-Ez0[id]) \
	      				     -dt/dz*(Ex0[idzf]-Ex0[id]);
	      
	      
	      	Bz1[id]=Bz0[id]+ dt/dy*(Ex0[idyf]-Ex0[id]) \
	      				     -dt/dx*(Ey0[idxf]-Ey0[id]);
		 
	      
	     
	     
	    }
	  }
	}

	#pragma omp parallel for collapse(3) 
	for(int ix=0; ix<Nx; ix++){
	  for(int iy=0; iy<Ny; iy++){		
	    for(int iz=0; iz<Nz; iz++){     
	      
	      int id=0;
		  int ixf, iyf, izf;
		  int idxf, idyf, idzf;
		  int ixb, iyb, izb;
		  int idxb, idyb, idzb;

	      ixf=ix+1; iyf=iy+1; izf=iz+1;
		  if(ix==Nx-1)ixf=ix;
		  if(iy==Ny-1)iyf=iy;
		  if(iz==Nz-1)izf=iz;
		  id=ix*Ny*Nz+iy*Nz+iz;
		  idxf=ixf*Ny*Nz+iy*Nz+iz;
		  idyf=ix*Ny*Nz+iyf*Nz+iz;
		  idzf=ix*Ny*Nz+iy*Nz+izf;

		  ixb=ix-1; iyb=iy-1; izb=iz-1;
		  if(ix==0)ixb=ix;
		  if(iy==0)iyb=iy;
		  if(iz==0)izb=iz;
		  id=ix*Ny*Nz+iy*Nz+iz;
		  idxb=ixb*Ny*Nz+iy*Nz+iz;
		  idyb=ix*Ny*Nz+iyb*Nz+iz;
		  idzb=ix*Ny*Nz+iy*Nz+izb;

		      
	      Ex1[id]=COE(it+0.5, ix+0.5, iy, iz)*Ex0[id]+CE(it+0.5, ix+0.5, iy, iz)*( dt/dy*(Bz1[id]-Bz1[idyb])\
	      				     -dt/dz*(By1[id]-By1[idzb])  -dt*jsx(it+0.5, ix+0.5, iy, iz) );
	      				     //-dt/8.0*sigma_chi(it+0.5, ix+0.5, iy, iz)*( Bx1[idxf]+ Bx1[(ix+1)*Ny*Nz+iy*Nz+iz-1]+ Bx1[(ix+1)*Ny*Nz+(iy-1)*Nz+iz]\
	      				     +Bx1[(ix+1)*Ny*Nz+(iy-1)*Nz+iz-1]+Bx1[ix*Ny*Nz+iy*Nz+iz] +Bx1[ix*Ny*Nz+iy*Nz+iz-1] +Bx1[ix*Ny*Nz+(iy-1)*Nz+iz]\
	      				     +Bx1[ix*Ny*Nz+(iy-1)*Nz+iz-1] )   );
	             
	      
	      Ey1[id]=COE(it+0.5, ix, iy+0.5, iz)*Ey0[id]+CE(it+0.5, ix, iy+0.5, iz)*( dt/dz*(Bx1[id]-Bx1[idzb]) \
	      				     -dt/dx*(Bz1[id]-Bz1[idxb])   -dt*jsy(it+0.5, ix, iy+0.5, iz)  );
	      				     //-dt/8.0*sigma_chi(it+0.5, ix, iy+0.5, iz)*( By1[ix*Ny*Nz+(iy+1)*Nz+iz]+ By1[ix*Ny*Nz+(iy+1)*Nz+iz-1]+ By1[(ix-1)*Ny*Nz+(iy+1)*Nz+iz]\
	      				     +By1[(ix-1)*Ny*Nz+(iy+1)*Nz+iz-1]+By1[ix*Ny*Nz+iy*Nz+iz] +By1[ix*Ny*Nz+iy*Nz+iz-1] +By1[(ix-1)*Ny*Nz+iy*Nz+iz]\
	      				     +By1[(ix-1)*Ny*Nz+iy*Nz+iz-1] )   );
	     
	      Ez1[id]=COE(it+0.5, ix, iy, iz+0.5)*Ez0[id]+CE(it+0.5, ix, iy, iz+0.5)*( dt/dx*(By1[id]-By1[idxb])\
	       				     -dt/dy*(Bx1[id]-Bx1[idyb])  -dt*jsz(it+0.5, ix, iy, iz+0.5) );
	       				     //-dt/8.0*sigma_chi(it+0.5, ix, iy, iz+0.5)*( Bz1[ix*Ny*Nz+iy*Nz+iz+1]+ Bz1[ix*Ny*Nz+(iy-1)*Nz+iz+1]+ Bz1[(ix-1)*Ny*Nz+iy*Nz+iz+1]\
	      				     +Bz1[(ix-1)*Ny*Nz+(iy-1)*Nz+iz+1]+Bz1[ix*Ny*Nz+iy*Nz+iz] +Bz1[ix*Ny*Nz+(iy-1)*Nz+iz] +Bz1[(ix-1)*Ny*Nz+iy*Nz+iz]\
	      				     +Bz1[(ix-1)*Ny*Nz+(iy-1)*Nz+iz] )  );

		   	      
	     
	     
	    }
	  }
	}



}

