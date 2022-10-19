void Maxwell::Output_EB(int it)
{
	double t, x, y, z;
	
	//reconstruct the data for E and B fields
	// double dataEx[Nx*Ny*Nz][4]={0}, dataEy[Nx*Ny*Nz][4]={0}, dataEz[Nx*Ny*Nz][4]={0};
	// double dataBx[Nx*Ny*Nz][4]={0}, dataBy[Nx*Ny*Nz][4]={0}, dataBz[Nx*Ny*Nz][4]={0};
	struct vec *dataEx= (struct vec *) calloc (Nx*Ny*Nz, sizeof(struct vec));
	struct vec *dataEy= (struct vec *) calloc (Nx*Ny*Nz, sizeof(struct vec));
	struct vec *dataEz= (struct vec *) calloc (Nx*Ny*Nz, sizeof(struct vec));
	struct vec *dataBx= (struct vec *) calloc (Nx*Ny*Nz, sizeof(struct vec));
	struct vec *dataBy= (struct vec *) calloc (Nx*Ny*Nz, sizeof(struct vec));
	struct vec *dataBz= (struct vec *) calloc (Nx*Ny*Nz, sizeof(struct vec));

	#pragma omp parallel for collapse(3)
    for(int i=0; i<Nx; i++  )
        for(int j=0; j<Ny; j++  )
            for(int k=0; k<Nz; k++  )
            {
                int id=i*Ny*Nz+j*Nz+k;
				//Ex
				x=x0+(i+0.5)*dx ; y=y0+j*dy; z=z0+k*dz;
                dataEx[id].x=x; dataEx[id].y=y; dataEx[id].z=z;
                dataEx[id].f=Ex0[id];

				//Ey
				x=x0+i*dx ; y= y0+(j+0.5)*dy; z=z0+k*dz;
                dataEy[id].x=x; dataEy[id].y=y; dataEy[id].z=z;
                dataEy[id].f=Ey0[id];

				//Ez
				x=x0+i*dx ; y= y0+j*dy; z=z0+(k+0.5)*dz;
                dataEz[id].x=x; dataEz[id].y=y; dataEz[id].z=z;
                dataEz[id].f=Ez0[id];


				//Bx
				x=x0+i*dx ; y=y0+(j+0.5)*dy; z=z0+(k+0.5)*dz;
                dataBx[id].x=x; dataBx[id].y=y; dataBx[id].z=z;
                dataBx[id].f=Bx1[id];
				//dataBx[id].f=(Bx1[id]+Bx0[id])/2.0; //this is for get the interger time

				//By
				x=x0+(i+0.5)*dx ; y= y0+j*dy; z=z0+(k+0.5)*dz;
                dataBy[id].x=x; dataBy[id].y=y; dataBy[id].z=z;
                dataBy[id].f=By1[id];
				//dataBy[id].f=(By1[id]+By0[id])/2.0;

				//Bz
				x=x0+(i+0.5)*dx ; y= y0+(j+0.5)*dy; z=z0+k*dz;
                dataBz[id].x=x; dataBz[id].y=y; dataBz[id].z=z;
                dataBz[id].f=Bz1[id];
				//dataBz[id].f=(Bz1[id]+Bz0[id])/2.0;
				
            }
	
	//output the data at center point
	if(Outdensity==0)
	{
		t=t0+it*dt;
		x=0.0; y=0.0; z=0.0;

		double ex=0.0, ey=0.0, ez=0.0;
		double bx=0.0, by=0.0, bz=0.0;

		//interpolation
		Interpolation interpolation;
		//iner field
		ex=interpolation.d3( dataEx, Nx, Ny, Nz, x, y, z);
		ey=interpolation.d3( dataEy, Nx, Ny, Nz, x, y, z);
		ez=interpolation.d3( dataEz, Nx, Ny, Nz, x, y, z);
		bx=interpolation.d3( dataBx, Nx, Ny, Nz, x, y, z);
		by=interpolation.d3( dataBy, Nx, Ny, Nz, x, y, z);
		bz=interpolation.d3( dataBz, Nx, Ny, Nz, x, y, z);
		//external field
		double Ex_ext=EB_ext(t, x, y, z, 1);
		double Ey_ext=EB_ext(t, x, y, z, 2);
		double Ez_ext=EB_ext(t, x, y, z, 3);
		double Bx_ext=EB_ext(t+0.5*dt, x, y, z, 4);
		double By_ext=EB_ext(t+0.5*dt, x, y, z, 5);
		double Bz_ext=EB_ext(t+0.5*dt, x, y, z, 6);
		//output the data
		f_Ex<< t <<"\t"<< x<<"\t" << y <<"\t" <<z <<"\t" << ex<<"\t" <<Ex_ext <<"\t" << ex+Ex_ext <<endl;
		f_Ey<< t <<"\t"<< x<<"\t" << y <<"\t" <<z <<"\t" << ey<<"\t" <<Ey_ext <<"\t" << ey+Ey_ext <<endl;
		f_Ez<< t <<"\t"<< x<<"\t" << y <<"\t" <<z <<"\t" << ez<<"\t" <<Ez_ext <<"\t" << ez+Ez_ext <<endl;

		f_Bx<< t+0.5*dt <<"\t"<< x<<"\t" << y <<"\t" <<z <<"\t" << bx<<"\t" <<Bx_ext <<"\t" << bx+Bx_ext <<endl;
		f_By<< t+0.5*dt <<"\t"<< x<<"\t" << y <<"\t" <<z <<"\t" << by<<"\t" <<By_ext<<"\t" << by+By_ext <<endl;
		f_Bz<< t+0.5*dt <<"\t"<< x<<"\t" << y <<"\t" <<z <<"\t" << bz<<"\t" <<Bz_ext <<"\t" << bz+Bz_ext <<endl;
	}
	

	//output the data at any space time
	if(Outdensity==1)
	{
		t=t0+it*dt;
		int nx=11, ny=11, nz=1;
		double xn0=-5.0, dnx=1.0;
		double yn0=-5.0, dny=1.0;
		double zn0=-0.0, dnz=0.2;
		x=0.0; y=0.0; z=0.0;

		double ex=0.0, ey=0.0, ez=0.0;
		double bx=0.0, by=0.0, bz=0.0;

		//interpolation
		Interpolation interpolation;
		for(int i=0; i<nx; i++  )
			for(int j=0; j<ny; j++  )
				for(int k=0; k<nz; k++  )
				{
					x=xn0+i*dnx; y=yn0+j*dny; z=zn0+k*dnz;
					//iner field
					ex=interpolation.d3( dataEx, Nx, Ny, Nz, x, y, z);
					ey=interpolation.d3( dataEy, Nx, Ny, Nz, x, y, z);
					ez=interpolation.d3( dataEz, Nx, Ny, Nz, x, y, z);
					bx=interpolation.d3( dataBx, Nx, Ny, Nz, x, y, z);
					by=interpolation.d3( dataBy, Nx, Ny, Nz, x, y, z);
					bz=interpolation.d3( dataBz, Nx, Ny, Nz, x, y, z);
					//external field
					double Ex_ext=EB_ext(t, x, y, z, 1);
					double Ey_ext=EB_ext(t, x, y, z, 2);
					double Ez_ext=EB_ext(t, x, y, z, 3);
					double Bx_ext=EB_ext(t+0.5*dt, x, y, z, 4);
					double By_ext=EB_ext(t+0.5*dt, x, y, z, 5);
					double Bz_ext=EB_ext(t+0.5*dt, x, y, z, 6);
					//output the data
					f_Ex<< t <<"\t"<< x<<"\t" << y <<"\t" <<z <<"\t" << ex<<"\t" <<Ex_ext <<"\t" << ex+Ex_ext <<endl;
					f_Ey<< t <<"\t"<< x<<"\t" << y <<"\t" <<z <<"\t" << ey<<"\t" <<Ey_ext <<"\t" << ey+Ey_ext <<endl;
					f_Ez<< t <<"\t"<< x<<"\t" << y <<"\t" <<z <<"\t" << ez<<"\t" <<Ez_ext <<"\t" << ez+Ez_ext <<endl;

					f_Bx<< t+0.5*dt <<"\t"<< x<<"\t" << y <<"\t" <<z <<"\t" << bx<<"\t" <<Bx_ext <<"\t" << bx+Bx_ext <<endl;
					f_By<< t+0.5*dt<<"\t"<< x<<"\t" << y <<"\t" <<z <<"\t" << by<<"\t" <<By_ext<<"\t" << by+By_ext <<endl;
					f_Bz<< t+0.5*dt <<"\t"<< x<<"\t" << y <<"\t" <<z <<"\t" << bz<<"\t" <<Bz_ext <<"\t" << bz+Bz_ext <<endl;
					
				}
	}
	
	
	free(dataEx);
	free(dataEy);
	free(dataEz);
	free(dataBx);
	free(dataBy);
	free(dataBz);

	

}
