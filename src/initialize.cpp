
Maxwell::Maxwell(int time_number)
{
	 
	//open the out files
	f_Ex.open("./data/Ex.dat");
	f_Ey.open("./data/Ey.dat");
	f_Ez.open("./data/Ez.dat");
	f_Bx.open("./data/Bx.dat");
	f_By.open("./data/By.dat");
	f_Bz.open("./data/Bz.dat");
	
	
	//calloc E,B fields
	Ex0 = (double *) calloc (Nx*Ny*Nz, sizeof(double));
	Ey0 = (double *) calloc (Nx*Ny*Nz, sizeof(double));
	Ez0 = (double *) calloc (Nx*Ny*Nz, sizeof(double));
	Bx0 = (double *) calloc (Nx*Ny*Nz, sizeof(double));
	By0 = (double *) calloc (Nx*Ny*Nz, sizeof(double));
	Bz0 = (double *) calloc (Nx*Ny*Nz, sizeof(double));
	Ex1 = (double *) calloc (Nx*Ny*Nz, sizeof(double));
	Ey1 = (double *) calloc (Nx*Ny*Nz, sizeof(double));
	Ez1 = (double *) calloc (Nx*Ny*Nz, sizeof(double));
	Bx1 = (double *) calloc (Nx*Ny*Nz, sizeof(double));
	By1 = (double *) calloc (Nx*Ny*Nz, sizeof(double));
	Bz1 = (double *) calloc (Nx*Ny*Nz, sizeof(double));
	
	
       	  
        
	
	
}

Maxwell::~Maxwell(void)
{
	//close the out files
	f_Ex.close();
	f_Ey.close();
	f_Ez.close();
	f_Bx.close();
	f_Bx.close();
	f_Bx.close();
	
	
	//free the pointer
	free(Ex0); 
	free(Ey0);
	free(Ez0);
	free(Bx0); 
	free(By0);
	free(Bz0);
	free(Ex1); 
	free(Ey1);
	free(Ez1);
	free(Bx1); 
	free(By1);
	free(Bz1);
	


}



void Maxwell::Initialize_EB()
{
	int t, i,j,k;
	
	
}

