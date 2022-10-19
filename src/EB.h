
class Maxwell
{
private:
	
	
	
	double *Ex0;
	double *Ey0;
	double *Ez0;
	double *Bx0;
	double *By0;
	double *Bz0;
	double *Ex1;
	double *Ey1;
	double *Ez1;
	double *Bx1;
	double *By1;
	double *Bz1;
	
	ofstream f_Ex;
	ofstream f_Ey;
	ofstream f_Ez;
	ofstream f_Bx;
	ofstream f_By;
	ofstream f_Bz;
	//coefficient
	double sigma_v(double t, double i, double j, double k);
	double sigma_chi(double t, double i, double j, double k);
	double COE(double t, double i, double j, double k);
	double CE(double t, double i, double j, double k);
	//source
	double jsx(double t, double i, double j, double k);
	double jsy(double t, double i, double j, double k);
	double jsz(double t, double i, double j, double k);
	
	//external field
	/*double Ex_ext(double t, double i, double j, double k);
	double Ey_ext(double t, double i, double j, double k);
	double Ez_ext(double t, double i, double j, double k);
	double Bx_ext(double t, double i, double j, double k);
	double By_ext(double t, double i, double j, double k);
	double Bz_ext(double t, double i, double j, double k);
	*/
	double EB_ext(double t, double i, double j, double k, int sn);
	double unitestep(double x);
	double charge(double r);
	double density(double x, double y, double z);
	double EBsingle_ext(double x, double y, double z, double xi, int path);
	
	//differential
	double min3(double x, double y, double z);
	double max3(double x, double y, double z);
	double min_mod(double step, double f_minus, double f_mid, double f_plus);
	double pde(double *Var, int it, int ix, int iy, int iz, int direction);
	int sign(double x);	

	//============================
	double f(double x, double y, double z);
	double F(string ebfield, double var, int direction );
	double Co(string ebfield, int it, int ix, int iy, int iz, int direction, int si);
	double Hf(string ebfield, double *var, int it, int ix, int iy, int iz, int direction, int si);
	double FF(string ebfield, double *var, double *Ex0, double *Ey0, double *Ez0, double *Bx0, double *By0, double *Bz0, int it, int ix, int iy, int iz);
	
	
	
	
	
	
	
	
public:
	Maxwell(int time_number);
	~Maxwell();
	void Initialize_EB();
	void Compute_EB(int time_step);
	//void Getall();
	void Output_EB(int t);
	void  Evolution(string ebfeiel, double *var, double *Ex0, double *Ey0, double *Ez0, double *Bx0, double *By0, double *Bz0, int it);
	friend class Interpolation;
};








