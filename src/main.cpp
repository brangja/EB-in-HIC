#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <cstring>
#include<omp.h>

using namespace std;
#include "interpolation.h"
#include "EB.h"
#include "parameter.cpp"
#include "function.cpp"
#include "initialize.cpp"
#include "evolution.cpp"
#include "output.cpp"



int main()
{
	
	
	Maxwell EB(Nt);
	EB.Initialize_EB();
	cout<<"initialzation is finished" <<endl;
	double t,time_open;

	if(Medium_open==1) time_open=time_hydro;
	else if(Medium_open==2) time_open=time_HIC;

	for (int it=0; it<Nt; it++)
	{
		t=t0+it*dt;
		if(t<time_open)
			EB.Output_EB(it);
		else 
		{
			EB.Compute_EB(it);
			EB.Output_EB(it);
		}
		cout<<"this is time-step: "<<it <<endl;
	}
	
	cout<<"program is finished" <<endl;
	
	
	return 0;
}

