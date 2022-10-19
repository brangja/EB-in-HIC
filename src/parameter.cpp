// #define x0 -10.0
// #define y0 -10.0
// #define z0 -8.0
// #define t0 -1.5
// #define Nx 101
// #define Ny 101
// #define Nz 161
// #define Nt 300
// #define dx 0.2
// #define dy 0.2
// #define dz 0.1
// #define dt (dz/10.0)
// #define time_hydro 0.0
// #define Beta 0.999956 //for the sqrt(S)=200Gev. beta=sqrt(1-(2mp/sqrt(S))^2)
// #define PI 3.14159265
// #define Impact_parameter 6.0
// #define alpha 0.0072993
// #define Radius 6.38
// #define Zp 79.0 // number of proton of heavy-ion


#define x0 -11.05
#define y0 -11.0
#define z0 -6.05
#define t0 -1.5
#define Nx 221
#define Ny 221
#define Nz 121
#define Nt 200
#define dx 0.1
#define dy 0.1
#define dz 0.1
#define dt (dx/2.0)
#define time_hydro 0.4 //unit of fm, the hydro start time
#define time_HIC 0.0 //unit of fm, the time of two heavy ions combinning each other
#define Beta 0.99995 //for the sqrt(S)=200Gev. beta=sqrt(1-(2mp/sqrt(S))^2)
#define PI 3.14159265
#define Impact_parameter 6.0
#define alpha 0.0072993
#define Radius 6.38
#define Zp 79.0 // number of proton of heavy-ion
#define Outdensity 1 //0/1 denotes print out the E,B at a given point/any points
#define sigma_LQCD (5.8/197.0)  //sigma_v=5.8 Mev in LQCD (convert to fm^{-1} unit)
#define sigmaTc (10.0*sigma_LQCD) //with sigmaTc=10 sigma_LQCD
#define Model   212 //numberofmodel //firt 1(2)--Bjorken(Hubble) model, second 0(1 or 2) ---Homogenous(spehere or Gauss) case, third 1(2)---vacuum(turn-on) case //000 and 001 is for constant case
#define Medium_open 2 //1(2) the medium is corresponding to vacuum (turn-on)






