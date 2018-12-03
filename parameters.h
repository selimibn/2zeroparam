/*Definitons to run Mathematica CForms*/
#define Power(x, y)	(pow((double)(x), (double)(y)))
#define Sqrt(x)		(sqrt((double)(x)))
#define Cos(x)		(cos((double)(x)))
#define Sin(x)		(sin((double)(x)))
#define Abs(x) 		(fabs((double)(x)))
#define Pi acos(-1.)

/*Parameters taken from PDG, in GeV*/
#define GF 1.1663787E-5 
#define vev 246.0
#define alphaEM_MZ 127.937 

#define mlt 1776.86E-3
#define mle 0.510998928E-3
#define mlm 105.6583715E-3

#define mqu 2.3E-3
#define mqd 4.8E-3
#define mqc 1.275
#define mqs 95.0E-3
#define mqt 173.21
#define mqb 4.18

//Masses at 
#define MW 80.385
#define MZ 91.1876
#define MBs0 5366.79E-3
#define MB0 5279.61E-3
#define MBc0 6275.1E-3
#define MK0 497.611E-3
#define MD0 1864.84E-3
#define MPi0 134.9766E-3
#define Meta 547.862E-3
#define Metap 957.78E-3

//Lifetime
#define tBs0 (1.466E-12/6.58E-25)  //Using 1s = (1/6.58)E25 GeV^-1
#define tB0  (1.519E-12/6.58E-25)

//Form factors
#define f0Pi 130.1E-3
#define f0B 216.0E-3
#define f0Bs 227.0E-3
#define f0D 205.0E-3
#define f0K 155.0E-3

#define hetaq 2.0E-3 //GeV^3
#define hetas -53.0E-3 //GeV^3
#define hetapq 1.6E-3 //GeV^3
#define hetaps 65.0E-3 //GeV^3

//Norms of CKM matrix elements [PDG 2017]
#define Vud 0.97434
#define Vus 0.22506
#define Vub 0.00357
#define Vcd 0.22492
#define Vcs 0.97351
#define Vcb 0.0411
#define Vtd 0.00875
#define Vts 0.0403
#define Vtb 0.99915

//Error of CKM matrix elements [PDG 2017]
#define ErrVud 0.00011
#define ErrVus 0.00050
#define ErrVub 0.00015
#define ErrVcd 0.00050
#define ErrVcs 0.00013
#define ErrVcb 0.0013
#define ErrVtd 0.00033
#define ErrVts 0.0013
#define ErrVtb 0.00005

//Jarlskog Invariant [PDG 2017]
#define JS 3.04E-5
#define ErrJS 0.21E-5


//Phases of CKM elements
/*
#define th12 (13.04*(Pi/180)) //Radians
#define th13 (0.201*(Pi/180)) //Radians
#define th23 (2.38*(Pi/180)) //Radians
#define delta -2.0460 //CP Phase (Calculated from the gamma measurement in JHEP 10 (2014) 097, 10.1007/JHEP10(2014)097)
*/
#define th12 0.227314//Radians
#define th13 0.00355001//Radians
#define th23 0.0414121 //Radians
#define delta -2.0460 //CP Phase (Calculated from the gamma measurement in JHEP 10 (2014) 097, 10.1007/JHEP10(2014)097)
/*
#define phud 0.0
#define phus 0.0
#define phub (-delta)
#define phcd atan((cos(th12)*sin(th23)*sin(th13)*sin(delta))/(sin(th12)*cos(th23)+cos(th12)*sin(th23)*sin(th13)*cos(delta)))

#define phcs atan((sin(th12)*sin(th23)*sin(th13)*sin(delta))/(sin(th12)*sin(th23)*sin(th13)*cos(delta)-cos(th12)*cos(th23)))

#define phcb 0.0

#define phtd atan((cos(th12)*cos(th23)*cos(th13)*sin(delta))/(cos(th12)*cos(th23)*sin(th13)*cos(delta)-sin(th13)*sin(th23)))

#define phts atan((sin(th12)*cos(th23)*sin(th13)*sin(delta))/(cos(th12)*sin(th23)-sin(th12)*cos(th23)*sin(th13)*cos(delta)))

#define phtb 0.0
*/


#define phud 0.0
#define phus 0.0
#define phub 2.0460
#define phcd 3.14103

#define phcs 0.0000302548

#define phcb 0.0

#define phtd 0.274514

#define phts 3.12381

#define phtb 0.0


