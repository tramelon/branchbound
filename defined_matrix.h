#include <ilcplex/ilocplex.h>
#include <ctime>
#include <algorithm>
#include <cstdlib> 
#include <fstream>
#include <math.h>
#include <tchar.h>
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <cstring>
#include <string>
#include <iostream>

typedef IloArray<IloIntVarArray>    IntVarMatrix;
typedef IloArray<IloNumVarArray>    NumVarMatrix;
typedef IloArray<NumVarMatrix>  ThreeDVarMatrix;
typedef IloArray<ThreeDVarMatrix> FourDVarMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<NumMatrix>  ThreeDMatrix;
typedef IloArray<ThreeDMatrix> FourDMatrix;
typedef IloArray<IloBoolVarArray>    BoolVarMatrix;
typedef IloArray<BoolVarMatrix>      ThreeDBoolVarMatrix;
typedef IloArray<ThreeDBoolVarMatrix>      FourDBoolVarMatrix;
typedef IloArray<FourDBoolVarMatrix>      FiveDBoolVarMatrix;
typedef IloArray<IloBoolArray> BoolMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<NumMatrix> ThreeDNumMatrix;
typedef IloArray<IloIntArray> IntMatrix;
typedef IloArray<NumMatrix>  ThreeDMatrix;

IloEnv env;
IloInt J=54;	//units
IloInt T=24;	//time
IloInt N=118;    // buses
IloInt L=186;    // trans lines
IloNum costuvw; //cost of start up and fixed cost
IloNum Md; // penalty for dual
IloNumArray CarD(env,N); //cardinality for uncertain demand
IloNumArray duals(env,L);
IloNumArray Fbounds(env,L);
IloInt m; // O(l),origin or line
IloInt r; // D(l),destination of line
IloNum maxdelta= 1.57;
double obj_val;
IloNum scale=1; // nodal demand scale
IloNum ldpen=50; // penalty for load shed
IloInt uncdbug = 5;  // demand uncertain budget
IloNum uncd =0.1; // uncertainty part / demand
IloNum in_tlrn = 0.01 ; // in gap tolorence
IloInt itrc_out = 0;
IloInt res_itrc = 20;// initial set of secenarios in mas problem
IloNum outgap = 1;
IloNum outLB = -IloInfinity;
IloNum outUB = IloInfinity;
double runningtime;
IloNum sellprice = 0;
IloInt it_out = 20; // # of out iterations
ThreeDNumMatrix Zoutstar(env,it_out);
ThreeDVarMatrix g(env,it_out);
ThreeDVarMatrix goverline(env,it_out);
ThreeDVarMatrix reserve(env,it_out);
ThreeDVarMatrix plt(env,it_out);
ThreeDVarMatrix dlnt(env,it_out);
ThreeDVarMatrix delta(env,it_out);
ThreeDBoolVarMatrix h(env,J); // start-up type
IloNumVar beta(env, -IloInfinity,IloInfinity);
NumMatrix Tsu(env, J);     // Tsu_gs
NumMatrix Csu(env, J);     // Csu_gs
NumMatrix Kjt(env, J);     
IloNumArray aj(env,J);    // no loas cost
IloNumArray rj(env,J);    // 
IloNumArray cj(env,J);    // linear variable cost
IloNumArray Gj(env,J);    // min output
IloNumArray G0(env,J);    // Initial output
IloNumArray Gjb(env,J);    // max output
IloNumArray Dj(env,J); 	// number of hours j required to be off at the start period
IloNumArray Uj(env,J); 	// number of hours j required to be up at the start period
IloNumArray uj0(env,J);     //initial status of gen
IloNumArray gj0(env,J);     //initial status of gen
IloNumArray RU(env,J); 	// ramping up
IloNumArray RD(env,J); 	// ramping down
IloNumArray UT(env,J);     // min up
IloNumArray UT0(env,J);    // online time before schedule
IloNumArray UTr(env,J);    // 
IloNumArray DT(env,J);     // min down
IloNumArray DT0(env,J);    // offline time before schedule
IloNumArray DTr(env,J);    // 
IloNumArray SD(env,J);     // max shut-down rate
IloNumArray SU(env,J);     // max start-up rate
IloNumArray Sj(env,J);     // start-up segments	// for now all 1 segments
	//  cout<<"Sj"<<Sj<<endl;
	//	IloNumArray et(env,T);     // purchase price
	//	IloNumArray qt(env,T);     // sale price
	//	IloNumArray dtmin(env,T);  // min demand
	//	IloNumArray dtmax(env,T);  // max demand
IloNumArray x(env,L);	// susceptance
IloNumArray Rt(env,T+1);	// spinning reserve at time t
IloNumArray Dt(env,T+1);	// total demand at time t
IloNumArray Pl(env,L);	// maximum flow on line
IloNumArray Gloc(env,J) ;	// generation location
//-------------unc demand
	//	NumVarMatrix  ucd(env,N);
	//	BoolVarMatrix card(env,N);
IntMatrix line(env,L);      // trans lines
NumMatrix dnt(env, N);     // dnt
NumMatrix dpnt(env, N);    //dpnt
IloInt lazcount=0;


