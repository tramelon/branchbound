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

#define Tt 24
#define Nbus 118
#define T1 24	//time
#define Nn 118

void readData (const char* filename, IloNumArray& Runcost,IloInt I);
void readIntData (const char* filename, IloIntArray& Inta, IloInt N);
void readIntMaData (const char* filename, IntMatrix& shift, IloInt T, IloInt S);
void readMaData (const char* filename, NumMatrix& shift, IloInt T, IloInt S);
void readBoolData (const char* filename, IloBoolArray& thisbool, IloInt N);
void readMyData (const char* filename, int* shift[][Tt+1], int T, int S);

extern IloEnv env;
extern int Zd[Nbus][Tt+1];
extern double node_gap0; // gap to redeem a node to be optimal
extern double outgap0; // out loop tolerable gap
extern double outmastol;// out mas problem tolerable gap
extern double in_tlrn;// in gap tolerable, not used in branch and bound routine
extern double inmastol;// in mas tolerable, not used in calub(), optional
extern double insubtol;// above, in sub tolerable
extern double costuvw;
extern double ldpen;
extern double ldpena;
extern int presol; // presolve option, 0: disable, 1: use presolving
extern int solver_opt; // solver option, 0: default(Dynamic Search), 1: tradition B&C
extern double uncda;// uncertainty portion
extern int uncbg;
extern long J;	//units
extern long T;	//time
extern long N;    // buses
extern long L;    // trans lines
extern double sellprice;// sell price of elctricity

/*
extern IloInt uncdbug ;
extern IloNum uncd; // uncertainty part
extern double Md; // penalty for dual
extern IloInt it_out;
extern IloNum sellprice;
extern IloNum ldpen;
extern ThreeDNumMatrix Zoutstar;
extern ThreeDVarMatrix g;
extern ThreeDVarMatrix goverline;
extern ThreeDVarMatrix reserve;
extern ThreeDVarMatrix plt;
extern ThreeDVarMatrix dlnt;
extern ThreeDVarMatrix delta;
extern ThreeDBoolVarMatrix h; // start-up type
extern IloNumVar beta;
extern NumMatrix Tsu;     // Tsu_gs
extern NumMatrix Csu;     // Csu_gs
extern NumMatrix Kjt;
extern IloNumArray aj;    // no loas cost
extern IloNumArray rj;    //
extern IloNumArray cj;    // linear variable cost
extern IloNumArray Gj;    // min output
extern IloNumArray G0;    // Initial output
extern IloNumArray Gjb;    // max output
extern IloNumArray Dj; 	// number of hours j required to be off at the start period
extern IloNumArray Uj; 	// number of hours j required to be up at the start period
extern IloNumArray uj0;     //initial status of gen
extern IloNumArray gj0;     //initial status of gen
extern IloNumArray RU; 	// ramping up
extern IloNumArray RD; 	// ramping down
//class IloNumArray<RD>;
//typedef IloNumArray RD;
extern IloNumArray UT;     // min up
extern IloNumArray UT0;    // online time before schedule
extern IloNumArray UTr;    //
extern IloNumArray DT;     // min down
extern IloNumArray DT0;    // offline time before schedule
extern IloNumArray DTr;    //
extern IloNumArray SD;     // max shut-down rate
extern IloNumArray SU;     // max start-up rate
extern IloNumArray Sj;     // start-up segments	// for now all 1 segments
extern IloNumArray x;	// susceptance
extern IloNumArray Rt;	// spinning reserve at time t
extern IloNumArray Dt;	// total demand at time t
extern IloNumArray Pl;	// maximum flow on line
extern IloNumArray Gloc;	// generation location
extern IntMatrix line(env,L);      // trans lines
extern NumMatrix dnt;     // dnt
extern NumMatrix dpnt;    //dpnt
*/