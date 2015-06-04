


//#ifndef DEFINEDMATRIX_H
//#define DEFINEDMATRIX_H
#include "definedmatrix.h"
//#endif DEFINEDMATRIX_H
using namespace std;

int Zd[Nbus][Tt+1];



IloEnv env;
double Md;
double ldpen = ldpena;
//IloInt it_out = 20;
//IloInt itrc_out = 0;
//ThreeDNumMatrix Zoutstar(env,it_out);
//ThreeDVarMatrix g(env,it_out);
//ThreeDVarMatrix goverline(env,it_out);
//ThreeDVarMatrix reserve(env,it_out);
//ThreeDVarMatrix plt(env,it_out);
//ThreeDVarMatrix dlnt(env,it_out);
//ThreeDVarMatrix delta(env,it_out);
//ThreeDBoolVarMatrix h(env,J); // start-up type
IloNumVar beta(env, -IloInfinity,IloInfinity);
NumMatrix Tsu(env, J);     // Tsu_gs
NumMatrix Csu(env, J);     // Csu_gs
NumMatrix Kjt(env, J);
IloNumArray aj(env,J);    // no load cost
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
IloNumArray x(env,L);	// susceptance
IloNumArray Rt(env,T+1);	// spinning reserve at time t
IloNumArray Dt(env,T+1);	// total demand at time t
IloNumArray Pl(env,L);	// maximum flow on line
IloNumArray Gloc(env,J) ;	// generation location
//IntMatrix line;      // trans lines
NumMatrix dnt(env, N);     // dnt
NumMatrix dpnt(env, N);    //dpnt

void readData (const char* filename, IloNumArray& Runcost,IloInt I)
{
   ifstream in(filename);
   for(IloInt i=0;i<I;i++)
	   in >> Runcost[i];
}
 void readIntData (const char* filename, IloIntArray& Inta, IloInt N)
{
   ifstream in(filename);
   for(IloInt i=0;i<N;i++)
	   in >> Inta[i];
}
void readBoolData (const char* filename, IloBoolArray& thisbool, IloInt N)
{
   ifstream in(filename);
   for(IloInt i=0;i<N;i++)
	   in >> thisbool[i];
}
void readMaData (const char* filename, NumMatrix& shift, IloInt T, IloInt S)
{
   ifstream in(filename);
   for(IloInt t=0;t<T;t++)
	   for (IloInt s=0;s<S;s++)
	   in >> shift[t][s];
}
 void readIntMaData (const char* filename, IntMatrix& shift1, IloInt T, IloInt S)
{
   ifstream in(filename);
   for(IloInt t=0;t<T;t++)
	   for (IloInt s=0;s<S;s++)
	   in >> shift1[t][s];
}
 void readMyData (const char* filename, int* shift[][25], int T, int S)
{
   ifstream in(filename);
   for(int t=0;t<T;t++)
	   for (int s=0;s<S;s++)
	   in >> *shift[t][s];
}