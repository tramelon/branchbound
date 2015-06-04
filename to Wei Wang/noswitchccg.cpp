// ccg for sub prblem

#ifndef DEFINEDMATRIX_H
#define DEFINEDMATRIX_H
#include "definedmatrix.h"
#endif DEFINEDMATRIX_H
#include "ilconcert\ilosys.h"
//#ifndef PARINIT_H
//#define PARINIT_H
//#include "parinit.h"
//#endif PARINIT_H
//#include "defined_matrix.h"
#include <iostream>
using namespace std;

int noswitchccg(double &UB, double &LB, int U[][T1+1], int V[][T1+1], int W[][T1+1],float Znode[][T1+1], int budget, double costuvw)
	{
	IloInt m; // O(l),origin or line
	IloInt r; // D(l),destination of line
	IloNum scale=1; // nodal demand scale
	IloNum maxdelta= 1.57;
	double obj_val;
	double ldpen=50; // penalty for load shed
	IloInt uncdbug = budget;  // demand uncertain budget
	double uncd = uncda; // uncertainty part / demand
	double in_tlrn = 0.005 ; // in gap torrence, 0.01
	double sol_gap = 0.005 ; // solver gap
	IloInt itrc_out = 0;
	IloInt res_itrc = 20;// initial set of scenarios in mas problem
	IloNum outgap = 1;
	IloNum outLB = -IloInfinity;
	IloNum outUB =  IloInfinity;
	double runningtime;
//	IloInt J=54;	//units
//	IloInt T=T1;	//time
//	IloInt N=Nn;    // buses
//	IloInt L=186;    // trans lines
	IloEnv env;

// NumMatrix dnt(env);
		IloNum Md; // penalty for dual
		IloNum sellprice = 0;
		IloInt it_out = 20; // # of out iterations
		ThreeDNumMatrix Zoutstar(env,it_out);
		 {for(IloInt it = 0; it < it_out; it++)
			   Zoutstar[it] = NumMatrix(env,N);
		  for(IloInt it = 0; it < it_out; it++)
				for(IloInt n = 0; n < N; n++)
					Zoutstar[it][n] = IloNumArray(env,T+1); }
		 for(IloInt it = 0; it < 2; it++)
			 for(IloInt n = 0; n < N; n++)
				 for(IloInt t = 0; t <= T; t++)
					 Zoutstar[it][n][t] = 0;
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
		IloNumArray aj(env,J);    // no loads cost
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
		NumMatrix dpnt(env, N);    //dpn
		IloInt lazcount=0;

			for(IloInt n=0;n<N;n++){
				dnt[n]=IloNumArray(env,T+1);
			}
			for(IloInt n=0;n<N;n++){
				dpnt[n]=IloNumArray(env,T+1);
			}
// kjt
			for(IloInt j=0;j<J;j++){
				Kjt[j]=IloNumArray(env,T+1);
			}
		NumMatrix Ustar(env, J);     // Ustar,output
			for(IloInt j=0;j<J;j++){
				Ustar[j]=IloNumArray(env,T+1);
			}
//		IntMatrix line(env,L);      // trans lines
		   for(IloInt l = 0; l < L; l++){
				line[l] = IloIntArray(env,2);
			}
		IloInt Betn = 6;
		IntMatrix cut(env, Betn);      // cut trans lines
		   for(IloInt l = 0; l < Betn; l++){
				cut[l] = IloIntArray(env,2);
			}
		IloBoolArray label(env, N);
		IloInt Betn1 = 4;
		IntMatrix cut1(env, Betn1);      // cut trans lines
		   for(IloInt l = 0; l < Betn1; l++){
				cut1[l] = IloIntArray(env,2);
			}
		IloBoolArray label1(env, N);
		IloInt Betn2 = 12;
		IntMatrix cut2(env, Betn2);      // cut trans lines
		   for(IloInt l = 0; l < Betn2; l++){
				cut2[l] = IloIntArray(env,2);
			}
		IloIntArray label2(env, N);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/Sj.txt", Sj, J);
			for(IloInt j=0;j<J;j++){
				Csu[j]=IloNumArray(env,Sj[j]);
			}
			for(IloInt j=0;j<J;j++){
				Tsu[j]=IloNumArray(env,Sj[j]);
			}
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/aj.txt", aj, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/rj.txt", rj, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/cj.txt", cj, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/Gj.txt", Gj, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/uj0.txt", uj0, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/G0.txt", G0, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/Gjb.txt", Gjb, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/Dj.txt", Dj, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/Uj.txt", Uj, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/RU.txt", RU, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/RD.txt", RD, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/UT.txt", UT, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/DT.txt", DT, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/UT0.txt", UT0, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/DT0.txt", DT0, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/SD.txt", SD, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/SU.txt", SU, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/Sj.txt", Sj, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/Gloc.txt", Gloc, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/Rt.txt", Rt, T+1);
		readMaData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/Csu.txt", Csu, J, 2);
		readMaData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/Tsu.txt", Tsu, J, 2);
		readMaData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/dnt.txt", dnt, N, T+1);
		readIntMaData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/transmission.txt", line, L, 2);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/flowlimit.txt",Pl, L);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/reactance.txt", x, L);
		Md=0;
		for(IloInt t = 1; t <= T; t++){
			Dt[t]=0;
			for(IloInt n = 0; n < N; n++){
				Dt[t]+=dnt[n][t];
				Md+=dnt[n][t];
				}
			}
		Md=IloMax(Dt);
		for(IloInt t = 0; t <= T; t++)
			for(IloInt n = 0; n < N; n++)
				dpnt[n][t]= uncd*dnt[n][t];

ofstream results("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/results1.txt");
ofstream Duals("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/duals1.txt");
ofstream flows("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/flows1.txt");
ofstream ccgresults("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/ccg1.txt",ios::app);

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ overall Sub problem $$$$$$$$$$$$$$$$$$$$$$$$$$//
	NumMatrix Zntstar(env,N);
	for(IloInt n=0;n<N;n++){Zntstar[n]=IloNumArray(env,T+1);}
	for(IloInt n=0;n<N;n++)
	  for(IloInt t = 1; t <= T; t++)
		  Zntstar[n][t]=0;
	  IloInt itrc = 20;
	ThreeDNumMatrix gcstar(env,itrc);
		 { for(IloInt it = 0; it < itrc; it++)
			  gcstar[it] = NumMatrix(env,J);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt j = 0; j < J; j++)
					gcstar[it][j] = IloNumArray(env,T+1); }
	NumMatrix this_g(env,J);
		{for(IloInt j=0;j<J;j++){this_g[j]=IloNumArray(env,T+1);}}
// ccg master problem variable
	  BoolVarMatrix znt(env,N);
	  for(IloInt n=0;n<N;n++){znt[n]=IloBoolVarArray(env,T+1);}
	  ThreeDVarMatrix pi5c(env,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  pi5c[it] = NumVarMatrix(env,L);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt l = 0; l < L; l++)
					pi5c[it][l] = IloNumVarArray(env,T+1, -IloInfinity, IloInfinity);
	  ThreeDVarMatrix pi6c(env,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  pi6c[it] = NumVarMatrix(env,N);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt n = 0; n < N; n++)
					pi6c[it][n] = IloNumVarArray(env,T+1, 0, IloInfinity);
	  ThreeDVarMatrix pi6ac(env,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  pi6ac[it] = NumVarMatrix(env,N);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt n = 0; n < N; n++)
					pi6ac[it][n] = IloNumVarArray(env,T+1, 0, IloInfinity);
	  ThreeDVarMatrix pi7c(env,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  pi7c[it] = NumVarMatrix(env,L);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt l = 0; l < L; l++)
					pi7c[it][l] = IloNumVarArray(env,T+1, 0, IloInfinity);
	 ThreeDVarMatrix pi8c(env,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  pi8c[it] = NumVarMatrix(env,L);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt l = 0; l < L; l++)
					pi8c[it][l] = IloNumVarArray(env,T+1, 0, IloInfinity);
	ThreeDVarMatrix pi9c(env,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  pi9c[it] = NumVarMatrix(env,N);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt n = 0; n < N; n++)
					pi9c[it][n] = IloNumVarArray(env,T+1, 0, IloInfinity);
	ThreeDVarMatrix pi10c(env,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  pi10c[it] = NumVarMatrix(env,N);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt n = 0; n < N; n++)
					pi10c[it][n] = IloNumVarArray(env,T+1, 0, IloInfinity);
	ThreeDVarMatrix pi11c(env,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  pi11c[it] = NumVarMatrix(env,N);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt n = 0; n < N; n++)
					pi11c[it][n] = IloNumVarArray(env,T+1, 0, IloInfinity);
		  ThreeDVarMatrix beta1c(env,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  beta1c[it] = NumVarMatrix(env,N);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt n = 0; n < N; n++)
					beta1c[it][n] = IloNumVarArray(env,T+1, 0, IloInfinity);
		  ThreeDVarMatrix beta2c(env,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  beta2c[it] = NumVarMatrix(env,N);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt n = 0; n < N; n++)
					beta2c[it][n] = IloNumVarArray(env,T+1, 0, IloInfinity);
		  ThreeDVarMatrix beta3c(env,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  beta3c[it] = NumVarMatrix(env,N);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt n = 0; n < N; n++)
					beta3c[it][n] = IloNumVarArray(env,T+1, 0, IloInfinity);

	IloNum lb_sub = -IloInfinity;
	IloNum ub_sub = IloInfinity;
	IloNum ingap = 1;
	itrc=0;
//***********ccg loop for sub problem****//
try{
time_t begin,end;
begin = clock();
	while(ingap>=in_tlrn)
	 {
	 itrc++;
//	 gcstar.add(this_g);
//***************************CCG sub problem*******************************//
	NumVarMatrix gc(env,J);  // generation level, gjt
		for(IloInt j = 0; j < J; j++)
			gc[j] = IloNumVarArray(env,T+1, 0, IloInfinity);
	NumVarMatrix pltc(env,L);  // power flow, plt
		for(IloInt l = 0; l < L; l++)
			pltc[l] = IloNumVarArray(env,T+1, -IloInfinity, IloInfinity);
	NumVarMatrix dlntc(env,N);  // load shed, dlnt
		for(IloInt n = 0; n < N; n++)
				dlntc[n] = IloNumVarArray(env,T+1, 0, IloInfinity);
	for(IloInt n = 0; n < N; n++)
			  for(IloInt t = 1; t <= T; t++)
				dlntc[n][t].setBounds(0,dnt[n][t]);
	NumVarMatrix deltac(env,N);  // phase angles, delta0
		for(IloInt n = 0; n < N; n++)
				deltac[n] = IloNumVarArray(env,T+1, -maxdelta, maxdelta);
IloModel subpc(env);
IloExpr subobjc(env);
   for(IloInt j=0; j<J; j++)
		for(IloInt t=1; t<=T; t++){
			subobjc+=cj[j]*gc[j][t];
		}
   for(IloInt t=1; t<=T; t++)
	   for(IloInt n= 0; n < N; n++)
		   subobjc+=ldpen*dlntc[n][t];
	subpc.add(IloMinimize(env,subobjc));
	subobjc.end();
//** constraints**//
for(IloInt j=0; j<J; j++){
	subpc.add(gc[j][0]==0);
	}
for(IloInt j=0; j<J; j++)
	for(IloInt t=1; t<=T; t++){
		subpc.add(gc[j][t-1]-gc[j][t] <= RD[j]*U[j][t]+SD[j]*W[j][t]); //Ramping down----traditional
		}
for(IloInt j=0; j<J; j++)
	for(IloInt t=1; t<=T; t++){
		subpc.add(gc[j][t]-gc[j][t-1] <= RU[j]*U[j][t-1]+SU[j]*V[j][t]); //Ramping Up----traditional
		}
for(IloInt j=0; j<J; j++)
	for(IloInt t= 1; t <= T; t++){
			subpc.add(gc[j][t]>=U[j][t]*Gj[j]); //......Generation lower bound limit
			subpc.add(gc[j][t]<=U[j][t]*Gjb[j]);
			}
//......Transmission Network Constraints.
	for(IloInt l = 0; l < L; l++)
		 for(IloInt t = 1; t <= T; t++)  {
				pltc[l][t].setBounds(-Pl[l],Pl[l]);
				subpc.add(pltc[l][t]>=-Pl[l]);
				subpc.add(pltc[l][t]<=Pl[l]);
		 }
//......for L lines, DC power flow.
  {for(IloInt t= 1; t <= T; t++)
	 for(IloInt l=0; l<L;l++)
	 {
		 m=line[l][0]; //origin bus
		 r=line[l][1]; //destination bus
		 subpc.add(pltc[l][t]*x[l]==100*(deltac[m][t]-deltac[r][t]));
	 }}
  {for(IloInt t= 1; t <= T; t++)
	 for(IloInt n=0; n<N;n++)
	 {
		 subpc.add(dlntc[n][t]<=dnt[n][t]+Zntstar[n][t]*dpnt[n][t]);
		 subpc.add(deltac[n][t]<=maxdelta);
		 subpc.add(deltac[n][t]>=-maxdelta);
	 }}
//......load Balancing constraints, for each bus.
   for(IloInt t= 1; t <= T; t++)
	for(IloInt n= 0; n < N; n++)
		{IloExpr exp60c(env);
		 IloExpr exp70c(env);
		 IloExpr exp80c(env);
		{for(IloInt i= 0; i < J; i++)
		{
			if(Gloc[i]==n){
				exp60c+=gc[i][t];
				}
		}}
	   {for(IloInt l = 0; l < L; l++)
		{
				if(line[l][1]==n)  // D(l)
					exp70c+=pltc[l][t];
				else if(line[l][0]==n) // O(l)
					exp80c+=pltc[l][t];
	   }}
		subpc.add(exp60c+exp70c-exp80c+dlntc[n][t] >= dnt[n][t]+Zntstar[n][t]*dpnt[n][t]); //load shed allowed
		exp60c.end();
		exp70c.end();
		exp80c.end();
	}//end n
// generation greater than demand
   for(IloInt t=1; t<=T; t++) {
		IloExpr expsr1(env);
		for(IloInt j=0; j<J; j++){
			expsr1+=gc[j][t];
			}
	//	subpc.add(expsr1 >= Dt[t]);
		expsr1.end();
		}

	IloCplex ccg_sub(subpc);
//********************CCG sub problem***************
//	ccg_sub.exportModel("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/uc_ccg_sub.lp") ;
	ccg_sub.setParam(IloCplex::MIPSearch,solver_opt);  //1, traditional branch & cut
	ccg_sub.setParam(IloCplex::EpGap, sol_gap);
	if(presol==0){
	ccg_sub.setParam(IloCplex::PreInd,false);
		}
	if (!ccg_sub.solve()) {
		env.error() << "Failed to solve CCG sub model__original form" << endl;
		throw(-1);
		}
	cout<<"ccg sub obj: "<<ccg_sub.getObjValue()<<endl;
	cout<<"obj: "<<costuvw+ccg_sub.getObjValue()<<endl;
//	cout<<" diff: "<<obj_val-costuvw- ccg_sub.getObjValue()<<endl;
//update
	for(IloInt t = 1; t <= T; t++)
		for(IloInt j= 0; j < J; j++){
			gcstar[itrc][j][t]=ccg_sub.getValue(gc[j][t]);
//			cout<<""<<gcstar[j][t]<<endl;
			}
	lb_sub = IloMax(lb_sub,ccg_sub.getObjValue()); // best feasible for LB
	ingap = (ub_sub -lb_sub)/ub_sub ;
	cout<<"ingap: "<<ingap<<endl;
	ccg_sub.end();
	subpc.end();

//**********************************************************************//
//*********** seperate  ccg master problem over time********************//
//**********************************************************************//
IloNumArray ccgmasobjv(env,T+1); // record obj func val
{for(IloInt t = 1; t <=T; t++)
	{
	IloNumVar eta(env,-IloInfinity,IloInfinity);
	IloModel maspc(env);
	//---------constraints---------------------------//
	//------- linearization constraints
	Md=ldpen;
	{for(IloInt n = 0; n < N; n++)	//dnt
	for(IloInt it=1; it<=itrc; it++){
			maspc.add(beta1c[it][n][t] <= pi6c[it][n][t]);
			maspc.add(beta2c[it][n][t] <= pi6ac[it][n][t]);
			maspc.add(beta3c[it][n][t] <= pi11c[it][n][t]);
			maspc.add(beta1c[it][n][t] <= Md*znt[n][t]);
			maspc.add(beta2c[it][n][t] <= Md*znt[n][t]);
			maspc.add(beta3c[it][n][t] <= Md*znt[n][t]);
			maspc.add(beta1c[it][n][t] >= Md*znt[n][t]+pi6c[it][n][t]-Md);
			maspc.add(beta2c[it][n][t] >= Md*znt[n][t]+pi6ac[it][n][t]-Md);
			maspc.add(beta3c[it][n][t] >= Md*znt[n][t]+pi11c[it][n][t]-Md);
			}}
//--------- dual constraints
		IloInt lm;
		IloInt lr;
	{for(IloInt l = 0; l < L; l++) // plt
		for(IloInt it=1; it<=itrc; it++){
			lm=line[l][0];
			lr=line[l][1];
			maspc.add((x[l]/100)*pi5c[it][l][t]-pi6c[it][lm][t]+pi6ac[it][lm][t]
			+pi6c[it][lr][t]-pi6ac[it][lr][t]+pi7c[it][l][t]-pi8c[it][l][t] == 0);
			}}
	{for(IloInt n = 0; n < N; n++)	// delta nt
		for(IloInt it = 1; it <= itrc; it++){
		IloExpr isexpr1(env);
		IloExpr isexpr2(env);
		 for(IloInt l = 0; l < L; l++)
		{
			if(line[l][1]==n)  // D(l)
			{
					isexpr1+= pi5c[it][l][t];
			}//end if
			if(line[l][0]==n)  // O(l)
			{
					isexpr2+= -pi5c[it][l][t];
			} //end if
		}//end for l
		maspc.add(isexpr2+isexpr1+pi9c[it][n][t]-pi10c[it][n][t]==0);
		isexpr1.end();
		isexpr2.end();
			}}
	{for(IloInt n = 0; n < N; n++)	//dnt
		{for(IloInt it = 1; it <= itrc; it++){             /*********************pi6ac[it][lr][t]*****************************/
			maspc.add(pi6c[it][n][t]-pi6ac[it][n][t]-pi11c[it][n][t]<=ldpen);
	}}}
	//--------------obj func
{for(IloInt it=1; it<=itrc; it++)
	{
	IloExpr submasobj(env);
	{for(IloInt l = 0; l < L; l++){
		submasobj += (-Pl[l]*(pi7c[it][l][t]+pi8c[it][l][t]));
			}}
	{for(IloInt n = 0; n < N; n++){
		submasobj += -maxdelta*(pi9c[it][n][t]+pi10c[it][n][t]);
		submasobj +=  dnt[n][t]*(pi6c[it][n][t]-pi6ac[it][n][t]-pi11c[it][n][t]);
		submasobj +=  dpnt[n][t]*(beta1c[it][n][t]-beta2c[it][n][t]-beta3c[it][n][t]);
			}}
		IloNum locgc;
	{for(IloInt n = 0; n < N; n++){
			locgc = 0 ;
			for(IloInt i= 0; i < J; i++){
				if(Gloc[i]==n)
					locgc+=gcstar[it][i][t];}
			submasobj += -locgc*(pi6c[it][n][t]-pi6ac[it][n][t]);
			}}
	IloNum totalgstar = 0;
	{for(IloInt j= 0; j < J; j++){
			totalgstar += cj[j]*gcstar[it][j][t];
			}}
//	IloRange ccgrange(env,-IloInfinity, eta-submasobj,totalgstar);
//	maspc.add(ccgrange);
	maspc.add(eta<=submasobj+totalgstar);
	submasobj.end();
//	ccgrange.end();
}}//end for it
	IloObjective ccgmobj = IloMaximize(env, eta);
	maspc.add(ccgmobj);

	IloExpr bugc(env);
	   {for(IloInt n = 0; n < N; n++){
			bugc += znt[n][t];}
		maspc.add(bugc<=uncdbug);} // budget on demand uncertainty
		bugc.end();

	IloCplex ccg_maspc(maspc);
//	ccg_maspc.exportModel("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/uc_ccgma.lp") ;
	cout<<"constraints and columns added "<<endl;
	ccg_maspc.setParam(IloCplex::MIPSearch,solver_opt);  //1, traditional branch & cut
		if(presol==0){
	ccg_maspc.setParam(IloCplex::PreInd,false);
		}
	ccg_maspc.setParam(IloCplex::EpGap, sol_gap); // inner master problem
	if (!ccg_maspc.solve()) {
		env.error() << "Failed to solve ccg mas model_original form" << endl;
		throw(-1);
		}
//****************update the uncertainty scenario***************
		{for(IloInt n= 0; n < N; n++){
			Zntstar[n][t]=ccg_maspc.getValue(znt[n][t]);
			Znode[n][t]=ccg_maspc.getValue(znt[n][t]);
		}}
	ccgmasobjv[t] = ccg_maspc.getObjValue();
	ccg_maspc.end();
//	maspc.remove(ccgmobj);
	ccgmobj.end();
	eta.end();
	maspc.end();
}
}// END for t, decomposed

	ub_sub = 0;
	{for(IloInt t = 1; t < T+1; t++){
		ub_sub += ccgmasobjv[t];
	}}
	ingap = (ub_sub -lb_sub )/ub_sub ;
	cout<<"ccg mas obj: "<<ub_sub<<endl;
	cout<<"obj: "<<costuvw+ub_sub<<endl;
	cout<<"ingap: "<<ingap<<endl;
	end=clock();
	runningtime=double(end-begin)/CLOCKS_PER_SEC;
	results<<"iteration: "<<itrc<<" ccg_lb: "<<lb_sub<<" ccg_ub: "<<ub_sub<<" gap: "<<ingap<<" time: "<<runningtime<<endl;
 }//end while

outUB = IloMin(costuvw+lb_sub,outUB);
outgap = (outUB-outLB)/outUB;
ccgresults<<"---outUB: "<<outUB<<" outGap: "<<outgap <<endl;
cout<<"good til"<<endl;
}// end try
	catch (IloException& ex){	  cerr << "Error: " << ex << endl;	}
	catch (...)				{	  cerr << "Error" << endl;			}
	env.end();
//update bounds
UB = ub_sub;
LB = IloMax(lb_sub,LB);
cout<<"Node info: "<<"UB--"<<UB<<" LB--"<<LB<<endl;
cout<<itrc<<endl;
return 0;
}//end ubcal