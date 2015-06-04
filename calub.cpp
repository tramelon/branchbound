//control unc budget line 47 uncdbug



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
using namespace std;

int calub(double &UB, double &LB, int Z[][T1+1], int U[][T1+1], int V[][T1+1], int W[][T1+1],float Znode[][T1+1], int budget, double costuvw)
	{
//	IloInt J=54;	//units
//	IloInt T=T1;	//time
//	IloInt N=Nn;    // buses
//	IloInt L=186;    // trans lines
IloEnv lazenv;
   NumMatrix laz_u(lazenv, J);
   NumMatrix laz_v(lazenv, J);
   NumMatrix laz_w(lazenv, J);
   for(IloInt j = 0; j < J; j++){
			laz_u[j] = IloNumArray(lazenv,T+1);
			laz_v[j] = IloNumArray(lazenv,T+1);
			laz_w[j] = IloNumArray(lazenv,T+1);
		}
   for(IloInt t = 0; t <= T; t++)
		for(IloInt j = 0; j < J; j++){
			laz_u[j][t] = U[j][t];  // pass start-up variables
			laz_v[j][t] = V[j][t];
			laz_w[j][t] = W[j][t];}

	IloInt m; // O(l),origin or line
	IloInt r; // D(l),destination of line
	IloNum scale=1; // nodal demand scale
	IloNum maxdelta= 1.57;
	double obj_val;
	IloNum ldpen=ldpena; // penalty for load shed
	IloInt uncdbug = budget;  // demand uncertain budget
	IloNum uncd = uncda; // uncertainty part / demand
//	double in_tlrn =  in_tlrn;//0.001 ; // in gap torrence,0.01
	double sol_gap = 0.005 ; // solver gap
	IloInt itrc_out = 0;
	IloInt res_itrc = 20;// initial set of scenarios in mas problem
	IloNum outgap = 1;
	IloNum outLB = -IloInfinity;
	IloNum outUB =  IloInfinity;
	double runningtime;
// NumMatrix dnt(env);
//		IloNum costuvw; //cost of start up and fixed cost
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
		NumMatrix dpnt(env, N);    //dpnt
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
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/Sj.txt", Sj, J);
			for(IloInt j=0;j<J;j++){
				Csu[j]=IloNumArray(env,Sj[j]);
			}
			for(IloInt j=0;j<J;j++){
				Tsu[j]=IloNumArray(env,Sj[j]);
			}
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/rj.txt", rj, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/cj.txt", cj, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/Gj.txt", Gj, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/uj0.txt", uj0, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/G0.txt", G0, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/Gjb.txt", Gjb, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/Dj.txt", Dj, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/Uj.txt", Uj, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/RU.txt", RU, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/RD.txt", RD, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/UT.txt", UT, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/DT.txt", DT, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/UT0.txt", UT0, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/DT0.txt", DT0, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/SD.txt", SD, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/SU.txt", SU, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/Sj.txt", Sj, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/Gloc.txt", Gloc, J);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/Rt.txt", Rt, T+1);
		readMaData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/Csu.txt", Csu, J, 2);
		readMaData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/Tsu.txt", Tsu, J, 2);
		readMaData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/dnt.txt", dnt, N, T+1);
		readIntMaData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/transmission.txt", line, L, 2);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/flowlimit.txt",Pl, L);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters1/reactance.txt", x, L);
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
cout<<"parameter initialization done"<<endl;
ofstream results("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/results2.txt");
ofstream Duals("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/duals2.txt");
ofstream flows("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/flows2.txt");
ofstream ccgresults("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/ccg2.txt",ios::app);

//$$$$$$$$$$$$$$$$$$$$ Overall Sub problem $$$$$$$$$$$$$$$$$$$$$$//
   NumMatrix gcstar(lazenv,J);
	for(IloInt j=0;j<J;j++){gcstar[j]=IloNumArray(lazenv,T+1);}
	for(IloInt j=0;j<J;j++)
	  for(IloInt t = 0; t <= T; t++)
		  gcstar[j][t]=0;
	IloInt itrc = 10;
	ThreeDNumMatrix Zntstar(lazenv,itrc);
		{ for(IloInt it = 0; it < itrc; it++)
			   Zntstar[it] = NumMatrix(lazenv,N);
		 for(IloInt it = 0; it < itrc; it++)
				for(IloInt n = 0; n < N; n++)
					Zntstar[it][n] = IloNumArray(lazenv,T+1); }
	for(IloInt it = 0; it < 2; it++)
		for(IloInt n = 0; n < N; n++)
			for(IloInt t = 0; t <= T; t++)
				Zntstar[it][n][t] = 0;
// ccg master problem variable
	NumVarMatrix gc(lazenv,J);
	  for(IloInt j=0;j<J;j++){gc[j]=IloNumVarArray(lazenv,T+1,0,IloInfinity);}
	ThreeDVarMatrix pc(lazenv,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  pc[it] = NumVarMatrix(lazenv,L);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt l = 0; l < L; l++)
					pc[it][l] = IloNumVarArray(lazenv,T+1, -IloInfinity, IloInfinity);
	ThreeDVarMatrix deltac(lazenv,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  deltac[it] = NumVarMatrix(lazenv,N);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt n = 0; n < N; n++)
					deltac[it][n] = IloNumVarArray(lazenv,T+1, -IloInfinity, IloInfinity);
	ThreeDVarMatrix dc(lazenv,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  dc[it] = NumVarMatrix(lazenv,N);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt n = 0; n < N; n++)
					dc[it][n] = IloNumVarArray(lazenv,T+1, 0, IloInfinity);
	ThreeDVarMatrix sc(lazenv,itrc);
		  for(IloInt it = 0; it < itrc; it++)
			  sc[it] = NumVarMatrix(lazenv,N);
		  for(IloInt it = 0; it < itrc; it++)
				for(IloInt n = 0; n < N; n++)
					sc[it][n] = IloNumVarArray(lazenv,T+1, 0, IloInfinity);
	IloNum lb_sub = -IloInfinity;
	IloNum ub_sub =  IloInfinity;
	IloNum ingap = 1;
	IloNumVar eta(lazenv,-IloInfinity,IloInfinity);
	IloModel maspc(lazenv);
	maspc.add(IloMinimize(lazenv,eta));
//** ccg master cconstraints**
for(IloInt j=0; j<J; j++){
	maspc.add(gc[j][0]==0);
	}
for(IloInt j=0; j<J; j++)
	for(IloInt t=1; t<=T; t++){
		maspc.add(gc[j][t-1]-gc[j][t] <= RD[j]*laz_u[j][t]+SD[j]*laz_w[j][t]); //Ramping down----traditional
		}
for(IloInt j=0; j<J; j++)
	for(IloInt t=1; t<=T; t++){
		maspc.add(gc[j][t]-gc[j][t-1] <= RU[j]*laz_u[j][t-1]+SU[j]*laz_v[j][t]); //Ramping Up----traditional
		}
for(IloInt j=0; j<J; j++)
	for(IloInt t= 1; t <= T; t++){
			maspc.add(gc[j][t]>=laz_u[j][t]*Gj[j]); //......Generation lower bound limit
			maspc.add(gc[j][t]<=laz_u[j][t]*Gjb[j]);
			}
try {
	IloCplex ccg_mas1(lazenv);
	ccg_mas1.extract(maspc);
	ccg_mas1.exportModel("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/uc_ccg_mas1.lp") ;
	cout<<"model exported"<<endl;
	ccg_mas1.end();
//*********************************
itrc=0;
ingap=1;
time_t begin,end;
begin = clock();
	while(ingap>=in_tlrn)
	 {
	 itrc++;
//......Transmission Network Constraints.
	for(IloInt l = 0; l < L; l++)
		 for(IloInt t = 1; t <= T; t++)  {
				maspc.add(pc[itrc][l][t]>=-Pl[l]);
				maspc.add(pc[itrc][l][t]<=Pl[l]);
		 }
  {for(IloInt t= 1; t <= T; t++)
	 for(IloInt n=0; n< N;n++)
	 {
		 maspc.add(dc[itrc][n][t]<=dnt[n][t]+Zntstar[itrc][n][t]*dpnt[n][t]);
		 maspc.add(deltac[itrc][n][t]<=maxdelta);
		 maspc.add(deltac[itrc][n][t]>=-maxdelta);
	 }}
  //......For L lines, DC power flow.
  { for(IloInt l=0; l<L;l++){
	  int m=line[l][0]; //origin bus
	  int r=line[l][1]; //destination bus
	  for(IloInt t= 1; t <= T; t++)
		 {
		 maspc.add(pc[itrc][l][t]*x[l]==100*(deltac[itrc][line[l][0]][t]-deltac[itrc][line[l][1]][t]));
		  }}}
//......load Balancing constraints, for each bus.
   for(IloInt t= 1; t <= T; t++)
	for(IloInt n= 0; n < N; n++)
		{IloExpr exp60c(lazenv);
		 IloExpr exp70c(lazenv);
		 IloExpr exp80c(lazenv);
		{for(IloInt i= 0; i < J; i++)
		{
			if(Gloc[i]==n){
				exp60c+=gc[i][t];
				}
		}}
	   {for(IloInt l = 0; l < L; l++)
		{
				if(line[l][1]==n)  // D(l)
					exp70c+=pc[itrc][l][t];
				else if(line[l][0]==n) // O(l)
					exp80c+=pc[itrc][l][t];
	   }}
		maspc.add(exp60c+exp70c-exp80c+dc[itrc][n][t]-sc[itrc][n][t] >= dnt[n][t]+Zntstar[itrc][n][t]*dpnt[n][t]); //load shed allowed
		maspc.add(sc[itrc][n][t] <= exp60c);
		maspc.add(dc[itrc][n][t] <= dnt[n][t]+Zntstar[itrc][n][t]*dpnt[n][t]);
		exp60c.end();
		exp70c.end();
		exp80c.end();
	}//end n
// generation greater than demand
   for(IloInt t=1; t<=T; t++) {
		IloExpr expsr1(lazenv);
		for(IloInt j=0; j<J; j++){
			expsr1+=gc[j][t];
			}
//		maspc.add(expsr1 >= Dt[t]);
		expsr1.end();
		}
IloExpr submasobjc(lazenv);
   for(IloInt j=0; j<J; j++)
		for(IloInt t=1; t<=T; t++){
			submasobjc+=cj[j]*gc[j][t];
		}
   for(IloInt t=1; t<=T; t++)
	   for(IloInt n= 0; n < N; n++){
		   submasobjc+=ldpen*dc[itrc][n][t];
		   submasobjc-=sellprice*sc[itrc][n][t];
		   }
maspc.add(eta>=submasobjc);
submasobjc.end();

IloCplex ccg_mas(lazenv);
ccg_mas.extract(maspc);
//	ccg_mas.exportModel("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/uc_ccg_mas.lp") ;
	ccg_mas.setParam(IloCplex::MIPSearch,solver_opt);  //1, traditional branch & cut
	if(presol==0){
	ccg_mas.setParam(IloCplex::PreInd,false);
		}
	ccg_mas.setParam(IloCplex::EpGap, sol_gap);
	if (!ccg_mas.solve()) {
		lazenv.error() << "Failed to solve ccg mas model" << endl;
		throw(-1);
		}
	cout<<"ccg mas obj: "<<ccg_mas.getObjValue()<<endl;
	cout<<"obj: "<<costuvw+ccg_mas.getObjValue()<<endl;
//	cout<<" diff: "<<obj_val-costuvw - ccg_mas.getObjValue()<<endl;
//  update
	for(IloInt t = 0; t <= T; t++)
		for(IloInt j= 0; j < J; j++){
			gcstar[j][t]=ccg_mas.getValue(gc[j][t]);
//		ccgresults<<""<<gcstar[j][t]<<endl;
			}
	lb_sub = ccg_mas.getObjValue(); // best feasible for LB
	ingap = abs((ub_sub -lb_sub)/ub_sub) ;
	end=clock();
	runningtime=double(end-begin)/CLOCKS_PER_SEC;
	ccgresults<<"iteration: "<<itrc<<" ccg_lb: "<<lb_sub<<" ccg_ub: "<<ub_sub<<" gap: "<<ingap<<" time:"<<runningtime<<endl;
	ccg_mas.end();
	if(ingap<in_tlrn) break;
//***************************CCG sub problem*******************************
	BoolVarMatrix zc(lazenv, N);	 //pi5lt
	  for(IloInt n=0;n<N;n++){zc[n]=IloBoolVarArray(lazenv,T+1);}
	NumVarMatrix pi5c(lazenv, L);	 //pi5lt
	  for(IloInt l=0;l<L;l++){pi5c[l]=IloNumVarArray(lazenv,T+1,-IloInfinity,IloInfinity);}
	NumVarMatrix pi6c(lazenv, N);	 //p6nt
	  for(IloInt n=0;n<N;n++){pi6c[n]=IloNumVarArray(lazenv,T+1,0,IloInfinity);}
	NumVarMatrix pi6ac(lazenv, N);	 //p6nt
	  for(IloInt n=0;n<N;n++){pi6ac[n]=IloNumVarArray(lazenv,T+1,0,IloInfinity);}
	NumVarMatrix pi7c(lazenv, L);	 //pi7lt
	  for(IloInt l=0;l<L;l++){pi7c[l]=IloNumVarArray(lazenv,T+1,0,IloInfinity);}
	NumVarMatrix pi8c(lazenv, L);	 //pi8lt
	  for(IloInt l=0;l<L;l++){pi8c[l]=IloNumVarArray(lazenv,T+1,0,IloInfinity);}
	NumVarMatrix pi9c(lazenv, N);	 //p9nt
	  for(IloInt n=0;n<N;n++){pi9c[n]=IloNumVarArray(lazenv,T+1,0,IloInfinity);}
	NumVarMatrix pi10c(lazenv, N);	 //p10nt
	  for(IloInt n=0;n<N;n++){pi10c[n]=IloNumVarArray(lazenv,T+1,0,IloInfinity);}
	NumVarMatrix pi11c(lazenv, N);	 //p11nt
	  for(IloInt n=0;n<N;n++){pi11c[n]=IloNumVarArray(lazenv,T+1,0,IloInfinity);}
	NumVarMatrix pi12c(lazenv, N);	 //p11nt
	  for(IloInt n=0;n<N;n++){pi12c[n]=IloNumVarArray(lazenv,T+1,0,IloInfinity);}
	NumVarMatrix beta1c(lazenv, N);	 //linearization
	  for(IloInt n=0;n<N;n++){beta1c[n]=IloNumVarArray(lazenv,T+1,0,IloInfinity);}
	NumVarMatrix beta2c(lazenv, N);	 //linearization
	  for(IloInt n=0;n<N;n++){beta2c[n]=IloNumVarArray(lazenv,T+1,0,IloInfinity);}
	NumVarMatrix beta3c(lazenv, N);	 //linearization
	  for(IloInt n=0;n<N;n++){beta3c[n]=IloNumVarArray(lazenv,T+1,0,IloInfinity);}

//* separate sub problem over time
IloNumArray	  ccgsubobjv(lazenv,T+1);
{for(IloInt t = 1; t <=T; t++)
	{
	IloModel subpc(lazenv);
	Md=ldpen;
	{for(IloInt n = 0; n < N; n++)	//dnt
	for(IloInt it=1; it<=itrc; it++){
			subpc.add(beta1c[n][t] <= pi6c[n][t]);
//			subpc.add(beta2c[n][t] <= pi6ac[n][t]);
			subpc.add(beta3c[n][t] <= pi11c[n][t]);
			subpc.add(beta1c[n][t] <= Md*zc[n][t]);
//			subpc.add(beta2c[n][t] <= Md*zc[n][t]);
			subpc.add(beta3c[n][t] <= Md*zc[n][t]);
			subpc.add(beta1c[n][t] >= Md*zc[n][t]+pi6c[n][t]-Md);
//			subpc.add(beta2c[n][t] >= Md*zc[n][t]+pi6ac[n][t]-Md);
			subpc.add(beta3c[n][t] >= Md*zc[n][t]+pi11c[n][t]-Md);
			}}
IloExpr subobjc(lazenv); //&&&& obj function
	{for(IloInt l = 0; l < L; l++){
		subobjc += (-Pl[l]*(pi7c[l][t]+pi8c[l][t]));
			}}
	{for(IloInt n = 0; n < N; n++){
		subobjc += -maxdelta*(pi9c[n][t]+pi10c[n][t]);
//		subobjc +=  dnt[n][t]*(pi6c[n][t]-pi6ac[n][t]-pi11c[n][t]);
//		subobjc +=  dpnt[n][t]*(beta1c[n][t]-beta2c[n][t]-beta3c[n][t]);
		subobjc +=  dnt[n][t]*(pi6c[n][t]-pi11c[n][t]);
		subobjc +=  dpnt[n][t]*(beta1c[n][t]-beta3c[n][t]);
			}}
IloNum locgc=0;
	{for(IloInt n = 0; n < N; n++){
			locgc = 0 ;
			for(IloInt i= 0; i < J; i++){
				if(Gloc[i]==n)
					locgc+=gcstar[i][t];}
//			subobjc += -locgc*(pi6c[n][t]-pi6ac[n][t]+pi12c[n][t]);
			subobjc += -locgc*(pi6c[n][t]+pi12c[n][t]);
			}}
	IloNum totalgstar = 0;
	{for(IloInt j= 0; j < J; j++){
			totalgstar += cj[j]*gcstar[j][t];
			}}
	subpc.add(IloMaximize(lazenv,subobjc+totalgstar));
	subobjc.end();
//--------- dual constraints
		IloInt lm;
		IloInt lr;
	{for(IloInt l = 0; l < L; l++) // plt
		{
			lm=line[l][0];
			lr=line[l][1];
//			subpc.add((x[l]/100)*pi5c[l][t]-pi6c[lm][t]+pi6ac[lm][t]
//			+pi6c[lr][t]-pi6ac[lr][t]+pi7c[l][t]-pi8c[l][t] == 0);
			subpc.add((x[l]/100)*pi5c[l][t]-pi6c[lm][t]+pi6c[lr][t]+pi7c[l][t]-pi8c[l][t] == 0);
			}}
	{for(IloInt n = 0; n < N; n++)	// delta nt
		{
		IloExpr isexpr1(lazenv);
		IloExpr isexpr2(lazenv);
		 for(IloInt l = 0; l < L; l++){
			if(line[l][1]==n)  // D(l)
			{
					isexpr1+= pi5c[l][t];
			}//end if
			if(line[l][0]==n)  // O(l)
			{
					isexpr2+= -pi5c[l][t];
			} //end if
		}//end for l
		subpc.add(isexpr2+isexpr1+pi9c[n][t]-pi10c[n][t]==0);
		isexpr1.end();
		isexpr2.end();
		}}
	{for(IloInt n = 0; n < N; n++)	//dnt
		{ //subpc.add(pi6c[n][t]-pi6ac[n][t]-pi11c[n][t]<=ldpen);
			subpc.add(pi6c[n][t]-pi11c[n][t]<=ldpen);}}

	{for(IloInt n = 0; n < N; n++)	//snt
		{//subpc.add(-pi6c[n][t]+pi6ac[n][t]-pi12c[n][t]<=-sellprice);
		 subpc.add(-pi6c[n][t]-pi12c[n][t]<=-sellprice);}}

	IloExpr bugc(lazenv);
	   {for(IloInt n = 0; n < N; n++){
			bugc += zc[n][t];}
		subpc.add(bugc<=uncdbug);} // budget on demand uncertainty
		bugc.end();

// update with fixing info
for(IloInt n = 0; n < N; n++)
	{	if(Z[n][t]==0){
		subpc.add(zc[n][t]==0);
		cout<<"Znt"<<n<<t<<"=0"<<endl;}
		if(Z[n][t]==1){
			subpc.add(zc[n][t]==1);
			cout<<"Znt"<<n<<t<<"=1"<<endl;
			}
		}

	IloCplex ccg_subpc(subpc);
	ccg_subpc.exportModel("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/uc_ccgsub.lp") ;
	ccg_subpc.setParam(IloCplex::MIPSearch,solver_opt);  //1, traditional branch & cut
	ccg_subpc.setParam(IloCplex::EpGap, sol_gap);
		if(presol==0){
	ccg_subpc.setParam(IloCplex::PreInd,false);
		}
	if (!ccg_subpc.solve()) {
		lazenv.error() << "Failed to solve sub model" << endl;
		throw(-1);
		}
	//update
		{for(IloInt n= 0; n < N; n++){
			Zntstar[itrc+1][n][t]=ccg_subpc.getValue(zc[n][t]);
		}}
	ccgsubobjv[t] = ccg_subpc.getObjValue();//taking 2% gap into consideration
	cout<<ccg_subpc.getObjValue()<<" "<<ccg_subpc.getCutoff()<<endl;
	ccg_subpc.end();
	subpc.end();
}} // END t, decomposed

double ub_sub_current = 0;
{for(IloInt t = 1; t < T+1; t++){
		ub_sub_current += ccgsubobjv[t];
	}}
ub_sub = IloMin(ub_sub,ub_sub_current); // best feasible for UB
	ingap = (ub_sub -lb_sub )/ub_sub ;
	cout<<"ccg sub obj: "<<ub_sub<<endl;
	cout<<"ccg LB: "<<lb_sub<<endl;
//	cout<<"obj: "<<stcost+ub_sub<<endl;
	cout<<"ingap: "<<ingap<<endl;
	end=clock();
	runningtime=double(end-begin)/CLOCKS_PER_SEC;
	ccgresults<<"iteration: "<<itrc<<" ccg_lb: "<<lb_sub<<" ccg_ub: "<<ub_sub<<" gap: "<<ingap<<" time: "<<runningtime<<endl;
 }//end while
 IloNum objtgc = 0;
 {for(IloInt t = 1; t < T+1; t++)
  for(IloInt j = 0; j < J; j++){
	  objtgc+=cj[j]*gcstar[j][t];
	  }}
  ccgresults<<"total gencost: "<<objtgc<<endl;
  ccgresults<<"load shed: "<<(lb_sub-objtgc)/ldpen<<endl;
  maspc.end();
  end=clock();
  runningtime=double(end-begin)/CLOCKS_PER_SEC;
  ccgresults<<"running time is "<<runningtime<<endl;

//************************* estimate Lower Bound *************
IloNum estlb = IloInfinity;
IloNum temUB = IloInfinity;
while(itrc>0)
	{
	itrc--;
	IloModel lbest(lazenv);
//** ccg master cconstraints**
for(IloInt j=0; j<J; j++){
	lbest.add(gc[j][0]==0);
	}
for(IloInt j=0; j<J; j++)
	for(IloInt t=1; t<=T; t++){
		lbest.add(gc[j][t-1]-gc[j][t] <= RD[j]*laz_u[j][t]+SD[j]*laz_w[j][t]); //Ramping down----traditional
		}
for(IloInt j=0; j<J; j++)
	for(IloInt t=1; t<=T; t++){
		lbest.add(gc[j][t]-gc[j][t-1] <= RU[j]*laz_u[j][t-1]+SU[j]*laz_v[j][t]); //Ramping Up----traditional
		}
for(IloInt j=0; j<J; j++)
	for(IloInt t= 1; t <= T; t++){
			lbest.add(gc[j][t]>=laz_u[j][t]*Gj[j]); //......Generation lower bound limit
			lbest.add(gc[j][t]<=laz_u[j][t]*Gjb[j]);
			}
//......Transmission Network Constraints.
	for(IloInt l = 0; l < L; l++)
		 for(IloInt t = 1; t <= T; t++)  {
				pc[0][l][t].setBounds(-Pl[l],Pl[l]);
				lbest.add(pc[0][l][t]>=-Pl[l]);
				lbest.add(pc[0][l][t]<=Pl[l]);
		 }
//......For L lines, DC power flow.
  {for(IloInt t= 1; t <= T; t++)
	 for(IloInt l=0; l<L;l++)
	 {
		 IloInt m,r;
		 m=line[l][0]; //origin bus
		 r=line[l][1]; //destination bus
		 lbest.add(pc[0][l][t]*x[l]==100*(deltac[0][m][t]-deltac[0][r][t]));
	 }}
  {for(IloInt t= 1; t <= T; t++)
	 for(IloInt n=0; n< N;n++)
	 {
		 lbest.add(dc[0][n][t]<=dnt[n][t]+Zntstar[itrc][n][t]*dpnt[n][t]);
		 lbest.add(deltac[0][n][t]<=maxdelta);
		 lbest.add(deltac[0][n][t]>=-maxdelta);
	 }}
//......load Balancing constraints, for each bus.
   for(IloInt t= 1; t <= T; t++)
	for(IloInt n= 0; n < N; n++)
		{IloExpr exp60c(lazenv);
		 IloExpr exp70c(lazenv);
		 IloExpr exp80c(lazenv);
		{for(IloInt i= 0; i < J; i++)
		{
			if(Gloc[i]==n){
				exp60c+=gc[i][t];
				}
		}}
	   {for(IloInt l = 0; l < L; l++)
		{
				if(line[l][1]==n)  // D(l)
					exp70c+=pc[0][l][t];
				else if(line[l][0]==n) // O(l)
					exp80c+=pc[0][l][t];
	   }}
		lbest.add(exp60c+exp70c-exp80c+dc[0][n][t]-sc[0][n][t] >= dnt[n][t]+Zntstar[itrc+2][n][t]*dpnt[n][t]); //load shed allowed
		lbest.add(sc[0][n][t] <= exp60c);
		lbest.add(dc[0][n][t] <= dnt[n][t]+Zntstar[itrc+2][n][t]*dpnt[n][t]);
		exp60c.end();
		exp70c.end();
		exp80c.end();
	}//end n

// generation greater than demand
IloExpr submasobjc(lazenv);
   for(IloInt j=0; j<J; j++)
		for(IloInt t=1; t<=T; t++){
			submasobjc+=cj[j]*gc[j][t];
		}
   for(IloInt t=1; t<=T; t++)
	   for(IloInt n= 0; n < N; n++){
		   submasobjc+=ldpen*dc[0][n][t];
		   submasobjc-=sellprice*sc[0][n][t];
		   }
lbest.add(IloMinimize(lazenv,submasobjc));
submasobjc.end();

IloCplex cplexlb(lbest);
//	cplexlb.exportModel("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/uc_ccg_mas.lp") ;
	cplexlb.setParam(IloCplex::MIPSearch,solver_opt);  //1, traditional branch & cut
		if(presol==0){
	cplexlb.setParam(IloCplex::PreInd,false);
		}
	cplexlb.setParam(IloCplex::EpGap, sol_gap);
	if (!cplexlb.solve()) {
		lazenv.error() << "Failed to solve lb estimation model" << endl;
		throw(-1);
		}
	cout<<"subp LB estimation: "<<cplexlb.getObjValue()<<endl;
	ccgresults<<"subp LB estimation: "<<cplexlb.getObjValue()<<endl;
	cout<<"obj: "<<costuvw+cplexlb.getObjValue()<<endl;
	ccgresults<<"obj: "<<costuvw+cplexlb.getObjValue()<<endl;
//update, pick the zstar with lowest LB
if(cplexlb.getObjValue()<estlb){
//	temUB = costuvw+cplexlb.getObjValue();
	for(IloInt t=1; t<=T; t++)
		for(IloInt n= 0; n < N; n++){
			Zoutstar[itrc_out][n][t]=Zntstar[itrc+2][n][t];
			Znode[n][t]=Zntstar[itrc+2][n][t]; // current best feasible solution
		}
	}
cplexlb.end();
lbest.end();
} //end innermost while
outUB = IloMin(costuvw+lb_sub,outUB);
outgap = (outUB-outLB)/outUB;
ccgresults<<"---outUB: "<<outUB<<" outGap: "<<outgap <<endl;
cout<<"good til"<<endl;
}// end try
	catch (IloException& ex){	  cerr << "Error: " << ex << endl;	}
	catch (...)				{	  cerr << "Error" << endl;			}
	lazenv.end();
//update bounds
UB = ub_sub;
LB = IloMax(lb_sub,LB);
cout<<"Node info: "<<"UB--"<<UB<<" LB--"<<LB<<endl;
return 0;
}//end ubcal