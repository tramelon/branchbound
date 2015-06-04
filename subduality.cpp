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

int subduality(double &UB, double &LB, int U[][T1+1], int V[][T1+1], int W[][T1+1],float Znode[][T1+1], int budget, double costuvw)
	{

	IloInt m; // O(l),origin or line
	IloInt r; // D(l),destination of line
	IloNum scale=1; // nodal demand scale
	IloNum maxdelta= 1.57;
	double obj_val;
	IloNum subobjval;
	double ldpen=ldpena; // penalty for load shed
	IloInt uncdbug = budget;  // demand uncertain budget
	double uncd = uncda; // uncertainty part / demand 
	double in_tlrn = 0.005 ; // in gap torrence, 0.01
	double sol_gap = 0.001 ; // solver gap
	IloInt itrc_out = 0;
	IloInt res_itrc = 20;// initial set of scenarios in mas problem
	IloNum outgap = 1;
	IloNum outLB = -IloInfinity;
	IloNum outUB =  IloInfinity;
	double runningtime;
//	IloInt J=J;	//units
//	IloInt T=T1;	//time
//	IloInt N=Nn;    // buses
//	IloInt L=186;    // trans lines
	IloEnv env;
	
try{
// NumMatrix dnt(env);
		IloNum Md; // penalty for dual
		IloNum sellprice = 0;
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
		for(IloInt j=0; j<J; j++){
			gj0[j]=0;
			}
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

BoolVarMatrix znt(env,N);
	  for(IloInt n=0;n<N;n++){znt[n]=IloBoolVarArray(env,T+1);}
	NumVarMatrix beta1(env, N);	 //linearization
	  for(IloInt n=0;n<N;n++){beta1[n]=IloNumVarArray(env,T+1,0,IloInfinity);}
	NumVarMatrix beta2(env, N);	 //linearization
	  for(IloInt n=0;n<N;n++){beta2[n]=IloNumVarArray(env,T+1,0,IloInfinity);}	
	NumVarMatrix beta3(env, N);	 //linearization
	  for(IloInt n=0;n<N;n++){beta3[n]=IloNumVarArray(env,T+1,0,IloInfinity);}
	NumVarMatrix q(env, N);	 //uncertain demand
	  for(IloInt n=0;n<N;n++){q[n]=IloNumVarArray(env,T+1,0,IloInfinity);}
	BoolVarMatrix z(env, N);	 //uncertain demand
	  for(IloInt n=0;n<N;n++){z[n]=IloBoolVarArray(env,T+1);}
	NumVarMatrix pi1(env, J);	 //p1jt
	  for(IloInt j=0;j<J;j++){pi1[j]=IloNumVarArray(env,T+1,0,IloInfinity);}   	
	NumVarMatrix pi2(env, J);	 //p2jt
	  for(IloInt j=0;j<J;j++){pi2[j]=IloNumVarArray(env,T+1,0,IloInfinity);}
	NumVarMatrix pi3(env, J);	 //p3jt
	  for(IloInt j=0;j<J;j++){pi3[j]=IloNumVarArray(env,T+1,0,IloInfinity);}
	NumVarMatrix pi4(env, J);	 //p4jt
	  for(IloInt j=0;j<J;j++){pi4[j]=IloNumVarArray(env,T+1,0,IloInfinity);}
	NumVarMatrix pi5(env, L);	 //pi5lt
	  for(IloInt l=0;l<L;l++){pi5[l]=IloNumVarArray(env,T+1,-IloInfinity,IloInfinity);}
	NumVarMatrix pi6(env, N);	 //p6nt
	  for(IloInt n=0;n<N;n++){pi6[n]=IloNumVarArray(env,T+1,-IloInfinity,IloInfinity);}	   
	NumVarMatrix pi7(env, L);	 //pi7lt
	  for(IloInt l=0;l<L;l++){pi7[l]=IloNumVarArray(env,T+1,0,IloInfinity);}
	NumVarMatrix pi8(env, L);	 //pi8lt
	  for(IloInt l=0;l<L;l++){pi8[l]=IloNumVarArray(env,T+1,0,IloInfinity);}
	NumVarMatrix pi9(env, N);	 //p9nt
	  for(IloInt n=0;n<N;n++){pi9[n]=IloNumVarArray(env,T+1,0,IloInfinity);}
	NumVarMatrix pi10(env, N);	 //p10nt
	  for(IloInt n=0;n<N;n++){pi10[n]=IloNumVarArray(env,T+1,0,IloInfinity);}
	NumVarMatrix pi11(env, N);	 //p11nt
	  for(IloInt n=0;n<N;n++){pi11[n]=IloNumVarArray(env,T+1,0,IloInfinity);}
	IloNumVarArray  pi12(env,T+1,0,IloInfinity);
	IloNumVarArray  pi14(env,T+1,0,IloInfinity);
	NumVarMatrix mq(env, N);	 // linearization 
	  for(IloInt n=0;n<N;n++){mq[n]=IloNumVarArray(env,T+1,0,IloInfinity);}	
//--------------obj func
	  IloExpr subobj(env);
	for(IloInt l = 0; l < L; l++)
		for(IloInt t = 1; t <= T; t++){
		subobj += (-Pl[l]*(pi7[l][t]+pi8[l][t]));
		}		
	for(IloInt n = 0; n < N; n++)
		for(IloInt t = 1; t <= T; t++){
//		subobj += (-maxdelta*(pi9[n][t]+pi10[n][t])+dnt[n][t]*(pi6[n][t]-pi11[n][t]));
		subobj += (-maxdelta*(pi9[n][t]+pi10[n][t])+dnt[n][t]*(pi6[n][t]-pi11[n][t]))+dpnt[n][t]*(beta1[n][t]-beta2[n][t]);
		}	
	for(IloInt j = 0; j < J; j++)
		for(IloInt t = 1; t <= T; t++){
		subobj += (-(RD[j]*U[j][t]+SD[j]*W[j][t])*pi1[j][t]-(RU[j]*U[j][t-1]+SU[j]*V[j][t])*pi2[j][t]);
		subobj += Gj[j]*U[j][t]*pi3[j][t]-Gjb[j]*U[j][t]*pi4[j][t];
		if(t==1){
			subobj += gj0[j]*(pi1[j][t] -pi2[j][t]);
		}
	}
	for(IloInt t = 1; t <= T; t++){
//		subobj+= Dt[t]*pi12[t];
		}
	for(IloInt n = 0; n < N; n++)
		for(IloInt t = 1; t <= T; t++){
			subobj+= dnt[n][t]*pi14[t]+dpnt[n][t]*beta3[n][t];
			}
	IloModel subp(env);
	subp.add(IloMaximize(env,subobj));
	subobj.end();
//---------constraints---------------------------//
//------- linearization constraints

	for(IloInt t = 1; t <= T; t++){
		IloExpr subbug(env);
		for(IloInt n = 0; n < N; n++){
			subp.add(znt[n][t]>=0);
			subbug += znt[n][t];} 
		subp.add(subbug<=uncdbug); // budget on demand uncertaint
		subbug.end();
		}
	Md=500;
	for(IloInt n = 0; n < N; n++)	//dnt
		for(IloInt t = 1; t <= T; t++){
			subp.add(beta1[n][t] <= pi6[n][t]);
			subp.add(beta2[n][t] <= pi11[n][t]);
			subp.add(beta3[n][t] <= pi14[t]);
			subp.add(beta1[n][t] <= Md*znt[n][t]);
			subp.add(beta2[n][t] <= Md*znt[n][t]);
			subp.add(beta3[n][t] <= Md*znt[n][t]);
			subp.add(beta1[n][t] >= Md*znt[n][t]+pi6[n][t]-Md);
			subp.add(beta2[n][t] >= Md*znt[n][t]+pi11[n][t]-Md);
			subp.add(beta3[n][t] >= Md*znt[n][t]+pi14[t]-Md);
			}	

//--------- dual constraints
	for(IloInt j = 0; j < J; j++)  // gjt
		for(IloInt t = 1; t <= T-1; t++){
			subp.add(pi1[j][t]-pi1[j][t+1]-pi2[j][t]+pi2[j][t+1]+pi3[j][t]-pi4[j][t]+pi6[Gloc[j]][t]/* + pi12[t]*/+ pi14[t]<= cj[j]);
		}
	for(IloInt j = 0; j < J; j++){
			subp.add(pi1[j][T]-pi2[j][T]+pi3[j][T]-pi4[j][T]+pi6[Gloc[j]][T]/*+ pi12[T]*/+ pi14[T]<= cj[j]);
		}
		IloInt lm;
		IloInt lr;
	for(IloInt l = 0; l < L; l++) // plt
		for(IloInt t = 1; t <= T; t++){
			lm=line[l][0];
			lr=line[l][1];
			subp.add((x[l]/100)*pi5[l][t]-pi6[lm][t]+pi6[lr][t]+pi7[l][t]-pi8[l][t] == 0);
		}
	for(IloInt n = 0; n < N; n++)	// delta nt
		for(IloInt t = 1; t <= T; t++){
		IloExpr iexpr1(env);
		IloExpr iexpr2(env);
		 for(IloInt l = 0; l < L; l++)
		{
			if(line[l][1]==n)  // D(l)
			{
					iexpr1+= pi5[l][t];
			}//end if
			if(line[l][0]==n)  // O(l)
			{
					iexpr2+= -pi5[l][t];
			} //end if
		}//end for l
		subp.add(iexpr2+iexpr1+pi9[n][t]-pi10[n][t]==0);
		iexpr1.end();
		iexpr2.end();
		}
	for(IloInt n = 0; n < N; n++)	//dnt
		for(IloInt t = 1; t <= T; t++){
			subp.add(pi6[n][t]-pi11[n][t]<=ldpen);
		}
	IloCplex cplex_sub(subp);
	
	if(presol==0){
	cplex_sub.setParam(IloCplex::PreInd,false);
		}
//	cplex_sub.exportModel("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/uc_sub.lp") ;
	cplex_sub.setParam(IloCplex::EpGap, 0.005);
	cplex_sub.setParam(IloCplex::MIPSearch,0);  //1, traditional branch & cut
	if (!cplex_sub.solve()) {
		env.error() << "Failed to solve model" << endl;
		throw(-1);
		}
		
	subobjval = cplex_sub.getObjValue();
	cout<<"dual obj: "<<cplex_sub.getObjValue()<<endl;
	cout<<"obj: "<<costuvw+cplex_sub.getObjValue()<<endl;
//	cout<<" diff: "<<obj_val-costuvw- cplex_sub.getObjValue()<<endl;
if(costuvw+subobjval <= outUB){
	outUB = IloMin(costuvw+subobjval,outUB);
	outgap = (outUB-outLB)/outUB;
	ccgresults<<"---outUB: "<<outUB<<" outGap: "<<outgap <<endl;
	if(uncdbug>=1){
	for(IloInt t=1; t<=T; t++)
		for(IloInt n= 0;n<N;n++){
			Znode[n][t]=cplex_sub.getValue(z[n][t]); // current best feasible solution
			}}
	}

	cplex_sub.end();
	subp.end();
}// end try
	catch (IloException& ex){	  cerr << "Error: " << ex << endl;	}
	catch (...)				{	  cerr << "Error" << endl;			}
	env.end();
//update bounds
UB = subobjval;
LB = subobjval;
return 0;
}//end ubcal