/*
Given a protection plan, 
*/

#ifndef DEFINEDMATRIX_H
#define DEFINEDMATRIX_H
#include "definedmatrix.h"
#endif DEFINEDMATRIX_H
#ifndef PARINIT_H
#define PARINIT_H
#include "parinit.h"
#endif PARINIT_H

//#include "defined_matrix.h"
using namespace std;

int callb()
	{

{IloInt m; // O(l),origin or line
IloInt r; // D(l),destination of line
IloNum scale=1; // nodal demand scale
IloNum maxdelta= 1.57;
		double obj_val;
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
//sellprice = 0;
//it_out = 20; // # of out iterations
//		NumMatrix dnt(env);
		IloNum costuvw; //cost of start up and fixed cost
		IloNum Md; // penalty for dual
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
		cout<<"parameter initialization done"<<endl;
}

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
	cplexlb.exportModel("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/uc_ccg_mas.lp") ;
	cplexlb.setParam(IloCplex::MIPSearch,0);  //1, traditional branch & cut
	cplexlb.setParam(IloCplex::EpGap, 0.001);
	if (!cplexlb.solve()) {
		lazenv.error() << "Failed to solve model" << endl;
		throw(-1);
		}
	cout<<"subp LB estimation: "<<cplexlb.getObjValue()<<endl;
	ccgresults<<"subp LB estimation: "<<cplexlb.getObjValue()<<endl;
	cout<<"obj: "<<stcost+cplexlb.getObjValue()<<endl;
	ccgresults<<"obj: "<<stcost+cplexlb.getObjValue()<<endl;
//update, pick the zstar with lowest LB
if(cplexlb.getObjValue()<estlb){
//	temUB = costuvw+cplexlb.getObjValue();
	for(IloInt t=1; t<=T; t++)
		for(IloInt n= 0; n < N; n++){  
			Zoutstar[itrc_out][n][t]=Zntstar[itrc+2][n][t];
		}
	}
cplexlb.end();
lbest.end();
} // end innermost while 
outUB = IloMin(costuvw+lb_sub,outUB);
outgap = (outUB-outLB)/outUB;
ccgresults<<"---outUB: "<<outUB<<" outGap: "<<outgap <<endl;
}

// scan and add lazy cuts
// updating for each run
for(IloInt t = 1; t <= T; t++){
	Dt[t]=0;
	for(IloInt n = 0; n < N; n++){
				Dt[t]+=dnt[n][t]+Zoutstar[itrc_out][n][t]*dpnt[n][t];
			    Md+=dnt[n][t]+Zoutstar[itrc_out][n][t]*dpnt[n][t];
				}
			}
Md=IloMax(Dt);
// upper bound  
	// find kjt: 1,..,UT 
	for(IloInt j=0; j<J; j++)
	    for(IloInt t= 1; t <= T; t++){
		    int k=1;
	        while( SD[j]+(k-1)*RD[j] < Gjb[j] && k+t <T && k <= UT[j] )
				k++;		
			Kjt[j][t]=k-1;// violates conditions, decrease  1
		} cout<<"this up added"<<endl;
	 for(IloInt j=0; j<J; j++)
		for(IloInt t=1; t<=T; t++){	 
		IloExpr Expub(env1);
			for(IloInt i=1; i <= Kjt[j][t];i++){
				Expub+=(SD[j]+(i-1)*RD[j])*w[j][t+i]-Gjb[j]*v[j][t+i];
				}
			add(g[itrc_out][j][t]<= Gjb[j]*u[j][t+Kjt[j][t]]+Expub).end();
		Expub.end();
		}   

	}