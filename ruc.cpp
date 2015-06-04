/*
* Tight MILP formulation of Unit Commitment Problem with Transmission Network
* Wei Yuan
* 2014-01-23
* version 1

* Note: t index from: 0 -> T, but only use 1 <= t <= T, input function must modify accordingly
* Network reduction based on: Jame Ostrowski
* unified bound based on min of phase angle and flow limit.
* Network partition based on tight transmission lines,2014-03-11
* call back function test,2014-04-19
* add uncertain demand,2014-4-22
* verify sub problem formulation,2014-4-24
* part of subproblem with fixed demand is done sub_callback_1.cpp to verify the formulation
* sub_callback_2.cpp tries to include the full subproblem with uncertain demand. 2014-05-06

* KKT conditions with dual objective function..........5/23/2014
* Sub problem using ccg
* switched order of max min min to min max min for a Upper Bound 6/12/2014
* estimate lower bound for max min.
* outer mas problem added to form loop 6/27/2014

* ccg for out sub problem is implemented in calub, previous use N->LB(in mas obj), now I use N-> to gurantee no fault
* subproblem solution using branch and bound, 10/02/2014 
*/

/*
dpnt: limit of uncertain part in demand, 25% of normal demand
*/

//#include <ilcplex/ilocplex.h>
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
#ifndef DEFINEDMATRIX_H
#define DEFINEDMATRIX_H
#include "definedmatrix.h"
#endif DEFINEDMATRIX_H
#ifndef CALUB_H
#define CALUB_H
#include "calub.h"
#endif CALUB_H
#ifndef PARINIT_H
#define PARINIT_H
//#include "parinit.h"
#endif PARINIT_H
using namespace std;
/*output location*/
ofstream results("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/resultsB3e.txt");
ofstream Duals("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/dualsB3e.txt");
ofstream flows("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/flowsB3e.txt");
ofstream ccgresults("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/ccgB3e.txt");
ofstream genu("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/genu3e.txt");
ofstream objval("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/objval3e.txt");
//ofstream brief("C:/Dropbox/RobustOpt/rucnk/results/brief_k1.txt");
//ofstream runtimes("C:/Dropbox/RobustOpt/rucnk/results/runtimes_k1.txt");
//ofstream bounds("C:/Dropbox/RobustOpt/rucnk/results/bounds_k1.txt");

ILOSTLBEGIN
//#include "defined_matrix.h"
//global variables
double BB_UB;   // update when a whole level of leaf node is processed
double BB_LB;   // best feasible solution bound, updated every node
double BB_gap0; // whole process gap
double BB_gap;  // branch and bound tree gap, in head node
double node_gap;
double node_gap0; // gap to redeem a node to be optimal
double outgap0 = 0.01; // out loop tolerable gap 
double outmastol = 0.005; // out mas problem tolerable gap
double in_tlrn = 0.01 ; // in gap tolerable, not used in branch and bound routine
double inmastol=0.001;// in mas tolerable, not used in calub(), optional
double insubtol=0.001;// above
double costuvw;
int uncbg = 3;  // budget for demand uncertainty
int warmst = 0; // 1 use warm start, 0, OT
int brahopt = 1; // 1, use branch and bound; 0, OT
int masdif = 0; //1, force master problem to a different solution, line 950


int leaf ;
int BB1;
int Level;  // level of bb tree
#define L 186
#define J 54
#define N1 118
#define T 24
ofstream BBres("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/BBres_rus.txt");
float Zb[N1][T+1]; // best feasible solution so far
//-------------------------------------index
//		IloInt Nbus=118;    // buses
//-------------Decision Variables for Sub problem
	//	NumMatrix dnt(env, N);     // dnt
			
struct node{
	 //used in the tree
	 node *parent, *l, *r; //record the parent and child nodes of current node
	 int level; //depth of current node, 0,1,...
	 int tag;// -1: pruned by bound 
		//0: not evaluated 
		//1:evaluated, but BB_gap>1% 
		//2: optimal, pruned 
		//3: BB_gap>1%, but can not be branched anymore,based on the constraints
	 // used in the chain
	 node *next; 
	 double lb;    
	 double ub; 
	 int Z[N1][T+1]; // var to be fixed
	 float Znode[N1][T+1]; // best feasible solution for current node	
};

node *head;// for chain, store the nodes to be evaluated
node *root;// for tree, store the nodes in the Branch and Bound Tree

void EnQueue(node *i) //enter the chain
{
   i->next=head->next;
   head->next=i;
}

void DeQueue(node *i)// leave the chain
{
   node *pre=head,*p=head->next;
   while(p!=i){
	  pre=p;
	  p=p->next;
   }
   pre->next=p->next;
}

void init()// initialization
{  
  leaf = 0;
  BB1=0;
  Level = 0;
  head=new(node);
  head->next=NULL;
  head->parent = NULL;
  head->l = NULL;
  head->r = NULL;
  head->level = 0; // head as level 0
  head->tag = 1; //evaluated
  head->ub = HUGE_VAL;
  head->lb = -HUGE_VAL;
  for(int n=0; n<N1; n++)// no var fixed
	  for(int t=0; t<=T;t++)
	  head->Z[n][t] = -1;
  BB_UB = HUGE_VAL;
  BB_LB = -HUGE_VAL;   
  BB_gap = 1;
  BB_gap0 = 0.001; // tolerable gap
  node_gap0=0.1*BB_gap0; // node tolerable gap
}

void NewNode(node *prnt) //generate new child nodes from prnt
{ 
	node *il=new(node);
	il->parent = prnt;
	node *ir=new(node);
	ir->parent = prnt;
// expand the tree correspondingly
	prnt->l = il;
	prnt->r = ir;
	il->level = prnt->level + 1;
	ir->level = prnt->level + 1;
	il->l = NULL;
	il->r = NULL;
	ir->l = NULL;
	ir->r = NULL;
	il->tag = 0;
	ir->tag = 0;
////determine the branching variable  z[a][b], first non -1 to branch;
	int a=0;
	int b=1;
	int l=0;
	int fix_index=1;
	int t=1;
	int bug = 0;
	while(t<=T && (fix_index==1))
		// find time(hour) with non-fixed z, satisfy budget constraints
		{	
		fix_index=0;
		bug = 0;
		for(int n=0; n<N1;n++){
			if(prnt->Z[n][t]==1||prnt->Z[n][t]==0)
				fix_index=1;
			if(prnt->Z[n][t]==1)
				bug++;
			}
		if(bug==uncbg) // meet the budget
			fix_index = 1;
		t++;
		}
	b=t-1;
	for(int n=0; n<N1; n++){
		  if(prnt->Z[n][b]==-1 && prnt->Znode[n][t]!=Zb[n][t]){ // first non-fixed bus, zb: best feas solution
			  a=n;
			  }
		  }
	BBres<<"Branch Var:"<<a<<" , "<<b<<endl;
//.......
//----------- update fix var info
	ir->Z[a][b]=0; //right child
	for(int n=0; n<N1; n++)
	  for(int t=0; t<=T;t++){
		if(n!=a||t!=b)
			ir->Z[n][t]=prnt->Z[n][t];
		  }
	ir->ub=prnt->ub;
	ir->lb=prnt->lb;
	if(ir->lb <= BB_UB)
		EnQueue(ir);
	il->Z[a][b]=1; //left child
	for(int n=0; n<N1; n++)
	  for(int t=0; t<=T;t++){
		if(n!=a||t!=b)
			il->Z[n][t]=prnt->Z[n][t];
		  }
	il->ub = prnt->ub;
	il->lb = prnt->lb;
//enter the chain
	if(il->lb <= BB_UB)
		EnQueue(il);
}

node *NextLiveNode() //return next live node in the chain, not used in program
{
		node *p=head->next,*choice=p;
	//  int lb=p->lb;
	//	while(p){
	//	   if(p->lb>lb){
	//		  choice=p;
	//	    }
	//	   p=p->next;
	//	}
	//	return(choice);
}



void update0(node *N) // see update(), update bounds when bounding info computed, for a single node and its child nodes
{
	if(N->l != NULL && N->r != NULL) // both child
	{
		if(N->l->ub > N->r->ub) // ub = max (ub_l,ub_r)
			N->ub = N->l->ub;
		else 
			N->ub = N->r->ub;	
		if(N->l->lb > N->r->lb) // lb = max (lb_l,lb_r) 
			N->lb = N->l->lb;
		else 
			N->lb = N->r->lb;
	}
	else if(N->l != NULL && N->r == NULL) // left child only
	{	//N->ub = N->ub; // upper bound stay the same 
		N->lb = N->l->lb; // lower bound can be updated anytime
	}
	else if(N->l == NULL && N->r != NULL)// right child only
		{//N->ub = N->ub; // upper bound stay the same 
		 N->lb = N->r->lb; // lower bound can be updated anytime
		}
	
	if((N->ub - N->lb)/N->ub < node_gap0)//optimal node, pruned, based on the updated bounds
	{
		N->tag = 2;
	}
}// end update0

void findo(node *N, int l)//see update(), recursively call findo, level l
{	
node *N1a, *N2;
if(N->l != NULL && N->r != NULL){ //with two child nodes
	if(N->l->tag == 1 && N->r->tag == 1) //1:both child evaluated, but BB_gap>1%
	{
		if(N->l->level == l) // leaf node, update child node until level l
		{
			N1a=N->l;
			N2=N->r;
			update0(N1a);
			update0(N2);
		}
		else  // no update required
		{
			findo(N->l,l); //nested, process child nodes
			findo(N->r,l);
		}
	}
	if(N->l->tag != 1 && N->r->tag == 1) // right child evaluated
	{
		if(N->r->level == l) // level l node, update0
		{
			N2=N->r;
			update0(N2);
		}
		else  // non-root node, findo
		{
			findo(N->r,l);
		}
	}
	
	if(N->l->tag == 1 && N->r->tag != 1) // left child evaluated
	{
		if(N->l->level == l){
			N1a=N->l;
			update0(N1a);
		}
		else{
			findo(N->l,l);
		}
	}}

if(N->l != NULL && N->r == NULL) //left child only
	{
		if(N->l->tag == 1 )
		{
			if(N->l->level == l)
			{
				N1a=N->l;
				update0(N1a); // root node
			}
			else
			{
				findo(N->l,l); //non-root node
			}

		}
	}

if(N->r != NULL && N->l == NULL) // right child only
	{
		if(N->r->tag == 1)
		{
			if(N->r->level == l)
			{
				N2=N->r;
				update0(N2);
			}
			else
			{
				findo(N->r,l);
			}
		}
	}
}// end findo

void prune(node *N) // cut off the nodes that do not need further evaluation, prune whole branch 
{
	if(N->lb > BB_UB) // node lb > ub, pruned
		N->tag = -1;
	if(N->l != NULL && N->r != NULL) // source down to child nodes
	{
		prune(N->l);
		prune(N->r);
	}
	if(N->l == NULL && N->r != NULL)
	{
		prune(N->r);
	}
	if(N->l != NULL && N->r == NULL)
	{	
		prune(N->l);
		}
}

float update() // update()+findo()+update0(): 
//update the lb and ub after all nodes in a level are evaluated; update() recursively calls findo() and update0()
{
	float g;
	for(int l=L-1; l>0; l--) // process every level from L-1 to 1
	{
		findo(head, l);
	}
	update0(head);
	if(BB_UB > head->ub) // update overall bounds based on the head node.
		BB_UB = head->ub;
	BB_LB=head->lb;
	g = abs((BB_UB - head->lb)/BB_UB); //gap
	BBres<<"BB_UB="<<BB_UB<<" BB_LB="<<BB_LB<<" BB_gap="<<g*100<<"%"<< " updating head node after whole tree is updated: " <<endl;
	if(g<BB_gap0)
		return g;
	else
	{
//		prune(head); // prune node if gap is not tolerable
		return g;
	}
}


int
main(int argc, char **argv)
{
//	IloEnv env;
	try 
	{	 //------------------------------------------------------------------------parameters
		IloNum Md; // penalty for dual
		IloNumArray CarD(env,N1); //cardinality for uncertain demand
		IloNumArray duals(env,L);
		IloNumArray Fbounds(env,L);
		IloInt m; // O(l),origin or line
		IloInt r; // D(l),destination of line
		IloNum maxdelta= 1.57;
		double obj_val;
		IloNum scale=1; // nodal demand scale
		IloNum ldpen=50; // penalty for load shed
		IloInt uncdbug = uncbg;  // demand uncertain budget
		IloNum uncd =0.1; // uncertainty part / demand
	
		double runningtime;
		IloNum sellprice = 0;
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
		IloNumArray UT0(env,J);    // on line time before schedule
		IloNumArray UTr(env,J);    // 
		IloNumArray DT(env,J);     // min down
		IloNumArray DT0(env,J);    // off line time before schedule
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
	//	NumVarMatrix  ucd(env,N1);
	//	BoolVarMatrix card(env,N1);
	//	IntMatrix line(env,L);      // trans lines
		IloInt lazcount=0;

//-------------Decision Variables for Sub problem
	//	NumMatrix dnt(env, N1);     // dnt
//		NumMatrix dnt(env);
//		NumMatrix Kjt(env, J);     // kjt
		NumMatrix dnt(env, N1);     // dnt
			for(IloInt n=0;n<N1;n++){
				dnt[n]=IloNumArray(env,T+1);
			}	
		NumMatrix dpnt(env, N1);     // dpnt, upper bound of unc part of demand.
			for(IloInt n=0;n<N1;n++){
				dpnt[n]=IloNumArray(env,T+1);
			}
		NumMatrix Kjt(env, J);     // kjt
			for(IloInt j=0;j<J;j++){
				Kjt[j]=IloNumArray(env,T+1);
			}
		NumMatrix Ustar(env, J);     // Ustar,output
			for(IloInt j=0;j<J;j++){
				Ustar[j]=IloNumArray(env,T+1);
			}		
		IntMatrix line(env,L);      // trans lines
		   for(IloInt l = 0; l < L; l++){
				line[l] = IloIntArray(env,2);
			}
		IloInt Betn = 6;
		IntMatrix cut(env, Betn);      // cut trans lines
		   for(IloInt l = 0; l < Betn; l++){
				cut[l] = IloIntArray(env,2);
			} 
		IloBoolArray label(env, N1);
		IloInt Betn1 = 4;
		IntMatrix cut1(env, Betn1);      // cut trans lines
		   for(IloInt l = 0; l < Betn1; l++){
				cut1[l] = IloIntArray(env,2);
			} 
		IloBoolArray label1(env, N1);
		IloInt Betn2 = 12;
		IntMatrix cut2(env, Betn2);      // cut trans lines
		   for(IloInt l = 0; l < Betn2; l++){
				cut2[l] = IloIntArray(env,2);
			} 
		IloIntArray label2(env, N1);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/Sj.txt", Sj, J);
		NumMatrix Csu(env, J);     // Csu_gs
			for(IloInt j=0;j<J;j++){
				Csu[j]=IloNumArray(env,Sj[j]);
			}
		NumMatrix Tsu(env, J);     // Tsu_gs
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
		readMaData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/dnt.txt", dnt, N1, T+1);	
		readIntMaData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/transmission.txt", line, L, 2);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/flowlimit.txt",Pl, L);
		readData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/reactance.txt", x, L);
		Md=0;
		for(IloInt t = 1; t <= T; t++){
			Dt[t]=0;
			for(IloInt n = 0; n < N1; n++){
				Dt[t]+=dnt[n][t];
				Md+=dnt[n][t];
				}
			}
		Md=IloMax(Dt);
		for(IloInt t = 0; t <= T; t++)
			for(IloInt n = 0; n < N1; n++)
				dpnt[n][t]= uncd*dnt[n][t];
			

//------------------------------------MAIN PROGRAM
IloInt it_out; // # of out iterations
it_out = 20;
ThreeDNumMatrix Zoutstar(env,it_out);
		 { for(IloInt it = 0; it < it_out; it++)
			   Zoutstar[it] = NumMatrix(env,N1);
		  for(IloInt it = 0; it < it_out; it++)
				for(IloInt n = 0; n < N1; n++)
					Zoutstar[it][n] = IloNumArray(env,T+1); }
for(IloInt it = 0; it < 2; it++)
		for(IloInt n = 0; n < N1; n++)
			for(IloInt t = 0; t <= T; t++) 
				Zoutstar[it][n][t] = 0;
//-------------------------------------------------------Master Problem Decision Variables		
	BoolVarMatrix u(env,J);  // on/off
		for(IloInt j = 0; j < J; j++)
			u[j] = IloBoolVarArray(env,T+1);
	BoolVarMatrix v(env,J);  // start up
		for(IloInt j = 0; j < J; j++)
			v[j] = IloBoolVarArray(env,T+1);
	BoolVarMatrix w(env,J);  // shut down
		for(IloInt j = 0; j < J; j++)
			w[j] = IloBoolVarArray(env,T+1);
	ThreeDBoolVarMatrix h(env,J); // start-up type
		for(IloInt j = 0; j < J; j++)
			h[j]=BoolVarMatrix(env,Sj[j]);
			for(IloInt j = 0; j < J; j++)
				for(IloInt s = 0; s < Sj[j]; s++) 
					h[j][s]=IloBoolVarArray(env,T+1);
ThreeDVarMatrix g(env,it_out);
for(IloInt it=0;it<it_out;it++)
	g[it]=NumVarMatrix(env,J); // generation level, gjt
	for(IloInt it=0;it<it_out;it++)
		for(IloInt j = 0; j < J; j++)
			g[it][j] = IloNumVarArray(env,T+1, 0, IloInfinity); 
ThreeDVarMatrix goverline(env,it_out);
for(IloInt it=0;it<it_out;it++)
	goverline[it]=NumVarMatrix(env,J); // generation level, gjt
	for(IloInt it=0;it<it_out;it++)
		for(IloInt j = 0; j < J; j++)
			goverline[it][j] = IloNumVarArray(env,T+1, 0, IloInfinity); 
ThreeDVarMatrix reserve(env,it_out);
for(IloInt it=0;it<it_out;it++)
	reserve[it]=NumVarMatrix(env,J); // generation level, gjt
	for(IloInt it=0;it<it_out;it++)
		for(IloInt j = 0; j < J; j++)
			reserve[it][j] = IloNumVarArray(env,T+1, 0, IloInfinity); 
ThreeDVarMatrix plt(env,it_out);
for(IloInt it=0;it<it_out;it++)
	plt[it]=NumVarMatrix(env,L); // power flow, plt
	for(IloInt it=0;it<it_out;it++)
		for(IloInt l = 0; l < L; l++)
			plt[it][l] = IloNumVarArray(env,T+1, -IloInfinity, IloInfinity); 
ThreeDVarMatrix dlnt(env,it_out);
for(IloInt it=0;it<it_out;it++)
	dlnt[it]=NumVarMatrix(env,N1); // load shed level, gjt
	for(IloInt it=0;it<it_out;it++)
		for(IloInt n = 0; n < N1; n++)
			dlnt[it][n] = IloNumVarArray(env,T+1, 0, IloInfinity); 
for(IloInt it=0;it<it_out;it++)
	for(IloInt n = 0; n < N1; n++)
			  for(IloInt t = 1; t <= T; t++) 
				dlnt[it][n][t].setBounds(0,dnt[n][t]);
ThreeDVarMatrix delta(env,it_out);
for(IloInt it=0;it<it_out;it++)
	delta[it]=NumVarMatrix(env,N1); // phase angle, gjt
	for(IloInt it=0;it<it_out;it++)
		for(IloInt n = 0; n < N1; n++)
			delta[it][n] = IloNumVarArray(env,T+1, -maxdelta, maxdelta); 	
/* // for add this iteration variables
		for(IloInt j = 0; j < J; j++)
			g[j] = IloNumVarArray(env,T+1, 0, IloInfinity); 
	NumVarMatrix goverline(env,J);  // generation level, gjt
		for(IloInt j = 0; j < J; j++)
			goverline[j] = IloNumVarArray(env,T+1, 0, IloInfinity); 
	NumVarMatrix reserve(env,J);  // spinning reserve
		for(IloInt j = 0; j < J; j++)
			reserve[j] = IloNumVarArray(env,T+1, 0, IloInfinity); 
	NumVarMatrix plt(env,L);  // power flow, plt
		for(IloInt l = 0; l < L; l++)
			plt[l] = IloNumVarArray(env,T+1, -IloInfinity, IloInfinity);
	NumVarMatrix dlnt(env,N1);  // load shed, dlnt
		for(IloInt n = 0; n < N; n++)
				dlnt[n] = IloNumVarArray(env,T+1, 0, IloInfinity);
	for(IloInt n = 0; n < N; n++)
			  for(IloInt t = 1; t <= T; t++) 
				dlnt[n][t].setBounds(0,dnt[n][t]);
	NumVarMatrix delta(env,N1);  // phase angles, delta0
		for(IloInt n = 0; n < N; n++)
				delta[n] = IloNumVarArray(env,T+1, -maxdelta, maxdelta);	
	NumVarMatrix ucd(env,N);  // demand undertainty
		for(IloInt n = 0; n < N; n++)
				ucd[n] = IloNumVarArray(env,T+1, 0, IloInfinity);	
	BoolVarMatrix card(env,N1);	// demand take UB (1)
		for(IloInt n = 0; n < N; n++)
				card[n] = IloBoolVarArray(env,T+1);	*/
 //------------------------------------Intermediate values
			NumMatrix ustar(env,J);  // start-up, best feasible solution so far
			for(IloInt i = 0; i < J; i++)
				ustar[i] = IloNumArray(env,T+1,0,1); 
			NumMatrix vstar(env,J);  // on-off, best feasible solution so far
			for(IloInt i = 0; i < J; i++)
				vstar[i] = IloNumArray(env,T+1,0,1);
			NumMatrix wstar(env,J);  // best feasible solution so far
			for(IloInt i = 0; i < J; i++)
				wstar[i] = IloNumArray(env,T+1,0,1); 
			NumMatrix ucstar(env,J);  //  current iteraion
			for(IloInt i = 0; i < J; i++)
				ucstar[i] = IloNumArray(env,T+1,0,1); 
			NumMatrix vcstar(env,J);  //  current
			for(IloInt i = 0; i < J; i++)
				vcstar[i] = IloNumArray(env,T+1,0,1); 
			NumMatrix wcstar(env,J);  //  current 
			for(IloInt i = 0; i < J; i++)
				wcstar[i] = IloNumArray(env,T+1,0,1); 
			int U[J][T+1], V[J][T+1],W[J][T+1]; // for calub func
	 for(int t = 0; t <= T; t++)
		for(int j = 0; j < J; j++){
			U[j][t] = 0;
			V[j][t] = 0;
			W[j][t] = 0;
			}    
			 IloNumArray dtstar(env,T);
		
		/************** Constraints******************/
//-----------------------------------------model
IloModel model1(env);
IloNumVar beta(env, -IloInfinity,IloInfinity);
model1.add(IloMinimize(env,beta));
//-----------Initial status
for(IloInt j=0; j<J; j++){
	for(IloInt it=0;it<it_out;it++)
	model1.add(g[it][j][0]==0);
	model1.add(u[j][0]==uj0[j]);
	model1.add(v[j][0]==0);
	model1.add(w[j][0]==0);
	gj0[j]=0;
	}
for(IloInt j=0; j<J; j++){ //initial min up/down
	UTr[j]=((UT[j]-UT0[j])*uj0[j]>0)?((UT[j]-UT0[j])*uj0[j]):0;
	DTr[j]=((DT[j]-DT0[j])*(1-uj0[j])>0)?((DT[j]-DT0[j])*(1-uj0[j])):0;
	for(IloInt t=1; t<=UTr[j]+DTr[j]; t++){ 
		model1.add(u[j][t]==uj0[j]);
		}
 /*	if(DT0[j]>=2)
	for(IloInt s=1; s<Sj[j]; s++) 
		for(IloInt t=Tsu[j][s+1]-DT0[j]+1; t<Tsu[j][s+1]; t++){ //initial start-up type	
			model1.add(h[j][s][t]==0);
		} */
	}// end for j
// logic relationship
	for(IloInt j=0; j<J; j++)
		for(IloInt t=1; t<=T; t++){ 
					model1.add(u[j][t]-u[j][t-1]==v[j][t]-w[j][t]);	
					}						
// min  Up time 
	for(IloInt j=0; j<J; j++)
		for(IloInt t=UT[j]+1; t <=T; t++){
			IloExpr Expjut(env);
			for(IloInt i=t-UT[j]+1; i<=t; i++){
			Expjut+=v[j][i];
			}
		model1.add(Expjut<=u[j][t]);
		Expjut.end();
		}	  	
// min  Down time 
	for(IloInt j=0; j<J; j++)
		for(IloInt t=DT[j]+1; t <=T; t++){
			IloExpr Expjdt(env);
			for(IloInt i=t-DT[j]+1; i<=t; i++){
			Expjdt+=v[j][i];
			}
		model1.add(Expjdt<=1-u[j][t]);
		Expjdt.end();
		}        

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% constraint for each scenario %%%%%%%%%%%%%%%%%%%%%%%%%%%%//				
IloInt itrc_out = 0;
IloNum outgap = 1;
IloNum outLB = -IloInfinity;
IloNum outUB = IloInfinity;

 while (outgap > outgap0){
	 itrc_out++;
// updating for each run
for(IloInt t = 1; t <= T; t++){
	Dt[t]=0;
	for(IloInt n = 0; n < N1; n++){
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
		}
	 for(IloInt j=0; j<J; j++)
		for(IloInt t=1; t<=T; t++){	 
		IloExpr Expub(env);
			for(IloInt i=1; i <= Kjt[j][t];i++){
				Expub+=(SD[j]+(i-1)*RD[j])*w[j][t+i]-Gjb[j]*v[j][t+i];
			//	Expub+=(i*RD[j])*w[j][t+i]-Gjb[j]*v[j][t+i];
				}
		 // model1.add(g[j][t] + reserve[j][t]<= Gjb[j]*u[j][t+Kjt[j][t]]+Expub);
			model1.add(g[itrc_out][j][t]<= Gjb[j]*u[j][t+Kjt[j][t]]+Expub);
		Expub.end();
		}    
// Generation lower bound limit
	for(IloInt j=0; j<J; j++)
		for(IloInt t= 1; t <= T; t++){
			model1.add(g[itrc_out][j][t]>=u[j][t]*Gj[j]);
		//	model1.add(g[j][t]<=u[j][t]*Gjb[j]);
			} 
 
// Ramping down ??? t=1,2,3 ???
 for(IloInt j=0; j<J; j++){
		if( RD[j] > SU[j]-Gj[j] && UT[j] >= 2){ 
			for(IloInt t=1; t<=T; t++){
				model1.add(g[itrc_out][j][t-1]-g[itrc_out][j][t] <= RD[j]*u[j][t]+SD[j]*w[j][t]-(RD[j]-SU[j]+Gj[j])*v[j][t-1]-(RD[j]+Gj[j])*v[j][t]);
			}
		}
		else{ 
			 for(IloInt t=1; t<=T; t++)
				 model1.add(g[itrc_out][j][t-1]-g[itrc_out][j][t] <= RD[j]*u[j][t]+SD[j]*w[j][t]);
			}
	}		  
	for(IloInt j=0; j<J; j++){
		if( RD[j] > SU[j]-Gj[j] && UT[j]>=3 && DT[j]>=2){ 
			for(IloInt t=1; t<=T-1; t++){
				model1.add(g[itrc_out][j][t-1]-g[itrc_out][j][t] <= RD[j]*u[j][t+1]+SD[j]*w[j][t]+RD[j]*w[j][t+1]-(RD[j]-SU[j]+Gj[j])*v[j][t-1]-(RD[j]+Gj[j])*v[j][t]-RD[j]*v[j][t+1]);
			//	model1.add(g[itrc_out][j][t-1]-g[itrc_out][j][t] <= RD[j]*u[j][t]+SD[j]*w[j][t]-(RD[j]-SU[j]+Gj[j])*v[j][t-1]-(RD[j]+Gj[j])*v[j][t]); //modified based on logic constraints
			}
		}
	}	 	 
	for(IloInt j=0; j<J; j++) //over two periods
		if( UT[j] >= 2)
			{ 
			for(IloInt t=2; t<=T; t++) {
				// from Jim's feedback:
				model1.add(g[itrc_out][j][t-2]-g[itrc_out][j][t] <=2*RD[j]*u[j][t]+SD[j]*w[j][t-1]+(SD[j]+RD[j])*w[j][t]+(-2*RD[j]+SU[j]-Gj[j])*v[j][t-2]-(2*RD[j]+Gj[j])*v[j][t-1]-(2*RD[j]+Gj[j])*v[j][t]);
				// new
			//	model1.add(g[itrc_out][j][t-2]-g[itrc_out][j][t] <=2*RD[j]*u[j][t]+(1*SD[j])*w[j][t-1]+2*SD[j]*w[j][t-2]+(SD[j]+RD[j])*w[j][t]+(-2*RD[j]+SU[j]-Gj[j])*v[j][t-2]-(2*RD[j]+Gj[j])*v[j][t-1]-(2*RD[j]+Gj[j])*v[j][t]);
				}
			} 	 		
 //Ramping down----traditional 
for(IloInt j=0; j<J; j++)
	for(IloInt t=1; t<=T; t++){
//		model1.add(g[itrc_out][j][t-1]-g[itrc_out][j][t] <= RD[j]*u[j][t]+SD[j]*w[j][t]); 
		} 

// Ramping up 
	for(IloInt j=0; j<J; j++){
		if( RU[j] > SD[j]-Gj[j] && UT[j] >= 2){ 
			for(IloInt t=1; t<=T-1; t++){
				model1.add(g[itrc_out][j][t]-g[itrc_out][j][t-1] <= RU[j]*u[j][t]-Gj[j]*w[j][t]-(RU[j]-SD[j]+Gj[j])*w[j][t+1]+(SU[j]-RU[j])*v[j][t]);
			//	model1.add(g[itrc_out][j][t]+reserve[itrc_out][j][t]-g[itrc_out][j][t-1] <= RU[j]*u[j][t]-Gj[j]*w[j][t]-(RU[j]-SD[j]+Gj[j])*w[j][t+1]+(SU[j]-RU[j])*v[j][t]);
			}
		}
		else{
			for(IloInt t=1; t<=T; t++){
				model1.add(g[itrc_out][j][t]-g[itrc_out][j][t-1] <= RU[j]*u[j][t-1]+SU[j]*v[j][t]);
			//	model1.add(g[j][t]+reserve[j][t]-g[j][t-1] <= RU[j]*u[j][t-1]+SU[j]*v[j][t]);
				}
			}
	}
	for(IloInt j=0; j<J; j++){
		if( RU[j] > SD[j]-Gj[j] && UT[j] >= 2 && DT[j]>=2){ 
			for(IloInt t=2; t<=T-2; t++){
				model1.add(g[itrc_out][j][t]-g[itrc_out][j][t-2] <= 2*RU[j]*u[j][t]-Gj[j]*w[j][t-1]-Gj[j]*w[j][t]+(SU[j]-RU[j])*v[j][t-1]+(SU[j]-2*RU[j])*v[j][t]);
			//	model1.add(g[itrc_out][j][t]+reserve[itrc_out][j][t]-g[itrc_out][j][t-2] <= 2*RU[j]*u[j][t]-Gj[j]*w[j][t-1]-Gj[j]*w[j][t]+(SU[j]-RU[j])*v[j][t-1]+(SU[j]-2*RU[j])*v[j][t]);
			}
		}
	}  	
//Ramping Up----traditional
for(IloInt j=0; j<J; j++)
	for(IloInt t=1; t<=T; t++){
	//	model1.add(g[itrc_out][j][t]-g[itrc_out][j][t-1] <= RU[j]*u[j][t-1]+SU[j]*v[j][t]);
	//	model1.add(g[itrc_out][j][t]<=goverline[j][t]);
	//	model1.add(goverline[itrc_out][j][t]<=Gjb[j]);
		}
//start-up cost 
 for(IloInt j=0; j<J; j++)
		for(IloInt s=0; s<Sj[j]-1; s++) // last period not included s [1,Sj) ->[0,Sj-1)
			for(IloInt t=Tsu[j][s+1]; t<=T; t++){
				IloExpr expsc(env);
				for(IloInt i=Tsu[j][s]; i<=Tsu[j][s+1]-1; i++){
					expsc+= w[j][t-i];
					}
				model1.add(h[j][s][t] <= expsc);
				expsc.end();
				}				
for(IloInt j=0; j<J; j++)
		for(IloInt t=1; t<=T; t++){
			IloExpr expsc1(env);
			for(IloInt s=0; s<Sj[j]; s++)
				expsc1+=h[j][s][t];
			model1.add(expsc1==v[j][t]);
			expsc1.end();
			}	
		
//  Spinning reserve	
for(IloInt t=1; t<=T; t++) {
		IloExpr expsr1(env);
		for(IloInt j=0; j<J; j++){
			//expsr1+=reserve[j][t];
			expsr1+=g[itrc_out][j][t];
			}
		//model1.add(expsr1 >=Rt[t]);//+0.01*Rt[t]
		//model1.add(expsr1 >= Dt[t]);//+0.01*Rt[t], Dt[t] is updated 
		expsr1.end();
		}	
//---------------------------Transmission Network Constraints----------------------------------------//
	for(IloInt l = 0; l < L; l++)
		 for(IloInt t = 1; t <= T; t++){
				plt[itrc_out][l][t].setBounds(-Pl[l],Pl[l]); 
		 }
//...... for L lines, DC power flow. 
  {for(IloInt t= 1; t <= T; t++)
	 for(IloInt l=0; l<L;l++)
	 {   
		 m=line[l][0]; //origin bus
		 r=line[l][1]; //destination bus
		 model1.add(plt[itrc_out][l][t]*x[l]==100*(delta[itrc_out][m][t]-delta[itrc_out][r][t]));
	 }}
//........load Balancing constraints, for each bus
   {for(IloInt t= 1; t <= T; t++) 
	for(IloInt n= 0; n < N1; n++)  
	{
		 IloExpr exp60(env);
		 IloExpr exp70(env);
		 IloExpr exp80(env);
		{for(IloInt i= 0; i < J; i++)
		{
			if(Gloc[i]==n){
				exp60+=g[itrc_out][i][t];
				}
		}}
	   {for(IloInt l = 0; l < L; l++)
		{
				if(line[l][1]==n)  // D(l)
					exp70+=plt[itrc_out][l][t];
				else if(line[l][0]==n) // O(l)
					exp80+=plt[itrc_out][l][t];
	   }}
		model1.add(exp60+exp70-exp80+dlnt[itrc_out][n][t]==dnt[n][t]+Zoutstar[itrc_out][n][t]*dpnt[n][t]); //load shed allowed
		//model1.add(exp60+exp70-exp80==dnt[n][t]); // no load shed allowed
		exp60.end();
		exp70.end();
		exp80.end();
	}}//end n	  
 
// highest cost
for(IloInt j=0; j<J; j++)  //initial min up/down
	for(IloInt t=5; t<=T; t++){ 
//		model1.add(u[j][t]==1);
//		model1.add(v[j][t]==1);
		}

//------------------objective function
   IloExpr objfunc(env);
   for(IloInt j=0; j<J; j++)
		for(IloInt t=1; t<=T; t++){
			objfunc+=cj[j]*g[itrc_out][j][t];
			objfunc+=aj[j]*u[j][t];
			for(IloInt s=0; s<Sj[j]; s++)
				objfunc+=Csu[j][s]*h[j][s][t];
		}
   for(IloInt t=1; t<=T; t++)
	   for(IloInt n= 0; n < N1; n++)  
		   objfunc+=ldpen*dlnt[itrc_out][n][t];
   model1.add(beta >= objfunc); 
   objfunc.end();	
//----------- ensure better LB
   model1.add(beta >= outLB);

if(itrc_out>1){
   NumVarMatrix pos(env,J);
	for(IloInt j=0;j<J;j++)
		pos[j]=IloNumVarArray(env,T+1,0,IloInfinity);
	for(IloInt j=0; j<J; j++)
		for(IloInt t=1; t<=T; t++){
			model1.add(pos[j][t] >= u[j][t]-ucstar[j][t]);
			model1.add(pos[j][t] >=-u[j][t]+ucstar[j][t]);
			model1.add(pos[j][t]<=1);
			}
		IloExpr right(env);
	for(IloInt j=0; j<J; j++)
		for(IloInt t=1; t<=T; t++){
			right+=pos[j][t];
			}
//	 model1.add(right>=5);
	pos.end();
	right.end();
	}

// force a different master solution
if(masdif == 1 && itrc_out >=2){
IloExpr ins(env);
for(IloInt j=0; j<J; j++)
	for(IloInt t=1; t<=T; t++)
	{
	ins+=IloAbs(Ustar[j][t]-u[j][t]);
	}
	 model1.add(ins>=1);
ins.end();
}


//--------------------------------------------Solving model--------------------------------------//
//-----------------------timer set-up	
	double totaltime, total_start,total_end;
	total_start=clock();
	IloCplex cplex_model1(model1);
//------------------- warm start-----------------//
//	if(itrc_out>1)
//	   cplex_model1.readMIPStart("outmas.mst");
if(warmst == 1){
	IloNumVarArray startVar(env);
	 IloNumArray startVal(env);
	 if(itrc_out>1){
	 for (int j = 0; j < J; ++j)
		 for (int t = 0; t < T; ++t) {
			 startVar.add(u[j][t]);
			 startVal.add(ucstar[j][t]);
			 startVar.add(v[j][t]);
			 startVal.add(vcstar[j][t]);
			 startVar.add(w[j][t]);
			 startVal.add(wcstar[j][t]);
			 }}
	 if(itrc_out>1)
	 cplex_model1.addMIPStart(startVar, startVal);
	 startVal.end();
	 startVar.end();
	} // end
//	cplex_model1.exportModel("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/uc_milp.lp");
//	cplex_model1.setParam(IloCplex::EpGap, 0.01);
	cplex_model1.setParam(IloCplex::EpGap, outmastol);
	cplex_model1.setParam(IloCplex::MIPSearch,0);  //1, traditional branch and cut, 0 automatic 
//	IloCplex::Callback  callback = cplex_model1.use(lazcut1(env,cplex_model1, u,v,w,Kjt));
//	IloCplex::Callback  callback2 = cplex_model1.use(b(env));
	if (!cplex_model1.solve()) {
		env.error() << "Failed to solve model" << endl;
		throw(-1);
		}
	else
	cout<<"Model Solved!" <<endl;
//	cplex_model1.writeMIPStart("outmas.mst");
// calculate mis-balance
	float total_im = 0;
	float total_flow = 0;
	float max_cutflow = 0;
	for(IloInt t=1; t<=T; t++){
		float g0=0;
		   for(IloInt j=0; j<J; j++){
			 if(label1[Gloc[j]]==0)
				g0+=cplex_model1.getValue(g[itrc_out][j][t]);
		 }
	double d0=0,d1=0;
	for(IloInt n= 0; n < N1; n++){
		if(label1[n]==0)
			d0+=dnt[n][t] ;
	}
//	cout<<"imbalance: "<<g0-d0<<endl;
	total_im+=abs(g0-d0);
	double Fcut=0;
	for(IloInt l=0;l<Betn1;l++){
		for(IloInt la=0;la<L;la++){
			m=line[la][0] ;
			r=line[la][1] ;
			if((m == cut1[l][0] && r == cut1[l][1]) || (r == cut1[l][0] && m == cut1[l][1])) {
			   Fcut+=IloAbs(cplex_model1.getValue(plt[itrc_out][la][t]));
			   total_flow += IloAbs(cplex_model1.getValue(plt[itrc_out][la][t]));
			   max_cutflow += Pl[la] ;
			}
		}
	}
//	cout<<"absolute flow: "<<Fcut<<endl	;
}
//	cout<<"total imbalance: "<<total_im<<endl;
//	cout<<"total imflow: "<<total_flow<<endl;
//	cout<<"max_cutflow: "<<max_cutflow<<endl;

/*	for(IloInt j=0; j<J; j++){ //over two periods
			for(IloInt t=1; t<=T; t++) {
				IloNum tt;
				//tt =( cplex_model1.getValue(g[j][t-2])-cplex_model1.getValue(g[j][t]) <= RD[j]*cplex_model1.getValue(u[j][t+1])+SD[j]*cplex_model1.getValue(w[j][t])+RD[j]*cplex_model1.getValue(w[j][t+1])-(RD[j]-SU[j]+Gj[j])*cplex_model1.getValue(v[j][t-1])-(RD[j]+Gj[j])*cplex_model1.getValue(v[j][t])-RD[j]*cplex_model1.getValue(v[j][t+1]));
				tt = (-cplex_model1.getValue(g[j][t-1])+cplex_model1.getValue(g[j][t]) + RD[j]*cplex_model1.getValue(u[j][t])+SD[j]*cplex_model1.getValue(w[j][t]));
				cout<<tt<< " , "<<endl;
				}
			cout<<endl;
		}*/
	results<<"cplex_model1"<<endl;
	obj_val=cplex_model1.getObjValue();
	outLB = obj_val;
	results<<"obj value: " << double(obj_val)<<endl;
	total_end=clock();
	totaltime = double(total_end-total_start)/CLOCKS_PER_SEC;
	outgap = (outUB-outLB)/outUB;
	ccgresults<<"---outLB: "<<double(obj_val)<<" outGap: "<<outgap<<" outMas time: "<<totaltime<<endl;
	results<<" total time "<<totaltime<<endl;
	{for(IloInt j=0;j<J;j++)
		for(IloInt t=1;t<=T;t++){
			 if(cplex_model1.getValue(u[j][t])<0.1)
						 Ustar[j][t]=0;
			 else
						 Ustar[j][t]=1;;
	}}
	for(IloInt j=0;j<J;j++){
		results<<"unit commitment:"<<j;
		for(IloInt t=1;t<=T;t++){
					 results<<", "<<Ustar[j][t];
				}
		results<<endl;
	}
	float printout;
	for(IloInt j=0;j<J;j++){
		results<<"Generator level:"<<j;
		for(IloInt t=1;t<=T;t++){
			printout=(cplex_model1.getValue(g[itrc_out][j][t])<0.00001)? 0 : cplex_model1.getValue(g[itrc_out][j][t]);
			results<<", "<<printout ;
				} 
		results<<endl;
	}
	for(IloInt l=0;l<L;l++){
		results<<"Power flow:"<<l;
		for(IloInt t=1;t<=T;t++){
					 results<<", "<<((abs(cplex_model1.getValue(plt[itrc_out][l][t]))<0.000001)?0:cplex_model1.getValue(plt[itrc_out][l][t]));
					 flows<<abs(cplex_model1.getValue(plt[itrc_out][l][t]))/Pl[l]<<"	";
				}
		results<<endl;
		flows<<endl;
	}
	double masloadshed=0;
	for(IloInt n=0;n<N1;n++){
		results<<"Load shed:"<<n;
		for(IloInt t=1;t<=T;t++){
					 results<<", "<<((cplex_model1.getValue(dlnt[itrc_out][n][t])<0.00001)?0:cplex_model1.getValue(dlnt[itrc_out][n][t]));
					 masloadshed+=cplex_model1.getValue(dlnt[itrc_out][n][t]);
				}
		results<<endl;
	}
	objval<<"load shedding in mas problem: "<<masloadshed<<" cost:"<<ldpen*masloadshed<<endl;
//--------------bounding info
	for(IloInt l=0;l<L;l++){
		Fbounds[l]=IloMin(200*maxdelta/x[l], Pl[l]); // bound for power flow
		}
	for(IloInt l=0;l<L;l++){
		duals[l]=0; // initialized to 0
		}
	for(IloInt l=0;l<L;l++){
		for(IloInt t=1;t<=T;t++){
					if( Fbounds[l]-abs(cplex_model1.getValue(plt[itrc_out][l][t]))<=0.0000001 )
						duals[l]=duals[l]+1;
				}
		Duals<< duals[l] <<endl;
	}

//--------------update solution
/*	for(IloInt j=0;j<J;j++)	{
		   cplex_model1.getValues(ucstar[j],u[j]);
		   cout<< "ucstar[j]: "<<ucstar[j]<<endl;
		   cplex_model1.getValues(vcstar[j],v[j]);
			cout<< "vcstar[j]: "<<vcstar[j]<<endl;
		   cplex_model1.getValues(wcstar[j],w[j]);
			cout<< "wcstar[j]: "<<wcstar[j]<<endl;
	}  */
	int difuc = 0;// number of different UC decisions
	for(IloInt j=0;j<J;j++)
		for(IloInt t=0;t<=T;t++){
		  if ( cplex_model1.getValue(u[j][t]) < 0.01 && U[j][t]>0.01 ){
			difuc++;
			genu<<" 0->1 ";
			}
		 else if(cplex_model1.getValue(u[j][t]) > 0.01 && U[j][t]<0.01){
			difuc++;
			genu<<" 1->0 ";
			 }
			}
	genu<<"# of different UC decisions: "<<difuc<<endl;
	
	for(IloInt j=0;j<J;j++)
		for(IloInt t=0;t<=T;t++){
		  if ( cplex_model1.getValue(u[j][t]) < 0.01){
			  ucstar[j][t] =0;
			  U[j][t]=0;
			}
		  else {ucstar[j][t] =1;U[j][t]=1;}
		  if ( cplex_model1.getValue(v[j][t]) < 0.01){
				vcstar[j][t] =0;
				V[j][t]=0;}
		  else {vcstar[j][t] =1;V[j][t]=1;}
		  if ( cplex_model1.getValue(w[j][t]) < 0.01){
			   wcstar[j][t] =0; 
			   W[j][t]=0;}
		  else {wcstar[j][t] =1;W[j][t]=1;}
	}// end for J
//**** verify solution differences
/*	genu<<"UC solution:"<<endl;
	for(IloInt j=0;j<J;j++)
			for(IloInt t=0;t<=T;t++){
				genu<<U[j][t]<<" ";
				}
	genu<<endl;	   */
	costuvw=0;
	for(IloInt j=0; j<J; j++)
		for(IloInt t=1; t<=T; t++){
			costuvw +=aj[j]*cplex_model1.getValue(u[j][t]);
			for(IloInt s=0; s<Sj[j]; s++)
				costuvw += Csu[j][s]*cplex_model1.getValue(h[j][s][t]);
		}
		cout<<"first part cost: "<<costuvw<<endl;
		objval<<"costuvw: "<<costuvw<<endl;
// free master model and algorithm
	cplex_model1.end();


// order switched sub problem 
init();
double ingap;
for(IloInt t=1; t<=T; t++)
	{for(IloInt n= 0; n < N1; n++)  
	cout<<head->Z[n][t]<<",";
cout<<endl;
}
double subtime, sub_start,sub_end;
sub_start=clock();
calub(head->ub,head->lb,head->Z,U,V,W,head->Znode,uncbg,costuvw);
sub_end=clock();
subtime = double(sub_end-sub_start)/CLOCKS_PER_SEC;
ingap = (head->ub-head->lb)/head->ub;

if(ingap>BB_gap0 && brahopt ==1 ){
	float B_gap;
	leaf=0;
	NewNode(head); // generate child nodes and update fixed variables
	node *N;
	N=head->next;
	while(1)
	{
		N=head->next;
		if(N == NULL)
		{
			cout<<"No nodes in linked list, tree level"<<Level<<"!"<<std::endl;
			break;
		}
		//	time(&f_stop);
		//	if(difftime(f_stop, f_start) > 36000 ) // in seconds
		//	{
		//	cout<<"Time Limit 36000s"<<endl;
		//	cout<<"**************"<<endl;
		//	break;
		//	}
//************* visit node in the linked list, cal(); 
		while(N != NULL)//evaluate all nodes of level L
		{
		// update gap with new upper bound .........
			calub(N->ub,N->lb,N->Z,U,V,W,N->Znode,uncbg,costuvw);
			cout<<"ub updated!"<<endl;
			leaf += 1;// record the # of nodes that are evaluated
		//	feasible solution of current node
			if(N->lb>BB_LB){
				for(int t=1; t<=T; t++)
					for(int n= 0; n < N1; n++){ 
						Zb[n][t] = N->Znode[n][t]; // update best feasible solution
						}
			//	BB_LB = N->lb;
				}
		//	BB_UB= max(N->ub,BB_UB); 
			BB_LB= max(N->lb,BB_LB); 
			BB_gap = abs((BB_UB-BB_LB)/BB_UB);
			cout<<"BB_UB:"<<BB_UB<<" BB_LB:"<<BB_LB<<" BB_gap:"<<BB_gap<<endl;
			BBres<<" BB_UB="<<BB_UB<<" BB_LB="<<BB_LB<<" BB_gap="<<BB_gap*100<<"%"<<" node updating, leaf:"<<leaf<<endl;

		
			node_gap= (N->ub-N->lb)/N->ub; 
			if(node_gap < node_gap0) // optimal leaf node
			{
				BBres<<"Node optimal"<<endl;
				N->tag = 2; // tag as optimal
				N = N->next;
			}
			else // non optimal,proceed to next node.
			{
				N->tag = 1;		
				N = N->next; // go to next node
			}
		} // end while when linked list is empty

//****************update the bounds for each node according to new evaluated nodes************
		BB_gap=update();
		// if overall BB_gap small enough, stop; else branch the active nodes and continue 
		

		BBres<<" BB_UB="<<BB_UB<<" BB_LB="<<BB_LB<<" BB_gap="<<BB_gap*100<<"%"<<" after updating tree, leaf:"<<leaf<<endl;
		if(BB_gap < BB_gap0) // optimal
		{
			cout<<"########### overall BB_gap small enough: "<<"BB_UB:"<<BB_UB<<" BB_LB:"<<BB_LB<<" BB_gap:"<<(BB_UB-BB_LB)*100/BB_LB<<"%"<<endl;
			cout<<"**************"<<endl;
			break; // break if whole b&b tree is optimal
		}
		else  //********b&b non-optimal, dequeue whole level nodes and add their child nodes, 
			  //************if tag = 1, add its left and right child to queue, dequeue itself.
		{
			node *h = head->next; // head node stays
			while(h != NULL)
			{		
				if( h->tag == 1 ) // evaluated, vars left to be fixed
					{
						NewNode(h); //put the two child nodes in the front of the queue
						DeQueue(h); // deque the 
					}
				else // pruned or cannot be branched any more
					{
			//		NewNode(h);
						DeQueue(h);
					}
			h=h->next;
			} // end while h
			Level += 1;
		} // end else *
	
	} // end while(1), whole loop
	cout<<"branch and bound performed"<<endl;
	ccgresults<<"branch and bound performed"<<endl;
}// end for bb

	


if(costuvw+head->lb < outUB){
for(IloInt t=1; t<=T; t++)
	for(IloInt n= 0; n < N1; n++){
		Zoutstar[itrc_out][n][t]=head->Znode[n][t]; // update best feasible solution of subproblem
		Zd[n][t] = head->Znode[n][t];
		}
	}

int dediff=0;
for(IloInt t=1; t<=T; t++)
		for(IloInt n= 0; n < N1; n++){
			if(itrc_out>=2){
				if(abs(Zoutstar[itrc_out-1][n][t]-Zoutstar[itrc_out][n][t]) >0.1  )
					dediff++;
				}
			}
		genu<<" sub difference: "<< dediff <<endl;

outUB = IloMin(costuvw+head->lb,outUB);
objval<<" current lb: "<<costuvw+head->lb<<" loadshed: "<<(head->lb)/ldpen<<endl;
outgap = (outUB-outLB)/outUB;
ccgresults<<"---outUB: "<<outUB<<" outGap: "<<outgap <<" sub prob time: "<<subtime<<endl;



} // end for outer while




	model1.end();
}// end of try block
	catch (IloException& ex){	  cerr << "Error: " << ex << endl;	}
	catch (...)				{	  cerr << "Error" << endl;			}
	env.end();
	return 0;
}//end of main

