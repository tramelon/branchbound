// budget for uncertainty  line 330 uncbg = 5

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


ILOSTLBEGIN
//#include "defined_matrix.h"

//global variables

float BB_UB;   // update when a whole level of leaf node is processed
float BB_LB;   // best feasible solution bound, updated every node
float BB_gap0; // whole process gap
float BB_gap;  // branch and bound tree gap, in head node
float node_gap;
float node_gap0; // gap to redem a node to be optimal
int uncbg = 2;  // budget for demand uncertainty
int leaf ;
int BB1;
int Level;  // level of bb tree
#define L 186
#define J 54
#define N1 118
#define T 24
ofstream BBres("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/results/BBres.txt");
float Zb[N1][T+1]; // best feasible solution so far

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
	 float lb;    
	 float ub; 
	 int Z[N1][T+1]; // var to be fixed
	 float Znode[N1][T+1]; // best feasible solution for current node
	
	 // define variables that store the multipliers/Benders cuts
	 // .....
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
		/*int lb=p->lb;
		while(p){
		   if(p->lb>lb){
			  choice=p;
		   }
		   p=p->next;
		}*/
		return(choice);
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

void main() // main function of Branch and Bound
	{				
	//main body
	IloEnv env;
	
try {
	int n0,n1;
	node *k;
	int b; //tag value
	time_t f_stop, f_start; 
// initial on-off status of generators	
	int U[J][T+1], V[J][T+1],W[J][T+1];
	 for(int t = 0; t <= T; t++)
		for(int j = 0; j < J; j++){
			U[j][t] = 0;
			V[j][t] = 0;
			W[j][t] = 0;
			}
	NumMatrix tempu(env,J);
	NumMatrix tempv(env,J);
	NumMatrix tempw(env,J);
	for(int j = 0; j < J; j++){
		tempu[j]=IloNumArray(env,T+1);
		tempv[j]=IloNumArray(env,T+1);
		tempw[j]=IloNumArray(env,T+1);
		} 
readMaData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/U0.txt",tempu, J, T+1);
readMaData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/V0.txt",tempv, J, T+1);
readMaData ("C:/Dropbox/Robust_UC_ISONE/data/118bus/iit/parameters/W0.txt",tempw, J, T+1);
	
		for(int j = 0; j < J; j++){
			for(int t = 0; t <= T; t++){
			U[j][t] = tempu[j][t];
			V[j][t] = tempv[j][t];
			W[j][t] = tempw[j][t];
	//		U[j][t] = 0;
	//		V[j][t] = 0;
	//		W[j][t] = 0;
				}
			}
		for(int j = 0; j < J; j++){
			for(int t = 0; t <= T; t++)
			   cout<<U[j][t]<<",";
			cout<<endl;
			}
				for(int j = 0; j < J; j++){
			for(int t = 0; t <= T; t++)
			   cout<<V[j][t]<<",";
			cout<<endl;
			}
						for(int j = 0; j < J; j++){
			for(int t = 0; t <= T; t++)
			   cout<<W[j][t]<<",";
			cout<<endl;
			}

	init();
	calub(head->ub,head->lb,head->Z,U,V,W,head->Znode,uncbg);
	BB_UB=head->ub;
	BB_LB=head->lb;
	BB_gap=(BB_UB-BB_LB)/BB_UB;


	NewNode(head);
	node *N;
	float B_gap;
	leaf=0;
	time(&f_start);
	while(1&&BB_gap>BB_gap0)
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
			calub(N->ub,N->lb,N->Z,U,V,W,N->Znode,uncbg);
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
}// end try
	catch (IloException& ex){	  cerr << "Error: " << ex << endl;	}
	catch (...)				{	  cerr << "Error" << endl;			}    
	cout<<"######################################################"<<endl;
	cout<<"# of nodes checked: "<<leaf<<endl;
	BB1 = leaf;

}

