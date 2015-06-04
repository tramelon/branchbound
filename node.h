
#ifndef NODE_H
#define NODE_H


#ifndef DEFINEDMATRIX_H
#define DEFINEDMATRIX_H
#include "definedmatrix.h"
#endif DEFINEDMATRIX_H




int leaf ;
int BB1;
int Level;  // level of bb tree

double BB_UB;   // update when a whole level of leaf node is processed
double BB_LB;   // best feasible solution bound, updated every node
double BB_gap0; // whole process gap
double BB_gap;  // branch and bound tree gap, in head node

struct node{
	 //used in the tree
	 node *parent, *l, *r; //record the parent and child nodes of current node
	 int level; //depth of current node, 0,1,...
	 int tag;// -1: pruned by bound
		//0: not evaluated
		//1:evaluated, but BB_gap>1%
		//2: optimal, pruned
		//3: BB_gap>1%, but can not be branched anymore,based on the constraints
	 node *next;
	 double lb;
	 double ub;
	 int Z[N1][T1+1]; // var to be fixed
	 float Znode[N1][T1+1]; // best feasible solution for current node
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
//determine the branching variable  z[a][b], first non -1 to branch;
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
} // prune

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
}// update




#endif