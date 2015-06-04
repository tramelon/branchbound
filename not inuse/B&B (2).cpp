struct node{
	 //used in the tree
     node *parent, *l, *r; //record the parent and child nodes of current node
	 int level; //depth of current node
	 int tag;// -1: pruned by bound 
					//0: not evaluated 
					//1:evaluated, but BB_gap>1% 
					//2: optimal, pruned 
					//3: BB_gap>1%, but can not be branched any more
	 // used in the chain
     node *next; 
	 //
	 
     float lb;    
     float ub;   

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
void init()//
{
  
  leaf = 0;
  BB1=0;
  L=1;

  head=new(node);
  head->next=NULL;
  head->parent = NULL;

  head->l = NULL;
  head->r = NULL;

  head->level = 0;
  head->tag = 1; 

  
  //spicify the initial values of ub, lb, multipliers/Benders cuts.... of the root node,

  

}

void NewNode(node *prnt) //generate new nodes from prnt
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
    // determine the branching variable  b
	//.......
	////////////////////////////////////////////////////////////////////////////////////////
    ir->P[b] = 0; // right child
	for(int j=0; j<nNode; j++)
		if(j != b)
			ir->P[j] = prnt->P[j];

	for(int i=0; i<nNode; i++)
		for(int j=0; j<nNode; j++)
		{
			ir->Alpha1[i][j] = prnt->Alpha1[i][j];
			ir->Alpha2[i][j] = prnt->Alpha2[i][j];
		}

	for(int i=0; i<nNode; i++)
		for(int j=0; j<nNode; j++)
			for(int k=0; k<nNode; k++)
			{
				ir->Gamma1[i][j][k] = prnt->Gamma1[i][j][k];
			    ir->Gamma2[i][j][k] = prnt->Gamma2[i][j][k];
			}
	ir->ub = prnt->ub;
	ir->lb = prnt->lb;

	if(ir->lb <= BB_UB)
		EnQueue(ir);

    il->P[b] = 1; // left child
	for(int j=0; j<nNode; j++)
		if(j != b)
			il->P[j] = prnt->P[j];

	int aux_nH = 0;
	for(int k=0;k<nNode;k++)
		if(il->P[k] == 1)
			aux_nH += 1;

	if(aux_nH == p)
	{
		for(int k=0;k<nNode;k++)
			if(il->P[k] != 1)
				il->P[k] = 0;

	}



	///////carry multipliers from their parent node
	for(int i=0; i<nNode; i++)
		for(int j=0; j<nNode; j++)
		{
			il->Alpha1[i][j] = prnt->Alpha1[i][j];
			il->Alpha2[i][j] = prnt->Alpha2[i][j];
		}

	for(int i=0; i<nNode; i++)
		for(int j=0; j<nNode; j++)
			for(int k=0; k<nNode; k++)
			{
				il->Gamma1[i][j][k] = prnt->Gamma1[i][j][k];
			    il->Gamma2[i][j][k] = prnt->Gamma2[i][j][k];
			}
	il->ub = prnt->ub;
	il->lb = prnt->lb;

	//enter the chain

	if(il->lb <= BB_UB)
		EnQueue(il);
}

node *NextLiveNode() //return next live node in the chain
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

void update0(node *N) // see update()
{
	if(N->l != NULL && N->r != NULL)
	{
		if(N->l->ub < N->r->ub)
			N->ub = N->l->ub;
		else 
			N->ub = N->r->ub;
		
		if(N->l->lb < N->r->lb)
			N->lb = N->l->lb;
		else 
			N->lb = N->r->lb;
	}
	if(N->l != NULL && N->r == NULL)
	{
			N->ub = N->l->ub;
			N->lb = N->l->lb;
	}
	if(N->l == NULL && N->r != NULL)
	{
			N->ub = N->r->ub;
			N->lb = N->r->lb;
	}

	if((N->ub - N->lb)/N->lb < BB_gap0)
	{
		N->tag = 2;
		//N->lb = N->ub;
	}



}
void findo(node *N, int l)//see update()
{
	
	node *N1, *N2;

    if(N->l != NULL && N->r != NULL){
	if(N->l->tag == 1 && N->r->tag == 1)
	{
		if(N->l->level == l)
		{
			N1=N->l;
			N2=N->r;
			update0(N1);
			update0(N2);

		}
		else
		{
			findo(N->l,l);
			findo(N->r,l);

		}

	}
	if(N->l->tag != 1 && N->r->tag == 1)
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
	if(N->l->tag == 1 && N->r->tag != 1)
	{
		if(N->l->level == l)
		{
			N1=N->l;
			update0(N1);

		}
		else
		{
			findo(N->l,l);

		}

	}}


	if(N->l != NULL && N->r == NULL)
	{
		if(N->l->tag == 1 )
		{
			if(N->l->level == l)
			{
				N1=N->l;
				update0(N1);
			}
			else
			{
				findo(N->l,l);

			}

		}
	}

	if(N->r != NULL && N->l != NULL)
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





}
void prune(node *N) // cut off the nodes that do not need further evaluation
{
	if(N->lb > BB_UB)
		N->tag = -1;
	if(N->l != NULL && N->r != NULL)
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
float update() // update()+findo()+update0(): update the lb and ub after all nodes in a level are evaluated; update() recursively calls findo() and update0()
{
	float g;
	for(int l=L-1; l>0; l--)
	{
		findo(head, l);

	}
	update0(head);
	if(BB_UB > head->ub)
		BB_UB = head->ub;


	g = (BB_UB - head->lb)/BB_UB;
	
	BB_LB=head->lb;
	cout<<"BB_UB="<<BB_UB<<" BB_LB="<<BB_LB<<" BB_gap="<<g*100<<"%"<<endl;
	cout<<"************************"<<endl;
	if(g<BB_gap0)
		return g;
	else
	{
		prune(head);
		return g;

	}

}


void BB_main() // main funciton of Branch and Bound
{

	
	int n0,n1;
	node *k;
	int b;
	
	//main body
	init();
	NewNode(head);
	node *N;

	float B_gap;
	

	
	while(1)
	{
		N=head->next;
		if(N == NULL)
		{
			cout<<"No nodes of Level:"<<L<<"!"<<endl;
			break;
		}


		while(N != NULL)//evaluate all nodes of level L
		{
			BB_gap = calculate(N);

			leaf += 1;// record the # of nodes that are evaluated

			time(&f_stop);
			if(difftime(f_stop, f_start) > 3600 )
			{
				break;
			}

			if(BB_gap < BB_gap0)
			{
				cout<<"leaf node #"<<leaf<<" evaluated, BB_gap:"<<BB_gap*100<<"%"<<endl;
				cout<<"ub: "<<N->ub<<"--"<<"lb: "<<N->lb<<endl;
				for(int k = 0; k<nNode; k++)
					if(N->P[k] != -1)
						cout<<k<<":"<<N->P[k]<<" ";
				cout<<endl;

				//N->lb = N->ub;
				N->tag = 2;

				cout<<"**************"<<endl;

				
				N = N->next;
			}
			else
			{
				cout<<"node #"<<leaf<<" evaluated, BB_gap:"<<BB_gap*100<<"%"<<endl;
				for(int k = 0; k<nNode; k++)
					if(N->P[k] != -1)
						cout<<k<<":"<<N->P[k]<<" ";
				cout<<endl;
			
				N->tag = 1;
			
				cout<<"ub: "<<N->ub<<"--"<<"lb: "<<N->lb<<endl;
            
			
				cout<<"**************"<<endl;
			
				N = N->next;
			}
		}

		//update the bounds for each node according to new evaluated nodes
		BB_gap=update();
		// if overall BB_gap small enough, stop; else brach the active nodes and continue 
		time(&f_stop);
		if(difftime(f_stop, f_start) > 3600)
		{
			cout<<"Time Limit 3600s"<<endl;
			cout<<"**************"<<endl;
			break;
		}
		if(BB_gap < BB_gap0)
		{
			cout<<"########### overall BB_gap small enough: "<<"BB_UB:"<<BB_UB<<" BB_LB:"<<BB_LB<<" BB_gap:"<<(BB_UB-BB_LB)*100/BB_LB<<"%"<<endl;
			cout<<"**************"<<endl;
			break;
		}
		else
		{
			node *h = head->next;
			while(h != NULL)
			{
				n0=0,n1=0;
				for(int i=0; i<nNode; i++)
				{
					if(h->P[i]==1)
						n1 += 1;
					if(h->P[i]==0)
						n0 += 1;
				}
			
			
				if(h->tag == 1 && n0 < nNode-p && n1 < p)
				{
					k=h->parent;
				    b=1;
					while(k != head)
					{
						if(k->tag != 1)
						{
							b=k->tag;
							break;
						}
						k = k->parent;

					}
				
					if(b==1)
					{
						NewNode(h); //put the two child nodes in the front of the queue
						DeQueue(h);
					}
					else
					{
						h->tag = b;
						DeQueue(h);

					}
				}
				else
				{
					if(n0 == nNode-p && n1 == p)
						h->tag = 3; 
					DeQueue(h);
				}

				h=h->next;


			}
			L += 1;

		}
	
	}

    
    cout<<"######################################################"<<endl;
	cout<<"# of nodes checked: "<<leaf<<endl;
	BB1 = leaf;
	
}
