

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
void EnQueue(node *i) //enter the chain
void DeQueue(node *i)// leave the chain
void init()// initialization
void NewNode(node *prnt) //generate new child nodes from prnt
node *NextLiveNode() //return next live node in the chain, not used in program
void update0(node *N) // see update(), update bounds when bounding info computed, for a single node and its child nodes
void findo(node *N, int l)//see update(), recursively call findo, level l
void prune(node *N) // cut off the nodes that do not need further evaluation, prune whole branch 
float update() // update()+findo()+update0(): 

