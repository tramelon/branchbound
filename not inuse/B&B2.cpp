

#include <stdio.h> 
#include <stdlib.h> 

#define MAXN       15 
#define MAXENTRIES 1000 
#define FALSE      0 
#define TRUE       1 

typedef struct node *link; 
struct node { 
  int nunsched;  // number of jobs still unscheduled 
  int *unsched;  // array to tell which jobs still unscheduled 
  int *sched;    // partial schedule 
  int tunsched;  // total processing time of unscheduled jobs 
  int lb;        // lower bound -- penalty of scheduled jobs 
}; 

//*************************************************************** 
// priority queue 
link pq[1+MAXENTRIES]; 
int pqsize;      // current size 

void pqinit() 
{ 
  // initalize priority queue -- put dummy entry at position 0 
  link t; 
  pqsize = 0; 
  t = new node; 
  t->lb = 0; t->nunsched = 0; 
  pq[0] = t; 
} 

void swap(link& p, link& q) 
{ 
  link temp; 
  temp = p; p = q; q = temp; 
} 

int pqempty() 
{ 
  return pqsize==0; 
} 

void pqfull() 
{ 
  printf("priority queue full -- program aborted\n"); 
  exit(0); 
} 

int lt(link p, link q) 
{ 
  return ( p->lb<q->lb 
    || (p->lb==q->lb && p->nunsched<q->nunsched) ) ; 
} 

void pqinsert(link p) 
{ 
  int i, j; 
  if (pqsize>=MAXENTRIES) pqfull(); 
  pq[++pqsize] = p; 
  i = pqsize; j = i / 2; 
  while ( lt(pq[i],pq[j]) ) { 
    swap(pq[i], pq[j]); 
    i = j; j = j / 2; 
  } 
} 

void deletemin(link & p) 
{ 
  // will only be called with pqsize>0 
  int i, j; 
  p = pq[1]; pq[1] = pq[pqsize--]; 
  i = 1; j = 2; 
  while ( j<=pqsize ) { 
    if (j<pqsize) { 
      if ( lt(pq[j+1],pq[j]) ) j++; 
    } 
    if ( !lt(pq[j],pq[i]) ) break; 
    swap(pq[i], pq[j]); 
    i = j; j = 2*i; 
  } 
} 
//*************************************************************** 

void main() 
{ 
  int proctime[1+MAXN],   // processing times 
      deadline[1+MAXN],   // deadlines 
      curbestsched[1+MAXN];  // current best schedule 

  int njobs,       // number of jobs, not to exceed MAXN 
      tottime,     // total processing time of all jobs 
      curbestpen;  // penalty of current best 

  int i, j; 
  link p, q; 

  while (printf("#jobs (0 to quit)? "), scanf("%d", &njobs), njobs > 0)  { 
    if (njobs > MAXN) { 
      printf("#jobs must be in range %d..%d\n", 1, MAXN); 
      continue; 
    } 

    tottime = 0; curbestpen = 0; 
    printf("Enter processing times: "); 
    for (i = 1; i<=njobs; i++) scanf("%d", &proctime[i]); 
    printf("Enter deadlines: "); 
    for (i = 1; i<=njobs; i++) { 
      scanf("%d", &deadline[i]); 
      tottime += proctime[i]; 
      curbestsched[i] = i; 
      if (tottime>deadline[i]) 
        curbestpen += tottime-deadline[i]; 
    } 

    // create root node & insert into priority queue 
    p = (link) malloc(sizeof(struct node)); 
    p->nunsched = njobs;  p->tunsched = tottime; 
    p->lb = 0; 
    p->unsched = new int[1+njobs]; 
    for (i = 1; i<=njobs; i++) p->unsched[i] = TRUE; 
    p->sched = NULL; 
    pqinit(); pqinsert(p); 

    while (!pqempty()) { 
      deletemin(p); 
      if (p->lb>=curbestpen) { 
        // only pruned nodes left in priority queue -- quit 
        break; 
      } 
      else { 
        for (i = 1; i<=njobs; i++) 
          if (p->unsched[i])  { 
            int thispen = p->tunsched-deadline[i]; 
            if (thispen<0)  thispen = 0; 
            int newpen = p->lb+thispen; 
            if (newpen>=curbestpen);  // throw away 
            else if (p->nunsched==1)  { 
              // complete sched -- update current best 
              for (j = 2; j<=njobs; j++) curbestsched[j] = p->sched[j]; 
              curbestsched[1] = i; curbestpen = newpen; 
            } 
            else { 
              q = (link) malloc(sizeof(struct node)); 
              q->nunsched = p->nunsched-1; 
              q->tunsched = p->tunsched-proctime[i]; 
              q->lb = newpen; 
              q->unsched = new int[1+njobs]; 
              q->sched = new int[1+njobs]; 
              for (j = q->nunsched+2; j<=njobs; j++) 
                q->sched[j] = p->sched[j]; 
              q->sched[q->nunsched+1] = i; 
              for (j = 1; j<=njobs; j++) q->unsched[j] = p->unsched[j]; 
              q->unsched[i] = FALSE; 
              pqinsert(q); 
            } 
          } 
      } 
    } 
    printf("optimal schedule:\n"); 
    for (i = 1; i<=njobs; i++) printf("%4d ", i); 
    printf("\n"); 
    for (i = 1; i<=njobs; i++) printf("%4d ", curbestsched[i]); 
    printf("\n"); 
    printf("penalty = %d\n", curbestpen); 

  } 

  printf("\n"); 
} 

