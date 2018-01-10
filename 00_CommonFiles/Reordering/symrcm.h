#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<time.h>

// A node struct for the Cuthill-McKee algorithm
typedef struct CMK_Node
{
	// the node's id (matrix row index)
	int id;
	// the node's degree
	int deg;
     	// minimal distance to the root of the spanning tree
     	int dist;
}CMK_NodeType;

// Predicate (queue empty)
#define Q_empty(Q, N, qh, qt)   ((qh) == (qt))

// Predicate (heap empty)
#define H_empty(H, h)   ((h) == 0)

// A simple, array-based binary heap (used as a priority queue for nodes)
// the left descendant of entry i
#define LEFT(i)         (((i) << 1) + 1)        // = (2*(i)+1)
// the right descendant of entry i
#define RIGHT(i)        (((i) << 1) + 2)        // = (2*(i)+2)
// the parent of entry i
#define PARENT(i)       (((i) - 1) >> 1)        // = floor(((i)-1)/2)

// Swap to variables independently of type
#define SWAP(a, b) do { typeof(a) temp = a; a = b; b = temp; } while (0)

void transform_CSR_to_CCS(int N,const int *ja, const int *ia, int *cidx,int *ridx);
void transpose (int N, const int *ridx, const int *cidx, int *ridx2, int *cidx2);
int calc_degrees(int N, const int *ridx, const int *cidx, int *D);
int find_starting_node (int N, const int *ridx, const int *cidx, const int *ridx2, const int *cidx2, int *D, int start);
void Q_enq (CMK_NodeType *Q, int N, int *qt, const CMK_NodeType *o);
CMK_NodeType Q_deq (CMK_NodeType *Q, int  N, int *qh);
void H_insert (CMK_NodeType *H, int *h, const CMK_NodeType *o);
CMK_NodeType H_remove_min (CMK_NodeType *H, int *h, int reorg);
void H_heapify_min (CMK_NodeType *A, int i, int size);


