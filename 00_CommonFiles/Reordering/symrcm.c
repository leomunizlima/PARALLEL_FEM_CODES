#include "reordering.h"

// An implementation of the Cuthill-McKee algorithm.
void REORDERING_SYMRCM(int n, int nz, int *ja, int *ia, int *p)
{
	// sizes of the heaps
	int s = 0;
  
	// head- and tail-indices for the queue
	int qt = 0;
	int qh = 0;
	CMK_NodeType v, w;

    	// dimension of the matrix
	int N = n;
	int nnz = nz;

	//ATTENTION: I'm admiting that the matrix is structurely symmetric otherwise transform_CSR_to_CCS is necessary

  	int *cidx; // = calloc(N+1, sizeof(int));
  	int *ridx; // = calloc(nnz, sizeof(int));

	
	//transform_CSR_to_CCS(N,ja,ia,cidx,ridx);
	cidx = ia; 
	ridx = ja; 

  	int *cidx2 = calloc(N+1, sizeof(int));
  	int *ridx2 = calloc(nnz, sizeof(int));

    	transpose (N, ridx, cidx, ridx2, cidx2);
  
	// the permutation vector
	int *P = p;
  
	// compute the node degrees
	int *D = calloc(N, sizeof(int));
	int max_deg = calc_degrees (N, ridx, cidx, D);
  
	// if none of the nodes has a degree > 0 (a matrix of zeros)
	// the return value corresponds to the identity permutation
	if (max_deg == 0)
	{
		int i;
		for (i = 0; i < N; i++)
			P[i] = i;
        	return;
      	}
	// a heap for the a node's neighbors. The number of neighbors is
	// limited by the maximum degree max_deg:
	CMK_NodeType *S = calloc(max_deg, sizeof(CMK_NodeType));
  
	// a queue for the BFS. The array is always one element larger than
	// the number of entries that are stored.
	CMK_NodeType *Q = calloc(N+1, sizeof(CMK_NodeType));
  
	// a counter (for building the permutation)
	int c = -1;
  
	// upper bound for the bandwidth (=quality of solution)
	// initialize the bandwidth of the graph with 0. B contains the
	// the maximum of the theoretical lower limits of the subgraphs
	// bandwidths.
	int B = 0;
  
	// mark all nodes as unvisited; with the exception of the nodes
	// that have degree==0 and build a CC of the graph.
  
    	int *visit = calloc(N, sizeof(int)); // visit[0..N-1] receive zero. It means visit_i is false initially
  
	do
	{
		// locate an unvisited starting node of the graph
		int i;
		for (i = 0; i < N; i++)
			if (! visit[i])
				break;
  
		// locate a probably better starting node
		v.id = find_starting_node (N, ridx, cidx, ridx2, cidx2, D, i);
  
		// mark the node as visited and enqueue it (a starting node
		// for the BFS). Since the node will be a root of a spanning
		// tree, its dist is 0.
		v.deg = D[v.id];
		v.dist = 0;
		visit[v.id] = 1; //true
		Q_enq (Q, N, &qt, &v);
  
		// lower bound for the bandwidth of a subgraph
		// keep a "level" in the spanning tree (= min. distance to the
		// root) for determining the bandwidth of the computed
		// permutation P
		int Bsub = 0;
		// min. dist. to the root is 0
		int level = 0;
		// the root is the first/only node on level 0
		int level_N = 1;
  
		while (! Q_empty (Q, N, qh, qt))
		{
			v = Q_deq (Q, N, &qh);
			i = v.id;
			c++;
  
			// for computing the inverse permutation P where
			// A(inv(P),inv(P)) or P'*A*P is banded
			// P(i) = c;
  
			// for computing permutation P where
			// A(P(i),P(j)) or P*A*P' is banded
			P[c]= i;
  
			// put all unvisited neighbors j of node i on the heap
			s = 0;
			int j1 = cidx[i];
			int j2 = cidx2[i];
			while (j1 < cidx[i+1] || j2 < cidx2[i+1])
			{
				if (j1 == cidx[i+1])
				{
					int r2 = ridx2[j2++];
					if (! visit[r2])
					{
						// the distance of node j is dist(i)+1
						w.id = r2;
						w.deg = D[r2];
						w.dist = v.dist+1;
						H_insert (S, &s, &w);
						visit[r2] = 1; //true
					}
				}
				else if (j2 == cidx2[i+1])
				{
					int r1 = ridx[j1++];
					if (! visit[r1])
					{
						w.id = r1;
						w.deg = D[r1];
						w.dist = v.dist+1;
						H_insert (S, &s, &w);
						visit[r1] = 1; //true
					}
				}
				else
				{
					int r1 = ridx[j1];
					int r2 = ridx2[j2];
					if (r1 <= r2)
					{
						if (! visit[r1])
						{
							w.id = r1;
							w.deg = D[r1];
							w.dist = v.dist+1;
							H_insert (S, &s, &w);
							visit[r1] = 1; //true
						}
						j1++;
				                if (r1 == r2)
							j2++;
					}
					else
					{
						if (! visit[r2])
						{
							w.id = r2;
							w.deg = D[r2];
							w.dist = v.dist+1;
							H_insert (S, &s, &w);
							visit[r2] = 1; //true
						}
						j2++;
					}
				}
			}
  
			// add the neighbors to the queue (sorted by node degree)
			while (! H_empty (S, s))
			{
				// locate a neighbor of i with minimal degree in O(log(N))
				v = H_remove_min (S, &s, 1);
  
				// entered the BFS a new level?
				if (v.dist > level)
				{
					// adjustment of bandwith:
					// "[...] the minimum bandwidth that
					// can be obtained [...] is the
					//  maximum number of nodes per level"
					if (Bsub < level_N)
						Bsub = level_N;
  
					level = v.dist;
					// v is the first node on the new level
					level_N = 1;
				}
				else
				{
					// there is no new level but another node on
					// this level:
					level_N++;
				}
  
				// enqueue v in O(1)
				Q_enq (Q, N, &qt, &v);
			}
  
			// synchronize the bandwidth with level_N once again:
			if (Bsub < level_N)
				Bsub = level_N;
		}
		// finish of BFS. If there are still unvisited nodes in the graph
		// then it is split into CCs. The computed bandwidth is the maximum
		// of all subgraphs. Update:
		if (Bsub > B)
			B = Bsub;
	}
	// are there any nodes left?
	while (c+1 < N);
 
	free(Q);
	free(S);
	free(D);
//	free(cidx);
//	free(ridx);
	free(cidx2);
	free(ridx2);
	free(visit);
 
	// compute the reverse-ordering
	s = N / 2 - 1;
	int i, j;
	for (i = 0, j = N - 1; i <= s; i++, j--)
		SWAP(P[i], P[j]);
  
}

void transform_CSR_to_CCS(int N,const int *ja, const int *ia, int *cidx,int *ridx)
{
	int i, j, count, nz = ia[N];
	ARRAY *IJ;

	IJ = calloc(nz, sizeof(ARRAY));
	
	count = 0;
	for (i=0; i<N; i++){
		for (j=ia[i]; j<ia[i+1]; j++){
			IJ[count].arr1 = i;
			IJ[count].arr2 = ja[j];
			count++;
		}
	}

	qsort(IJ,nz,sizeof(ARRAY),COMPARE_array);

	for (i=0; i<nz;i++){
		cidx[IJ[i].arr2]++;
		ridx[i] = IJ[i].arr1;
	}

	for (i=N-1;i>=0;i--)
		cidx[i+1] = cidx[i];

	cidx[0]=0;		
	for (i=0; i<N; i++)
		cidx[i+1] += cidx[i];

	free(IJ);
}

void transpose (int N, const int *ridx, const int *cidx, int *ridx2, int *cidx2)
{
	int i, j, k, q, nz = cidx[N];
	int *w;
	
	w = calloc(N+1, sizeof(int));
  
	for (i = 0; i < N; i++)
		w[i] = 0;
	for (i = 0; i < nz; i++)
		w[ridx[i]]++;
	nz = 0;
	for (i = 0; i < N; i++)
	{
		cidx2[i] = nz;
		nz += w[i];
		w[i] = cidx2[i];
	}
	cidx2[N] = nz;
	w[N] = nz;
  
	for (j = 0; j < N; j++)
	{
		for (k = cidx[j]; k < cidx[j + 1]; k++)
		{
			q = w[ridx[k]]++;
			ridx2[q] = j;
		}
	}

	free(w);
}

int calc_degrees(int N, const int *ridx, const int *cidx, int *D)
{
	int i, j, k, l, max_deg = 0;
  
	for (i = 0; i < N; i++)
		D[i] = 0;
  
	for (j = 0; j < N; j++)
	{
		for (i = cidx[j]; i < cidx[j+1]; i++)
          	{
            		k = ridx[i];

			// there is a nonzero element (k,j)
			D[k]++;
			if (D[k] > max_deg)
				max_deg = D[k];

			// if there is no element (j,k) there is one in
			// the symmetric matrix:
			if (k != j)
			{
				int found = 0;
				for (l = cidx[k]; l < cidx[k + 1]; l++)
                  		{
					if (ridx[l] == j)
                      			{
			                        found = 1;
						break;
					}
					else if (ridx[l] > j)
						break;
				}	
  
				if (! found)
				{
					// A(j,k) == 0
					D[j]++;
					if (D[j] > max_deg)
						max_deg = D[j];
                  		}
			}
		}
	}
	return max_deg;
}

int find_starting_node (int N, const int *ridx, const int *cidx, const int *ridx2, const int *cidx2, int *D, int start)
{
	CMK_NodeType w;
	CMK_NodeType *Q = calloc(N+1, sizeof(CMK_NodeType));
	
	int *visit = calloc(N, sizeof(int)); // visit[0..N-1] receive zero. It means visit_i is false initially
  
	int qh = 0;
	int qt = 0;

	CMK_NodeType x;
	x.id = start;
	x.deg = D[start];
	x.dist = 0;
	Q_enq (Q, N, &qt, &x);
	visit[start] = 1; //true
  
	// distance level
	int level = 0;
	// current largest "eccentricity"
	int max_dist = 0;
  
	for (;;)
	{
		while (!Q_empty (Q, N, qh, qt))
		{
			CMK_NodeType v = Q_deq (Q, N, &qh);
  
			if (v.dist > x.dist || (v.id != x.id && v.deg > x.deg))
				x = v;
  
	           	int i = v.id;
  
		        // add all unvisited neighbors to the queue
			int j1 = cidx[i];
			int j2 = cidx2[i];
			while (j1 < cidx[i+1] || j2 < cidx2[i+1])
			{
				if (j1 == cidx[i+1])
				{
					int r2 = ridx2[j2++];
					if (! visit[r2])
					{
						// the distance of node j is dist(i)+1
						w.id = r2;
						w.deg = D[r2];
						w.dist = v.dist+1;
						Q_enq (Q, N, &qt, &w);
						visit[r2] = 1; //true
  
						if (w.dist > level)
							level = w.dist;
					}
				}
				else if (j2 == cidx2[i+1])
                  		{
					int r1 = ridx[j1++];
					if (! visit[r1])
					{
						// the distance of node j is dist(i)+1
						w.id = r1;
						w.deg = D[r1];
						w.dist = v.dist+1;
						Q_enq (Q, N, &qt, &w);
						visit[r1] = 1; //true
 
						if (w.dist > level)
							level = w.dist;
					}
				}
				else
				{
					int r1 = ridx[j1];
					int r2 = ridx2[j2];
					if (r1 <= r2)
					{
						if (! visit[r1])
						{
							w.id = r1;
							w.deg = D[r1];
							w.dist = v.dist+1;
							Q_enq (Q, N, &qt, &w);
							visit[r1] = 1;//true

							if (w.dist > level)
								level = w.dist;
                          			}
						j1++;
						if (r1 == r2)
							j2++;
					}
                    			else
                      			{
						if (! visit[r2])
						{
							w.id = r2;
							w.deg = D[r2];
							w.dist = v.dist+1;
							Q_enq (Q, N, &qt, &w);
							visit[r2] = 1; //true
  
                            				if (w.dist > level)
                              					level = w.dist;
                          			}
						j2++;
					}
				}
			}
		} // finish of BFS
  
		if (max_dist < x.dist)
		{
			max_dist = x.dist;

			int i;	
  			for (i = 0; i < N; i++)
				visit[i] = 0; //false

			visit[x.id] = 1; //true
			x.dist = 0;
			qt = qh = 0;
			Q_enq (Q, N, &qt, &x);
		}
		else
			break;
	}

	free(Q);
	free(visit);

	return x.id;
}

// A simple queue.
// Queues Q have a fixed maximum size N (rows,cols of the matrix) and are
// stored in an array. qh and qt point to queue head and tail.
   
// Enqueue operation (adds a node "o" at the tail)
void Q_enq (CMK_NodeType *Q, int N, int *qt, const CMK_NodeType *o)
{
	Q[*qt] = *o;
	*qt = (*qt + 1) % (N + 1);
}

// Dequeue operation (removes a node from the head)
CMK_NodeType Q_deq (CMK_NodeType *Q, int  N, int *qh)
{
	CMK_NodeType r = Q[*qh];
     	*qh = (*qh + 1) % (N + 1);
	return r;
}

// Heap operation insert. Running time is O(log(n))
void H_insert (CMK_NodeType *H, int *h, const CMK_NodeType *o)
{
	int i = (*h)++;
  
	H[i] = *o;

	if (i == 0)
		return;

	do
	{
		int p = PARENT(i);
		if (H[i].deg < H[p].deg)
		{
			SWAP(H[i], H[p]);
			i = p;
		}
		else
			break;
	}
	while (i > 0);
}
 
// Heap operation remove-min. Removes the smalles element in O(1) and
// reorganizes the heap optionally in O(log(n))
CMK_NodeType H_remove_min (CMK_NodeType *H, int *h, int reorg/*=1*/)
{
	CMK_NodeType r = H[0];
	H[0] = H[--(*h)];
	if (reorg)
		H_heapify_min (H, 0, *h);

	return r;
}

// Builds a min-heap (the root contains the smallest element). A is an array
// with the graph's nodes, i is a starting position, size is the length of A.
void H_heapify_min (CMK_NodeType *A, int i, int size)
{
	int j = i;
	for (;;)
	{
		int l = LEFT(j);
		int r = RIGHT(j);
  
		int smallest;
		if (l < size && A[l].deg < A[j].deg)
			smallest = l;
		else
			smallest = j;

		if (r < size && A[r].deg < A[smallest].deg)
			smallest = r;
  
		if (smallest != j)
		{
			SWAP(A[j], A[smallest]);
			j = smallest;
		}
		else
			break;
	}
}


