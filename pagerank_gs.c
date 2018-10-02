#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Source: https://www.ece.ucsb.edu/~hespanha/published/2018ACC_0753_MS.pdf

double tolerance(double**,double*,double,int);

void main(int argc, char *argv[])
{
	struct timeval startwtime, endwtime;
	double seq_time;
	
	int i, j;	
	
	double damp = 0.85; //damping factor
	
	/*Input arguments*/
	int n;			 //matrix dimension
	double conv_tol; //convergance tolerance
	int max_iter;	 //maximum number of iterations to be made
	char *path; 	//path to file or filename if in the same folder
	
	double tol = 10; // initialized with a value >> conv_tol
	int iterations = 0; //actual iterations
	
	double b = (1 - damp) / n;
		
	if (argc != 5)
	{
		printf( "Wrong number of arguments.\n1.Dimension\n2.Convergance tolerance\n3.Maximum iterations\n4.Path to adjacency list\n");
		exit( 1 );
	}
	
	n = atoi(argv[1]);
	conv_tol = atof(argv[2]);
	max_iter = atoi(argv[3]);
	path = argv[4];
	
	int *count_ones = (int*) calloc(n,sizeof(int));
		
	//Adjacency matrix memory allocation
	double **B = (double**) calloc(n,sizeof(double*));
    for(i=0; i<n; i++)
		B[i] = (double*) calloc(n,sizeof(double));
	
	
	/****************************************************************************************************/
	/* Read adjacency matrix from file. Matrix in the form of: "Node-i: Node-1 Node-2 ... Node-n -1" 
		*That means Node_num has an outgoing link to Node-1 Node-2 ... Node-n
		*Each node has a line.
		*If node is a dangling node then: "Node-i: -1"
	/****************************************************************************************************/
	FILE *flist;
	char  list_file[1000];	
	sprintf(list_file,"%s",path);
	flist = fopen(list_file,"r");
	for(i=0; i<n; i++)
	{
		fscanf(flist,"%*d: %d",&j); // ignore #: at start of line e.g. "12: 13" will be read as "13"
		while (j != -1)
		{
			B[i][j] = 1;
			fscanf(flist,"%d",&j);
		}
	}
	fclose(flist);
		
	
	//Counts outgoing nodes for each node. If the node has no outgoing links (dangling node), assume that it connects to every other node and change accordingly B and count_ones.
	for(j=0; j<n; j++)
	{
			
		for(i=0; i<n; i++)
			if(B[i][j] == 1) count_ones[j]++;
		
		if (!count_ones[j])
		{
			for(i=0; i<n; i++)
				B[i][j] = 1;
			count_ones[j] = n;
		}
	}

	
	//Normalized propability matrix
	double **S = (double**) calloc(n,sizeof(double*));
	for(i=0; i<n; i++)
		S[i] = (double*) calloc(n,sizeof(double));
	
	//S = (B * damp) / outgoing_links
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			if(count_ones[j]) S[i][j] = (B[i][j] * damp) / count_ones[j];
				
	//I-S
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
		{
			if(i == j) S[i][j] = 1.0 - S[i][j];
			else S[i][j] = 0.0 - S[i][j];
		}
		
	
	double *pagerank = (double*) malloc(n*sizeof(double));
	double *pagerank_old = (double*) malloc(n*sizeof(double));
	
	//Initialize each node pagerank with 1/n
	for(i=0; i<n; i++)
		pagerank[i] = (1/(double)n);
	
	for(i=0; i<n; i++)
		pagerank_old[i] = pagerank[i];
	
	double sum = 0 , sum_new = 0 ;
	
	gettimeofday( &startwtime, NULL );// start timer

	/**************************************/
	/***	 Gauss Seidel itertaions 	***/
	/**************************************/
	while(iterations < max_iter && tol > conv_tol)
	{
		iterations++;
					
		for(i=0; i<n; i++)
		{
			for(j=0; j<n; j++)
			{
					if(j<i)
						sum_new += S[i][j] * pagerank[j];

					if(j>i)
						sum += S[i][j] * pagerank_old[j];
			}

			pagerank[i] = 1/S[i][i] * ( b - sum_new - sum);

			sum = 0.0;
			sum_new = 0.0;
		}
				
		tol = tolerance(S, pagerank, b, n);
		
		for(i=0; i<n; i++)
			pagerank_old[i] = pagerank[i];
	}
	
	gettimeofday( &endwtime, NULL );// end timer
	seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);

	printf("Time: %f\n", seq_time);
	printf("Iterations: %d\n", iterations);
	printf("Error: %f\n", tol);
}

//Calculates the Euclidean norm of Ax-b (convergance condition)
double tolerance(double **A, double *x, double b, int dim)
{
  double sum, tol = 0;
  for(int i=0; i<dim; i++)
  {
    sum = 0;
    for(int j=0;j<dim;j++)
      sum += A[i][j]*x[j];
    sum -= b;

    tol += sum*sum;
  }

  return sqrt(tol);
}
