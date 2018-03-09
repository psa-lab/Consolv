/************************************************/
/*   Consolv was developed by Michael Raymer,   */
/*   Paul Sanschagrin, and Leslie Kuhn at the   */
/*   Protein Structural Analysis and Design     */
/*   Laboratory,  Department of Biochemistry    */
/*   (517) 353-8745    kuhn@agua.bch.msu.edu    */
/*          www.bch.msu.edu/labs/kuhn           */
/*                                              */
/*   in collaboration with William Punch and    */
/*   Erik Goodman of the Genetic Algorithms     */
/*   Research and Applications Group, at        */
/*   Michigan State University.                 */
/*                                              */
/*         (c) 1998 Board of Trustees           */
/*          Michigan State University           */
/*                                              */
/*   THIS IS FREE SOFTWARE, AND IS DISTRIBUTED  */
/*   ACCORDING TO THE TERMS OF THE ACCOMPANYING */
/*   LICENSE AGREEMENT. INFORMATION ABOUT ITS   */
/*   COPYING AND USE IS CONTAINED IN THE FILE   */
/*   'license.txt' WHICH ACCOMPANIES THIS       */
/*   DISTRIBUTION.  THERE IS NO WARRANTY,       */
/*   EXPRESSED OR IMPLIED, ASSOCIATED WITH      */
/*   THIS SOFTWARE, AND ALL RISK FOR ITS USE    */
/*   RESTS WITH THE USER.                       */
/************************************************/

/***********************************************************************/
/* bbknn.c - Branch and bound K-NN pattern classification algorithm.   */
/*                                                                     */
/* Reference:  This code is an implementation of the algorithm         */
/*             described in                                            */
/*               K. Fukunaga & P. M. Narendra, "A Branch and Bound     */
/*               Algorithm for Computing k-Nearest Neighbors,"  IEEE   */
/*               Transactions on Computers, v.24, pp.750-753,          */
/*               July, 1975.                                           */
/*                                                                     */
/*  Initial version: Michael Raymer, 3/29/96                           */
/*                    raymermi@cps.msu.edu                             */
/*                                                                     */
/*  Some code segments (specifically, parts of the  read_patterns      */
/*  subroutine) were taken from a knn classifier written by            */
/*     P.J. Flynn and modified by C. S. Lai,                           */
/*     Pattern Recognition and Image Processing Laboratory,            */
/*     Michigan State University.  This knn contained the following    */
/*     reference:                                                      */ 
/*      Reference: "Introduction to Statistical Pattern Recognition",  */
/*                  P.A. Devijver and J. Kittler.                      */
/*                                                                     */
/* 12/96 - Version numbering started, initial numbered version is      */
/*         4.2 to match the bbknnf version numbering.                  */
/*                                                                     */
/* Version history:                                                    */
/*   v4.2 - Added CBV_SCALING.  Controlled by the a #DEFINE in the     */
/*          defines.h file.  Added version string printing in output.  */
/*          M. Raymer, 12/17/96.                                       */
/***********************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include "defines.h"
#include "protos.h"

#define BBKNN_VERSION "bbknnf.c, version 4.2"

/***********************************************************************/
/* SUBROUTINES AND FUNCTIONS                                           */
/***********************************************************************/

/***********************************************************************/
/* read_patterns - read a set of patterns and labels from a file       */
/*                                                                     */
/*   This routine will read a set of patterns from an input file,      */
/*   allocating space for the patterns as necessary.  It returns       */
/*   two arrays:                                                       */
/*      patterns - a 2 dimensional (n x d) array of doubles containing */
/*                 the patterns read.                                  */
/*      labels - a 1 dimensional array of integers containing the      */
/*                 class labels of the patterns.                       */
/*   The routine also returns:                                         */
/*      num_patterns - the number of patterns read                     */
/*      num_features - the number of features for each pattern         */
/*                                                                     */
/*   Inputs include:                                                   */
/*      fptr - a pointer to the file to read from                      */
/*                                                                     */
/*   Input file format:                                                */
/*    - Lines starting with a '#' are ignored.                         */
/*    - The first non-blank, non-comment line should be the number     */
/*      of features.  This line should be of the form:                 */
/*       Features: 4                                                   */
/*    - The rest of the lines should simply contain feature values,    */
/*      separated by white space, followed by an integer label.        */
/*                                                                     */
/*    So, for a 3-feature set of patterns from two classes, the        */
/*    input file might look like this:                                 */
/*                                                                     */
/*      # Sample input file: 3 features, 4 patterns, 2 classes         */
/*      Features: 3                                                    */
/*      12.2   1.7  3.5  1                                             */
/*      0.8    2.5  1.0  1                                             */
/*      8.8   -2.0  0.9  2                                             */
/*      1.2    2.2  1.3  2                                             */   
/*                                                                     */
/*   The return value from the routine is either SUCCESS or            */
/*   FAILURE.                                                          */
/***********************************************************************/

int read_patterns(fptr, patterns, labels, num_patterns, num_features,
		  labelled)
     FILE *fptr;		/* File to read from (INPUT)             */
     double **patterns;		/* The patterns read (OUTPUT)            */
     int **labels;		/* The labels read (OUTPUT)              */
     int *num_patterns;		/* The number of patterns read (OUTPUT)  */
     int *num_features;		/* The number of features read (OUTPUT)  */
     int labelled;		/* FLAG:  Is the input labelled? (INPUT) */
{
  char inbuf[MAXLINE+1];	/* Buffer for input lines */
  char jbuf[MAXLINE+1];		/* Junk buffer */
  char token[MAXLINE+1];	/* Buffer for parsing tokens */
  int matches;			/* Number of matches during parsing */
  int nmax;			/* The max #patterns we have room for so far */
  int nread;			/* The number of patterns read so far */
  int abort;			/* A read-abort flag */
  int i;			/* Loop counter */

  /* Skip opening comments and blank lines */
  inbuf[0] = '\0';
  while (inbuf[0] == '#' || isblank(inbuf)) {
    if (feof(fptr)) {
      fprintf(stderr, "Error reading patterns: File empty.\n");
      return FAILURE;
    }
    else
      fgets(inbuf, MAXLINE, fptr);
  }

  /* The first non-comment, non-blank line should contain the number */
  /* of patterns.                                                    */
  matches = sscanf(inbuf, "%s %d", token, num_features);
  upcase(token);
  
  /* If we found some other string, print an error */
  if (strncmp(token, "FEATURES", 8)) {
    fprintf(stderr, "Error reading patterns: Could not find FEATURES tag.\n");
    fprintf(stderr, "Found '%s' instead.\n", token);
    return FAILURE;
  }

  /* If we can't parse the number of features, also print */
  /* an error.                                            */
  else if (matches < 2) {
    fprintf(stderr,
	    "Error reading patterns: could not parse number of features.\n");
    return FAILURE;
  }

  /* Debugging printout */
  if (DEBUG_LEVEL >= 1)
    printf ("Identified %d features.\n", *num_features);

  /* Allocate room for MBLOCK patterns & labels */
  nmax = MBLOCK;
  *patterns = (double *) malloc (MBLOCK * (*num_features) * sizeof(double));
  if (labelled)
    *labels = (int *) malloc (MBLOCK * sizeof(int));
  else
    *labels = (int *) NULL;

  /* If we can't allocate enough memory, print an error */
  if (!(*patterns)) {
    fprintf(stderr, "Error reading patterns: insufficient memory.\n");
    return FAILURE;
  }
  if (labelled && !(*labels)) {
    fprintf(stderr, "Error reading labels: insufficient memory.\n");
    return FAILURE;
  }

  /* Read the patterns, enlarging the pattern and label arrays as necessary */
  nread = 0;
  abort = FALSE;
  while (!feof(fptr)) {
  
    /* Check if the arrays must be enlarged... */
    if (nread == nmax) {
      nmax += MBLOCK;
      /* Enlarge and check for memory errors */
      *patterns = (double *)
	realloc (*patterns, nmax * (*num_features) * sizeof(double));
      if (labelled)
	*labels = (int *) realloc (*labels, nmax * sizeof(int));
      if (!(*patterns)) {
	fprintf(stderr, "Error reading patterns: insufficient memory.\n");
	return FAILURE;
      }
      if (!(*labels) && labelled) {
	fprintf(stderr, "Error reading labels: insufficient memory.\n");
	return FAILURE;
      }
    }

    /* Read the current pattern's features */
    for (i=0; i < *num_features; i++)
      if (fscanf(fptr, "%lf", (*patterns)+(((*num_features)*nread) + i)) < 1) {
	if (i > 0) {
	  fprintf(stderr,"Error reading patterns: Pattern %d is incomplete.\n",
		  nread+1);
	  return FAILURE;
	}
	else {
	  abort = TRUE;
	  break;
	}
      }

    /* Read the current pattern's label */
    if (labelled && !abort)
      if (fscanf(fptr, "%d", (*labels)+nread) < 1) {
	fprintf(stderr, "Error reading patterns: Pattern %d has no label.\n",
		nread+1);
	return FAILURE;
      }

    /* Clear to end of line, if necessary */
    if (!feof(fptr))
      fgets(jbuf, MAXLINE, fptr);

    /* Go to the next pattern */
    if (!abort)
      nread++;
  }

  /* Shrink the arrays to the proper size */
  *patterns = (double *) realloc (*patterns,
				  nread * (*num_features) * sizeof(double));
  if (labelled)
    *labels = (int *) realloc (*labels, nread * sizeof(int));
  
  *num_patterns = nread;

  return SUCCESS;
}

/***********************************************************************/
/* sqdist() - compute the squared euclidean distance between two       */
/*            d-dimensional patterns.                                  */
/*                                                                     */
/*  This routine takes two d-dimensional vectors of doubles and finds  */
/*  the squared euclidean distance between them.                       */
/***********************************************************************/
double sqdist(pat1, pat2, d)
     double *pat1;		/* First pattern */
     double *pat2;		/* Second pattern */
     int d;			/* The dimensionality of the patterns */
{
  double curdist;		/* The squared distance */
  double dst;			/* Distance along each dimension */
  int i;			/* A loop counter */
 
  curdist = 0.0;
  for (i=0; i<d; i++) {
    dst = pat1[i] - pat2[i];
    curdist += dst * dst;
  }

  return curdist;
}

/***********************************************************************/
/* compute_mean() - compute the multivariate mean of a set of samples. */
/*                                                                     */
/* This routine takes in a node, which must have valid values for      */
/* samples and num_samples, and computes the sample mean vector for    */
/* the node.                                                           */
/***********************************************************************/
void compute_mean(root, patterns, d)
     treenode *root;		/* The node to compute the sample mean for */
     double *patterns;		/* The training patterns */
     int d;			/* The number of dimensions */
{
  int i, j;			/* Loop counters */

  /* Delete the old mean vector, if there is one */
  if (root->sample_mean) {
    free(root->sample_mean);
    root->sample_mean = (double *) NULL;
  }

  /* Allocate a new mean vector */
  root->sample_mean = (double *) malloc(d * sizeof(double));
  for (i=0; i < d; i++)
    (root->sample_mean)[i] = 0.0;

  /* Compute the values:                                         */
  /* For each element, i, of the samples array, patterns[i]      */
  /* contains the d values for that sample.  Each value is added */
  /* to the sample_mean vector, which is finally divided by      */
  /* num_samples.                                                */
  for (i=0; i < root->num_samples; i++)
    for (j=0; j < d; j++)
      (root->sample_mean)[j] += patterns[((root->samples)[i] * d) + j];
  for (i=0; i < d; i++)
    (root->sample_mean)[i] /= root->num_samples;

}

/***********************************************************************/
/* compute_radius - find the radius of a cluster                       */
/*                                                                     */
/* This routine takes a node, which must have valid values for the     */
/* the samples array, num_samples, and a valid mean vector, and finds  */
/* the radius of the node, which is the maximum distance from the      */
/* mean to any point in the cluster.                                   */
/***********************************************************************/
void compute_radius(root, patterns, d)
     treenode *root;		/* The node to compute the radius of */
     double *patterns;		/* The patterns array */
     int d;			/* The dimensionality of the data */
{
  double curdist;		/* Current distance */
  int i;			/* loop counter */

  root->radius = 0.0;

  /* Find the largest squared-euclidean distance from the mean */
  for (i=0; i < root->num_samples; i++) {
    curdist = sqdist(&(patterns[(root->samples)[i]*d]), root->sample_mean, d);
    if (curdist > root->radius)
      root->radius = curdist;
  }

  /* Store the non-squared distance */
  root->radius = sqrt(root->radius);
}

/***********************************************************************/
/* dump_cluster - dump the contents of a cluster.  Used only for       */
/*                debugging purposes.                                  */
/*                                                                     */
/*   root - the root node of the subtree to dump.                      */
/*   patterns - the patterns array.                                    */
/*   d - the dimensionality (# features) in the patterns.              */
/***********************************************************************/
void dump_cluster(root, patterns, d)
     treenode *root;
     double *patterns;
     int d;
{
  /* Local variables */
  int i, j;			/* loop counters */

  printf("Dumping cluster\n");
  printf("---------------\n");
  printf("Samples: %d\tRadius: %lf\n", root->num_samples, root->radius);
  printf("Mean: ");
  for(i=0; i<d; i++)
    printf("%lf\t", root->sample_mean[i]);
  printf("\nSamples:\n");
  for(i=0; i < root->num_samples; i++) {
    printf("%5d: ",i);
    for(j=0; j<d; j++)
      printf("%lf  ", patterns[(root->samples[i]*d)+j]);
    printf("\n");
  }
  printf("\n");
}

/***********************************************************************/
/* three_means() - Find a 3 means clustering of the input data         */
/*                                                                     */
/* The input to this routine is a tree node.  The first time the       */
/* routine is called, this node will contain every member of the       */
/* training set.                                                       */
/*                                                                     */
/* The routine generates 3 new tree nodes, which are assigned to the   */
/* left, right, and middle children of the input node.  The members    */
/* of each node are found by doing a 3-means clustering of the data    */
/* in the input node.                                                  */
/*                                                                     */
/* The method of selection of the 3 initial cluster means is specified */
/* by the 'method' parameter.  The UNIFORM method attempts to find a   */
/* more uniform clustering by selecting:                               */
/*  -the pattern closest to the origin,                                */
/*  -the pattern farthest from the origin, and                         */
/*  -the centroid of the entire cluster                                */
/* as the initial cluster means.  Although this can produce more even- */
/* ly distributed clusters, it requires n distance calculations and    */
/* (n*d) additions.  The FAST method simply selects the first, last    */
/* and middle (n/2) pattern as the initial cluster centers.            */
/***********************************************************************/
void three_means(root, patterns, d, dist_to_mean, method)
     treenode *root;		/* The node to compute three children for */
     double *patterns;		/* The pattern array */
     int d;			/* The dimensionality of the data */
     double *dist_to_mean;	/* Distance from each pattern to the mean   */
				/* of the cluster it's assigned to (OUTPUT) */
     int method;		/* FLAG: how to choose the initial centers */
{
  int i, j;			/* Loop counters */
  int cur_sample;		/* The current sample (loop index) */
  int iterations;		/* Iteration counter for 3-means loop */
  int sample_moved;		/* Flag: did a sample change clusters? */
  double curdist;		/* The current distance */
  double mindist;		/* The minimum distance so far from a    */
				/* pattern to the origin or cluster mean */
  double maxdist;		/* Maximum "   "   "    "                */
  double dist[3];		/* The distance from a pattern to each   */
				/* cluster mean */
  int *cluster;			/* The cluster number of each sample */
  int num_members[3];		/* The number of members of each cluster */
  int min;			/* The closest cluster id */
  int first, middle, last;	/* The first/mid/last sample for the node */
  int n;			/* Shorthand for root->num_samples */
  double radius[3];		/* The radius of each cluster */

  static double *origin = NULL;	/* The origin */

  /* The 3 d-dimensional cluster center vectors */
  static double *center[3] = {NULL, NULL, NULL};

  /* If the cluster center vectors have not been allocated, do so. */
  if (!center[0])
    for (i=0; i<3; i++)
      center[i] = (double *) malloc (d * sizeof(double));

  /* Same for the origin (allocate if necessary) */
  if (!origin) {
    origin = (double *) malloc (d * sizeof(double));
    for (i=0; i<d; i++)
      origin[i] = 0;
  }

  /**********************************************************************/
  /* Compute starting cluster centers:                                  */
  /*                                                                    */
  /* If the 'method' is UNIFORM, then the cluster centers are chosen as */
  /* follows:                                                           */
  /* center[0] = point closest to origin                                */
  /* center[1] = cluster mean                                           */
  /* center[2] = point farthest from origin                             */
  /**********************************************************************/

  if (method == UNIFORM) {

    /* Initialize the minimum and maximum clusters */
    mindist = -1.0;
    maxdist = -1.0;

    /* Initialize all cluster centers */
    for (i=0; i<3; i++)
      for(j=0; j<d; j++)
	center[i][j] = 0.0;

    /* Find the minimum and maximum */
    for (i=0; i<root->num_samples; i++) {
      curdist = sqdist(&(patterns[(root->samples)[i]*d]), origin, d);
      if ((curdist < mindist) || (mindist == -1.0)) {
	mindist = curdist;
	for (j=0; j<d; j++)
	  center[0][j] = patterns[((root->samples)[i]*d)+j];
      }
      if ((curdist > maxdist) || (maxdist == -1.0)) {
	maxdist = curdist;
	for (j=0; j<d; j++)
	  center[2][j] = patterns[((root->samples)[i]*d)+j];
      }
    }

    /* Set the middle cluster center */
    for (i=0; i<d; i++)
      center[1][i] = (root->sample_mean)[i];
  }

  /*******************************************************************/
  /* If method is FAST, then simply select the 1st, n-th and (n/2)th */
  /* patterns as the initial cluster centers.                        */
  /*******************************************************************/
  else if (method == FAST) {
    /* Simply choose the 1st, nth, and (n/2)th patterns */
    n = root->num_samples;

    if (n > 0) {
      first = root->samples[0];
      last = root->samples[n-1];
      middle = root->samples[n/2];
      for (i=0; i<d; i++) {
	center[0][i] = patterns[first*d+i];
	center[1][i] = patterns[middle*d+i];
	center[2][i] = patterns[last*d+i];
      }
    }
  }

  /* Error check */
  else {
    fprintf(stderr, "Unknown method for 3-means computation. Exiting.\n");
    exit (-1);
  }

  /* Debugging printout */
  if (DEBUG_LEVEL >= 3) {
    printf("3-Means computation:  ");
    if (method == UNIFORM) printf ("(method = uniform)\n");
    else printf("(method = FAST)\n");
    printf("Center0\t\tCenter1\t\tCenter2\n");
    for (i=0; i<d; i++)
      printf("%lf\t%lf\t%lf\n", center[0][i], center[2][i], center[1][i]);
    printf("\n");
  }

  /****************************************************************/
  /* Now that we have the starting cluster centers, run a 3-means */
  /* until no pattern changes clusters, or until ITERATION_LIMIT  */
  /* is reached.                                                  */
  /****************************************************************/

  /* Allocate an array of members for each cluster */
  cluster = (int *) malloc((root->num_samples) * sizeof(int));
  
  iterations = 0;
  sample_moved = TRUE;
  while ((iterations < ITERATION_LIMIT) && sample_moved) {

    /* in loop initializations */
    if (iterations > 0) sample_moved = FALSE;
    iterations++;
    for(i=0; i<3; i++)
      num_members[i] = 0;

    /* Step 1:  Assign each pattern to the nearest cluster center */
    /* ---------------------------------------------------------- */
    for (cur_sample = 0; cur_sample < root->num_samples; cur_sample++) {

      /* find the distance to each cluster center */
      for (i=0; i < 3; i++)
	dist[i] = sqdist(center[i],
			 &(patterns[(root->samples)[cur_sample]*d]), d);

      /* Assign the sample to the closest one, note whether it moves or not */
      min=0;
      for (i=1; i < 3; i++)
	if (dist[i] < dist[min])
	  min = i;
      
      if (min != cluster[cur_sample]) {
	cluster[cur_sample] = min;
	sample_moved = TRUE;
      }

      /* Keep track of the distance to the mean of the cluster each */
      /* pattern gets assigned to...                                */
      dist_to_mean[root->samples[cur_sample]] = dist[min];

      num_members[min]++;

    }

    /* Step 2:  Re-calculate cluster means */
    /* ----------------------------------- */

    /* initialize to zero */
    for (i=0; i<3; i++)
      for (j=0; j<d; j++)
	center[i][j] = 0.0;

    /* sum up the distance along each dimension */
    for (cur_sample = 0; cur_sample < root->num_samples; cur_sample++) {
      for (i=0; i<d; i++)
	center[cluster[cur_sample]][i] += 
	  patterns[((root->samples)[cur_sample]*d)+i];
    }

    /* divide by the number of members in the cluster */
    for (i=0; i<3; i++)
      for (j=0; j<d; j++)
	if (num_members[i] > 0)
	  center[i][j] /= num_members[i];
    
  }

  /* After the 3-means has converged or stopped, find the radii */
  /* of the resulting clusters...                               */
  /* So far we have been working with squared-euclidean-distance  */
  /* as we find the radii, convert the radii and the dist-to-mean */
  /* array to non-squared euclidean distance.                     */

  for (i=0; i<3; i++) {
    radius[i] = 0.0;

    for( j=0; j < (root->num_samples); j++)
      if (cluster[j] == i) {
	dist_to_mean[root->samples[j]] = sqrt(dist_to_mean[root->samples[j]]);
	curdist = dist_to_mean[root->samples[j]];
	if (curdist > radius[i])
	  radius[i] = curdist;
      }
  }
	  

  /* Debug printout */
  if (DEBUG_LEVEL >= 2) {
    printf("3-Means converged/halted after %d iterations\n", iterations);
    printf("(Iteration limit:  %d)\n", ITERATION_LIMIT);
    printf("Members: 0=%d, 1=%d, 2=%d\n\n", num_members[0], num_members[1],
	   num_members[2]);
  }
  if (DEBUG_LEVEL >= 3) {
    printf("Cluster0\tCluster1\tCluster2\n");
    for (i=0; i<d; i++)
      printf("%lf\t%lf\t%lf\n", center[0][i], center[1][i], center[2][i]);
    printf("\n");
  }

  /****************************************************************/
  /* Last step:                                                   */
  /* Build three child nodes, and link them from the current node */
  /****************************************************************/
  for (i=0; i<3; i++) {

    /* Allocate the child and fill in the fields: */
    root->child[i] = (treenode *) malloc (sizeof(treenode));
    (root->child[i])->radius = radius[i];

    /* The children of the child are NULL */
    for (j=0; j<3; j++)
      (root->child[i])->child[j] = (treenode *) NULL;

    /* Allocate an array of d doubles for the mean vector and fill it in */
    (root->child[i])->sample_mean = (double *) malloc (d * sizeof(double));
    for (j=0; j<d; j++)
      (root->child[i])->sample_mean[j] = center[i][j];

    /* Allocate an array of integers for the members of the child */
    /* and fill them in...                                        */
    (root->child[i])->samples = (int *) malloc (num_members[i] * sizeof(int));
    (root->child[i])->num_samples = 0;

    /* for each member of the root node */
    for (j=0; j < root->num_samples; j++) {

      /* If it should be assigned to this child, add it to this child's */
      /* member array, and increment the num_samples for this child.    */
      if (cluster[j] == i) {
	((root->child[i])->samples)[(root->child[i])->num_samples] =
	  (root->samples)[j];
	((root->child[i])->num_samples)++;
      }
    }

    /* Check our work by making sure the number of samples assigned */
    /* to this child matches the num_members for this cluster...    */
    if ((root->child[i])->num_samples != num_members[i]) {
      fprintf(stderr, "Internal fatal error:\n");
      fprintf(stderr, "Inconsistent number of samples during 3-means.\n");
      exit (-1);
    }
  }

  /* Clean up and exit */
  free(cluster);
}

/***********************************************************************/
/* tree_start() - build a root node for a hierarchical partition tree  */
/*                                                                     */
/* This routine simply returns a pointer to a valid root node for the  */
/* hierarchical partition tree, based on the set of training patterns  */
/* passed in.                                                          */
/***********************************************************************/
treenode* tree_start(train_patterns, num_train_pats, d)
     double *train_patterns;	/* The training patterns array */
     int num_train_pats;	/* Number of elements in train_patterns */
     int d;			/* Dimensionality of the data */
{
  treenode *root;		/* The root of the tree to build */
  int i;			/* loop counter */

  /* allocate space for the root node */
  root = (treenode *) malloc (sizeof (treenode));

  /* Load its samples array */
  root->samples = (int *) malloc (num_train_pats * sizeof(int));
  for(i=0; i<num_train_pats; i++) root->samples[i] = i;
  root->num_samples = num_train_pats;

  /* Fill out the rest of the root data & return the new node */
  root->sample_mean = NULL;
  compute_mean(root, train_patterns, d);
  compute_radius(root, train_patterns, d);
  root->child[0] = (treenode *) NULL;
  root->child[1] = (treenode *) NULL;
  root->child[2] = (treenode *) NULL;

  return(root);
}

/***********************************************************************/
/* tree_build() - recursively build a hierarchical partitioning of the */
/*                training samples.                                    */
/*                                                                     */
/* The branch and bound knn algorithm depends upon the hierarchical    */
/* partition tree, which contains the index of each training pattern,  */
/* arranged in a tree of partitions.  Each level of the tree contains  */
/* a group of disjoint sets which comprise a partition of the training */
/* data.  By searching the tree using a few branch-and-bound cutoff    */
/* rules, we can speed the neighbor search considerably.  For more     */
/* details on the tree, see the reference noted in the first comment   */
/* block in this file.                                                 */
/*                                                                     */
/* The tree itself is a set of treenodes, starting with the 'root'     */
/* node.  (See defines.h for the structure of a treenode.              */
/* Also important is the dist_to_mean array, which contains the        */
/* distance from each training pattern to the mean (or center) of the  */
/* cluster associated with the leaf node for that pattern.             */
/*                                                                     */
/* Both the tree pointed to by 'root', and the dist_to_mean array are  */
/* built by this routine, and they will be passed to every routine     */
/* associated with neighbor searching from here on out.                */
/*                                                                     */
/* NOTE:  The root node should be created with tree_start(), and the   */
/*        space for the dist_to_mean array should be malloc'ed BEFORE  */
/*        calling this routine!                                        */
/***********************************************************************/
void tree_build(depth, root, patterns, d, dist_to_mean, method)
     int depth;			/* The depth to build the tree to */
     treenode *root;		/* The root of the subtree to build */
     double *patterns;		/* The pattern array */
     int d;			/* The dimensionality of the data */
     double *dist_to_mean;	/* The distance from each pattern to the */
				/* ...nearest cluster mean (OUTPUT) */
     int method;		/* Initial mean selection method for */
				/* the three-means algorithm */
{
  int i;			/* Loop counter */

  /* The termination criteria for this recursive routine is */
  /* halt when the depth is zero.                           */
  if (depth > 0) {

    /* If the depth is not zero, split the current node into three  */
    /* children, then call tree_build recursively for each child... */
    
    /* debug printout */
    if (DEBUG_LEVEL >= 2) {
      printf("\nCalling 3-means at level %d\n", depth);
      printf("---------------------------\n");
    }
    three_means(root, patterns, d, dist_to_mean, method);
    for(i=0; i<3; i++)
      tree_build(depth-1, root->child[i], patterns, d, dist_to_mean, method);
  }
}

/***********************************************************************/
/* search_tree - peform a branch and bound search on the hierarchical  */
/*               partition tree to find the nearest neighbors of the   */
/*               test_pattern.                                         */
/*                                                                     */
/* For more information on the branch and bound search, the partition  */
/* tree, and the cutoff rules, see the reference cited in the first    */
/* comment block in this file.                                         */
/***********************************************************************/

double search_tree(test_pattern, self, root, train_patterns, d,
		 dist_to_mean, k, nn, nndist, cut_val, dc)
     double *test_pattern;	/* The pattern to find the neighbors of */
     int self;			/* The pattern being voted on */
     treenode *root;		/* Root of the tree to search */
     double *train_patterns;	/* The training patterns array */
     int d;			/* The number of dimensions */
     double *dist_to_mean;	/* The dist. from ea. pattrn to cluster mean */
     int k;			/* The number of neighbors to find */
     int *nn;			/* The nearest neighbor array (OUTPUT) */
     double *nndist;		/* The distance to each nn (OUTPUT) */
     double cut_val;		/* The cutoff value (B) */
     int *dc;			/* The # of distance calculations needed */
{
  /* Local variables */
  double test_to_mean;		/* Distance from the test pattern to */
				/* the mean of the current cluster   */
  int sample;			/* The index of the sample */
  int cur_patnum;		/* The number of the current pattern */
  double *cur_pat;		/* The actual current pattern */
  double cur_dist_to_mean;	/* The dist_to_mean of the cur_pat */
  int cur_nn;			/* The current near neighbor */
  double cur_dist;		/* The distance to the current pattern */
  int i;			/* Loop counter */
  
  test_to_mean = sqrt(sqdist(test_pattern, root->sample_mean, d));

  /* Test for rule 1 -- See reference paper for details. */
  /* Note that this is one of only two places where we   */
  /* use the true euclidean distance.  Everywhere else   */
  /* we compute and store squared-euclidean-dist.        */
  if ((cut_val == -1.0) ||
      (test_to_mean <= cut_val + root->radius)) {

    /* If we pass and this is an internal node, evaluate the children */
    if (root->child[0])
      for(i=0; i<3; i++)
	cut_val = search_tree(test_pattern, self, root->child[i], train_patterns,
			      d, dist_to_mean, k, nn, nndist, cut_val, dc);
    
    /* If we pass rule 1, and this is a leaf node, evaluate the patterns */
    /* in the cluster associated with this leaf node.                    */
    else {
      for(sample=0; sample < root->num_samples; sample++) {
	cur_patnum = root->samples[sample];
	cur_pat = &(train_patterns[cur_patnum*d]);
	cur_dist_to_mean = dist_to_mean[cur_patnum];


	/* if SELF_VOTE is FALSE, then we assume the training           */
	/* and test patterns are the same, and do not allow self-voting */
	if (SELF_VOTE || cur_patnum != self) {

	  /* Test for rule 2 -- See reference paper for details.  */
	  /* This is the second place where we use true euclidean */
	  /* distance as opposed to squared-euclidean-distance.   */
	  if ((cut_val == -1.0) ||
	      (test_to_mean <= cut_val + cur_dist_to_mean)) {
	  
	    /* Only now do we do a distance check */
	    (*dc)++; /* Keep track of the # of dist calc.s we have done */
	    cur_dist = sqrt(sqdist(cur_pat, test_pattern, d));

	    /* If this pattern is closer than some near neighbor... */
	    for(cur_nn=0; cur_nn<k; cur_nn++)
	      if((nndist[cur_nn] == -1.0) || (cur_dist < nndist[cur_nn])) {
		
		/* ... make room in the list for it ... */
		for(i=k-2; i>=cur_nn; i--) {
		  nn[i+1] = nn[i];
		  nndist[i+1] = nndist[i];
		}
		
		/* ... and insert it into the nn list. */
		nn[cur_nn] = cur_patnum;
		nndist[cur_nn] = cur_dist;
		
		/* Break out of the cur_nn loop so that we don't add */
		/* this pattern to the nn list more than once.       */
		break;
	      }				/* If nndist...       */
	  }				/* Test rule 2        */
	}                               /* SELF_VOTE test     */
      }					/* For each pattern   */
    }					/* Leaf node ("else") */
  }					/* Test rule 1        */
  
  /* The new cut_val is the distance to the farthest nn */
  return nndist[k-1];
} 

/***********************************************************************/
/* predict() - find the nearest neighbors of the test pattern, and     */
/*             classify it.                                            */
/*                                                                     */
/* First the hierarchical partition tree is searched to quickly find   */
/* the k-nearest-neighbors of the test pattern.  Then the class of the */
/* neighbors are tallied, and the majority is the prediction for the   */
/* test pattern.                                                       */
/* NOTE:  a tie-breaking scheme is yet to be added.                    */
/*                                                                     */
/* The number of distance calculations is returned in distcalcs.       */
/*                                                                     */
/* Vote Tracking - if return_votes is a non-null pointer, then it is   */
/*                 assumed to point to an ALREADY-ALLOCATED array of   */
/*                 integers of size >= k, and the votes will be stored */
/*                 in this array.                                      */
/***********************************************************************/
int predict(test_pattern, self, train_patterns, num_train_pats, train_labels,
	    maxclass, class_count, root, dist_to_mean, d, k, distcalcs,
	    raw_votes)

     double *test_pattern;	/* The test pattern to classify */
     int self;			/* Number of the pattern being classified */
     double *train_patterns;	/* The training patterns */
     int num_train_pats;	/* The number of training patterns */
     int *train_labels;		/* The training pattern labels */
     int maxclass;		/* The largest training label */
     int *class_count;		/* The # patterns in each class */
     treenode *root;		/* The root of the partition tree */
     double *dist_to_mean;	/* The distance-to-closest-mean array */
     int d;			/* Dimensionality (# features) of the data */
     int k;			/* The k parameter (# nn's to vote) */
     int *distcalcs;		/* OUTPUT: # of dist. calculations done */
     int *raw_votes;		/* OUTPUT: if !NULL, the voting results */
{
  /* Static variables */
  static int *nn = NULL;	/* The nearest neighbor array */
  static double *votes = NULL;	/* The voting array */
  static double *nndist = NULL;	/* The distance to each near neighbor */

  /* Other local variables */
  int winner;			/* the class that got the most votes */
  int i, j;			/* loop counters */

  /* The first time the routine is called, we allocate space for */
  /* the various necessary arrays...                             */
  if (!nn) {
    nn = (int *) malloc (k * sizeof(int));
    nndist = (double *) malloc (k * sizeof(double));
    votes = (double *) malloc ((maxclass+1) * sizeof(double));
  }

  /* find the neighbors of the test_pattern */
  for(i=0; i<k; i++) nndist[i] = -1.0;
  search_tree(test_pattern, self, root, train_patterns, d, dist_to_mean, k, nn,
	      nndist, -1.0, distcalcs);
  
  /* Count up the raw and scaled votes for each near neighbor     */
  /* For scaled votes, divide each vote by the number of training */
  /* patterns in that class.                                      */
  for (i=0; i<=maxclass; i++)
    votes[i] = 0.0;
  if (raw_votes)
    for (i=0; i<=maxclass; i++)
      raw_votes[i] = 0;
  
  for (i=0; i<k; i++) {
#ifdef CBV_SCALING
    votes[train_labels[nn[i]]] += 1.0/(double)class_count[train_labels[nn[i]]];
#else
    votes[train_labels[nn[i]]]++;
#endif
    /* Return the voting information, if requested */
    if (raw_votes)
      raw_votes[train_labels[nn[i]]]++;
  }

  /* Find the winner */
  winner = 0;
  for (i=0; i<=maxclass; i++)
    if(votes[i] > votes[winner])
      winner = i; 

  /* Debugging printouts */
  if (DEBUG_LEVEL >= 3) {
    printf("Test pattern: ");
    for (i=0; i<d; i++) printf("%lf\t", test_pattern[i]);
    printf("\nNeighbors: ");
    for (i=0; i<k; i++) {
      printf("%d", nn[i]);
      if (i<k-1) printf(", ");
    }
    printf("\tDistances: ");
    for (i=0; i<k; i++) {
      printf("%.4lf", nndist[i]);
      if(i<k-1) printf(", ");
    }
    printf("\n");
  }
  
  /* Tiebreaking scheme:  nearest neighbor wins */
  for (i=0; i <= maxclass; i++) {

    /* If we find a class that got the same number of votes as the */
    /* current winner, we need to run a tie-breaker */

    if ((i != winner) && (votes[i] == votes[winner])) {

      /* Since the nearest neighbor list is sorted by distance,   */
      /* we can run through it and choose the class which matches */
      /* the first neighbor...                                    */

      for (j=0; j<k; j++) {

	/* If the first neighbor found matches the current winner, */
	/* leave things alone.                                     */
	if (train_labels[nn[j]] == winner)
	  break;

	/* if the first neighbor found matches the challenger (i), */
	/* then i becomes the current winner.                      */
	else if (train_labels[nn[j]] == i) {
	  winner = i;
	  break;
	}
      }
    }
  }

  return(winner);
}

/**********************************************************************/
/* usage() -  Print out a usage message and a user-supplied string to */
/*            standard out.                                           */
/**********************************************************************/
void usage(outstr, progname)
     char *outstr;
     char *progname;
{
  /* Output a usage message */

  printf("\n");

  if (outstr)
    if (strlen(outstr) > 0)
      printf("%s\n", outstr);

  printf("Usage:\t%s [-v] <k-value>\n", progname);
  printf("\t%s [-v] <k-value> training-file test-file y|n [weight-file]\n",
	 progname);
  printf(" (If the input files are specified on the command line, the\n");
  printf("  y|n argument indicates whether the test data is labelled.)\n");

  printf("\n");
}
    
/***********************************************************************/
/* MAIN ROUTINE:                                                       */
/***********************************************************************/

int main(argc, argv) 
     int argc;
     char *argv[];
{

  /* -------------------------- Local Variables: -------------------------- */
  double *train_patterns;	/* The training pattern array */
  double *test_patterns;	/* The test pattern array */
  double *dist_to_mean;		/* Dist. from each trning pat to clust mean */
  double *test_pattern;		/* The current test pattern */
  double *weights;		/* The feature-axis weights (if any) */
  int *train_labels;		/* The class labels of the training data */
  int *test_labels;		/* Class labels of test data, if known */
  int *predicted;		/* The prediction for each test pattern */
  int *votes_array;		/* Array for tracking votes */
  int *votes;			/* Index into votes_array */
  int *class_count_trn;		/* # of training patterns of each class */
  int k;			/* The k-parameter for the knn */
  int num_train_pats;		/* The number of training patterns */
  int num_test_pats;		/* The number of testing patterns */
  int num_features;		/* The number of training features */
  int num_test_features;	/* Should equal num_features */
  int result;			/* A result flag */
  int curpat;			/* The index of the current pattern */
  int distcalcs;		/* # of dist. calculations required */
  int distcalcsum;		/* sum of dc for each iteration  */
  int correct_count;		/* The number of correct predictions */
  int maxclass;			/* Maximum class label of training data */
  int minclass;			/* Minimum class label of training data */
  int labelled;			/* Flag: is the training data labelled? */
  int weighted;			/* Flag: use a set of axis weights? */
  int interactive;		/* Flag: interactive mode or not? */
  int show_votes;		/* Flag: print the voting results? */
  int next_arg;			/* The next command-line argument */
  int i, j;			/* Loop counters */
  char filename[MAXFILE+1];	/* A file name string buffer */
  char inbuf[MAXFILE+1];	/* A buffer for user input */
  treenode *root;		/* The root of the partition tree */
  FILE *fptr;			/* A file pointer */
  clock_t cpu_start;		/* These last variables are all for...  */
  clock_t cpu_stop;		/* ...checking the efficiency of the    */
  double cpu_time;		/* ...algorithm.  The clock_t variables */
  time_t wall_start;		/* ...are for checking CPU time, while  */
  time_t wall_stop;		/* ...the time_t variable are for wall  */
  double wall_time;		/* ...clock time.                       */
  /*------------------------------------------------------------------------*/

  /* Check the usage */
  if (argc < 2) {
    usage("", argv[0]);
    return( 1 );
  }

  /* Print a version string */
  printf ("%s\n", BBKNN_VERSION);

  next_arg = 1;

  /* Check if the -v flag is present:                                    */
  /* If the -v flag is present, and any other arguments follow it, then  */
  /* the filenames must be included on the command line, so we are not   */
  /* in interactive mode.  If the -v flag is not present, and there are  */
  /* other arguments after the k-value, then the filenames must be       */
  /* included on the command line, so we are likewise not in interactive */
  /* mode.                                                               */
  if (argc > 2) {
    if (!strcmp(argv[next_arg], "-v")) {
      show_votes = TRUE;
      next_arg++;
      /* -v + other arguments */
      if (argc > 3)
	interactive = FALSE;
      /* -v only */
      else
	interactive = TRUE;
    }
    /* no -v, more than 2 arguments */
    else {
      show_votes = FALSE;
      interactive = FALSE;
    } 
  }
  else {
    /* less than 2 arguments */
    show_votes = FALSE;
    interactive = TRUE;
  }

  /* Get the k-value from the command line */
  result = sscanf(argv[next_arg++], "%d", &(k));
  if ((result != 1) || (k <= 0)) {
    usage("K-value must be a positive integer", argv[0]);
    return(1);
  }

  /* Start the CPU clock */
  clock();

  /* Open the training pattern file*/
  if (interactive) {
    strcpy(filename, "train.dat");
    printf("Enter filename for training data [%s]\n", filename);
    fgets(inbuf, MAXFILE, stdin);
    if (!isblank(inbuf))
      sscanf(inbuf, "%s", filename);
  }
  else
    strcpy(filename, argv[next_arg++]);
  
  fptr = fopen(filename, "r");
  if (!fptr) {
    fprintf(stderr, "Unable to open input file: %s\n", filename);
    return(1);
  }

  /* Read the training patterns */
  result = read_patterns(fptr, &train_patterns, &train_labels,
			 &num_train_pats, &num_features, TRUE);
  if (result == FAILURE) {
    fprintf(stderr, "Unable to read training patterns.\n");
    return(1);
  }
  fclose(fptr);

  /* Allocate space for the dist_to_mean array */
  dist_to_mean = (double *) malloc(num_train_pats * sizeof(double));

  /* Determine the maximum class label */
  maxclass = 0;
  minclass = -1;
  for (i=0; i < num_train_pats; i++) {
    if (train_labels[i] > maxclass)
      maxclass = train_labels[i];
    if ((train_labels[i] < minclass) || (minclass == -1))
      minclass = train_labels[i];
  }

  /* Count the number of training patterns of each type */
  class_count_trn = (int *) malloc ((maxclass+1) * sizeof(int));
  for (i=0; i<= maxclass; i++)
    class_count_trn[i] = 0;
  for (i=0; i < num_train_pats; i++)
    class_count_trn[train_labels[i]]++;

  /* Now load the testing patterns... */
  /* Open the test pattern file*/
  if (interactive) {
    strcpy(filename, "test.dat");
    printf("Enter filename for test data [%s]\n", filename);
    fgets(inbuf, MAXFILE, stdin);
    if (!isblank(inbuf))
      sscanf(inbuf, "%s", filename);
  }
  else
    strcpy(filename, argv[next_arg++]);
  
  fptr = fopen(filename, "r");
  if (!fptr) {
    fprintf(stderr, "Unable to open input file: %s\n", filename);
    return(1);
  }

  /* Ask if the test data is labelled */
  if (interactive) {
    inbuf[0] = '\0';
    while ((inbuf[0] != 'y') && (inbuf[0] != 'n')) {
      printf("\nIs the test data labelled?  ");
      fflush(stdout);
      scanf("%s", inbuf);
      if(isupper(inbuf[0]))
	inbuf[0] = tolower(inbuf[0]);
    }
  }
  else
    strcpy(inbuf, argv[next_arg++]);

  if (inbuf[0] == 'y') labelled = TRUE;
  else labelled = FALSE;
  
  /* Read the test patterns */
  result = read_patterns(fptr, &test_patterns, &test_labels,
			 &num_test_pats, &num_test_features, labelled);
  if (result == FAILURE) {
    fprintf(stderr, "Unable to read test patterns.\n");
    return(1);
  }
  fclose(fptr);

  /* Make sure the number of features matches */
  if (num_test_features != num_features) {
    fprintf(stderr, "Number of training and test features does not match.\n");
    return(1);
  }

  /* If vote tracking enabled, allocate vote-tracking array */
  if (show_votes)
    votes_array = (int *) malloc (num_test_pats * (maxclass+1) * sizeof(int));

  /* See if the user desires a weighted-knn run... */
  weighted = FALSE;
  if (!interactive && (argc >= next_arg+1)) {
    weighted = TRUE;
    strcpy (filename, argv[next_arg++]);
  }
  else if (interactive) {
    inbuf[0] = '\0';
    while ((inbuf[0] != 'y') && (inbuf[0] != 'n')) {
      printf("\nUse a feature-weights file?  ");
      fflush(stdout);
      scanf("%s", inbuf);
      if(isupper(inbuf[0]))
	inbuf[0] = tolower(inbuf[0]);
    }
    fflush(stdin);
    if (inbuf[0] == 'y') { 
      weighted = TRUE;
      strcpy(filename, "weight.dat");
      printf("Enter filename for feature weights [%s]\n", filename);
      fgets(inbuf, MAXFILE, stdin);
      if (!isblank(inbuf))
	sscanf(inbuf, "%s", filename);
    }
  }
  
  /* If feature-weights are provided, load them */
  if (weighted) {
    weights = (double *) malloc (num_features * sizeof (double));
    
    fptr = fopen(filename, "r");
    if (!fptr) {
      fprintf(stderr, "Unable to open input file: %s\n", filename);
      return(1);
    }

    for (i=0; i< num_features; i++)
      if (fscanf(fptr, "%lf", &(weights[i])) < 1) {
	fprintf (stderr, "Not enough weights in %s\n", filename);
	return(1);
      }
    fclose(fptr);
  }

  printf("Patterns loaded, building partition tree...\n");

  /* Scale both training and testing patterns, if necessary */
  if (weighted) {
    wall_start = time(NULL);
    cpu_start = clock();
    
    for (curpat = 0; curpat < num_train_pats; curpat++)
      for (i=0; i<num_features; i++)
	train_patterns[curpat*num_features +i] *= weights[i];

    for (curpat = 0; curpat < num_test_pats; curpat++)
      for (i=0; i<num_features; i++)
	test_patterns[curpat*num_features +i] *= weights[i];

    cpu_stop = clock();
    wall_stop = time(NULL);
    wall_time = wall_stop - wall_start;
    cpu_time = cpu_stop - cpu_start;
    printf("Pattern scaling:  CPU=%.2lfs, Wall=%.2lfs\n",
	   cpu_time/(double)CLOCKS_PER_SEC, wall_time);
  }
  
  /* Build the root node */
  root = tree_start(train_patterns, num_train_pats, num_features);
  
  /* Recursively build the rest of the hierarchy  */
  /* The first parameter is the desired tree depth - this parameter really */
  /* depends on the number of training samples.  In the reference, they    */
  /* used a tree depth of 7, which is what I have initially used here.     */
  /* Some empirical experimentation will be necessary to find a good rule  */
  /* of thumb on the appropriate size for the tree depth.                  */
  /*                                                                       */
  /* The last parameter specifies the initial center selection method, */
  /* see the three_means() inline documentation for details...         */
  wall_start = time(NULL);
  cpu_start = clock();
  tree_build(7, root, train_patterns, num_features, dist_to_mean, FAST); 
  cpu_stop = clock();
  wall_stop = time(NULL);
  wall_time = wall_stop - wall_start;
  cpu_time = cpu_stop - cpu_start;
  printf("Tree build:  CPU=%.2lfs, Wall=%.2lfs\n",
	 cpu_time/(double)CLOCKS_PER_SEC, wall_time);
  
  /* Classification: */
  /* --------------- */
  printf("Test data loaded, classifying...\n");
  predicted = (int *) malloc (num_test_pats * sizeof(int));

  /* Keep track of how long and how many distance calculations */
  /* it takes...                                               */
  distcalcsum = 0;
  wall_start=time(NULL);
  cpu_start=clock();
  
  /* Classify each test pattern */
  for (curpat=0; curpat < num_test_pats; curpat++) {
    test_pattern = &(test_patterns[curpat*num_features]);
    if (show_votes)
      votes = &(votes_array[curpat*(maxclass+1)]);
    else
      votes = (int *) NULL;
    distcalcs = 0;
    predicted[curpat] = predict(test_pattern, curpat, train_patterns,
				num_train_pats, train_labels, maxclass,
				class_count_trn, root, dist_to_mean,
				num_features, k, &distcalcs, votes);
    distcalcsum += distcalcs;
  }

  /* Check the execution times */
  cpu_stop = clock();
  wall_stop = time(NULL);
  cpu_time = (cpu_stop - cpu_start)/(double)CLOCKS_PER_SEC;
  wall_time = wall_stop - wall_start;

  /* Print some performance data */
  printf("Classification completed.  CPU=%.2lfs, Wall=%.2lfs\n",
	 cpu_time, wall_time);
  printf("%d distance calculations required.\n\n", distcalcsum);

  /* Print detailed results information */
  correct_count = 0;
  
  if (show_votes) {
    printf("votes:  ");
    for (i=minclass; i<=maxclass; i++)
      printf("[class %d] ", i);
    printf("\n");
  }

  for (i=0; i<num_test_pats; i++) {
    printf("Pattern %4d -- ", i);
    if (labelled)
      printf("observed: %d  ", test_labels[i]);
    printf("predicted: %d", predicted[i]);
    
    if (show_votes) {
      printf("  votes:");
      for (j=minclass; j<=maxclass; j++)
	printf(" %d", votes_array[(i*(maxclass+1))+j]);
    }

    printf("\n");

    /* If we already know the test labels, track summary data */
    if (labelled) {
      if (predicted[i] == test_labels[i])
	correct_count++;
    }
  }
 
  /* Print the summary data, if we were tracking it */
  if (labelled)
    printf("%d/%d correct predictions (%.2lf%%)\n", correct_count,
	   num_test_pats, ((double)correct_count/(double)num_test_pats)*100.0);

  return(0);
}
