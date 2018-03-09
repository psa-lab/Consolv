/***********************************************************************/
/* defines.h - type definitions and constants for the bbknn classifier */
/***********************************************************************/

#ifndef _DEFINES_H
#define _DEFINES_H

/*------------------------*/
/* Necessary header files */
/*------------------------*/
#include "util.h"

/*----------------------------*/
/* Generally useful constants */
/*----------------------------*/
#define FAILURE  0
#define SUCCESS -1

#define TRUE -1
#define FALSE 0

/*-------------*/
/* Flag values */
/*-------------*/

/* Values for the three-means 'method' flag */
#define UNIFORM 1
#define FAST 2

/*----------------------------*/
/* Adjustable BBKNN constants */
/*----------------------------*/

/* The maximum number of iterations for the 3means clustering */
#define ITERATION_LIMIT 20

/* The max number of samples in a leaf node of the partition hierarchy */
#define MAX_LEAF_SAMPLES 50

/* Maximum number of characters/line in the input files: */
#define MAXLINE 1000

/* Maximum number of characters in a file name */
#define MAXFILE 256

/* Memory is allocated in blocks of MBLOCK*sizeof(double) units */
#define MBLOCK 1000

/* The DEBUG_LEVEL controls the level of diagnostic message output */
/* 0 = none, 1 = a little, 2 = a lot, 3 = dump everything          */
#define DEBUG_LEVEL 0

/*-------------------------*/
/* Typedefs and structures */
/*-------------------------*/
typedef struct _treenode 
{
  int *samples;			/* The set of samples for this node */
  int num_samples;		/* The number of samples for this node */
  double *sample_mean;		/* The mean vector for this set of samples */
  double radius;		/* The max dist from a point to the mean */
  struct _treenode *child[3];	/* The three successors to this node */
} treenode;


/* Do not adjust the following definitions */
#define SELF_VOTE TRUE
#undef CBV_SCALING

#endif
