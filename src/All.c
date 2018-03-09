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

/***************************************************************************/
/* all.c - Calculate the atomic density, atomic hydrophilicity, Bond Value,*/
/*          number of H-Bonds, mobility and sum of BVals of nonwater atoms */
/*          near each water for a set of waters in the context of a protein*/
/*                                                                         */
/* Synopsis                                                                */
/* --------                                                                */
/*  This program takes two input files, one is a PDB file (*.xfm)containing*/
/*  a protein, and another is a smaller file in PDB (*.xfm.hits)           */
/*  format which contains the HETATM records of the waters to be           */
/*  considered.  For each water present in the waters file the atomic      */
/*  hydrophilicity will be computed as follows:                            */
/*   1) Find all atoms within HY_CUTOFF angstroms of the of the water      */
/*      molecule                                                           */
/*   2) For each atom compute the atomic hydrophilicity value --           */
/*        neutral oxygen: 0.528 (hydrations/occurrence)                    */
/*        negatively charged oxygen: 0.511                                 */
/*        positively charged nitrogen: 0.444                               */
/*        neutral nitrogen: 0.351                                          */
/*        all carbon, sulfur atoms: 0.078                                  */
/*   3) Compute the atomic hydrophilicity for this water molecule as the   */
/*      sum of the values found in step 2.                                 */
/*                                                                         */
/*  Atomic density is computed as:					   */
/*   1) Find all atoms withen HY_CUTTOFF angstroms of the water molecule   */
/*   2) For each water atom, check its residue ---			   */
/*      If residue is H2O,HOH,DOD,D2O,WAT DO NOT count it as a neigbour    */
/*	If the residue is something else, it IS a valid neighbour, hence   */
/*      count it.							   */
/*   3) Repeat the step 2 for each Water Atom in .HITS file i.e.           */
/*	for(each atom in .hits file )					   */
/*	  								   */
/*	    check(each atom in .xfm, for cutoff and residue type)	   */
/*	    if( atom in .xfm != Water moleucle )			   */
/*		atomic_density ++;					   */
/*	   								   */
/*  Hbond # is computed as follows:					   */
/*   1) Find all atoms and withen HBOND_CUTTOF angstroms of the water      */
/*      molecule. Hbonds to waters and to proteins are tabulated separately*/
/*   2) For each water atom, check its residue ---			   */
/*      If residue is H2O,HOH,DOD,D2O, DO NOT count it as a candidate for  */
/*      for a hbond. The only candidates are: N*,X* (halides),O*,S*	   */
/*      STEP 2 IS NOT APPLICABLE for HBonds to water                       */
/*   3) Repeat the step 2 for each Water Atom in .HITS file. Its similar to*/
/*      step 3 in Atomic Density above					   */
/*   Note: HBonds of lengths upto 4.0 with Sulphur are allowed             */
/*         & No hbonds are allowed betwen a water and Nitrogen in PRO      */
/*         Hbonds to protein atoms are calculated separately from HBonds to*/
/*         waters                                                          */
/*									   */
/*  Mobility is calculated as follows:					   */
/*   1) While read the FREE PDB file, if the current record read is a HETAM*/
/*      Check whether its water or no. If it is water, add the Bval and    */
/*      and mobility values to an accumulating variable and increment the  */
/*	counter.							   */
/*   2) While actually processing, the scaled mobility is given as	   */
/*									   */
/*	Mob=(Current_Bval/Avg_Bval)/(Current_Mob/Avg_Mobility)		   */
/*                                                                         */
/*  Sum_of_Bvals: It is the sum of Bvals of all nonwater atoms withen      */
/*    HY_CUTOFF distance of each water molecule.                           */
/*									   */
/*  Avg_of Bvals: It is the Avg of Bvals of all nonwater atoms withen      */
/*    HY_CUTOFF distanec of each water molecule. ie                        */
/*      Avg = Sum_of_Bvals/Adn                                             */
/*                                                                         */
/* Programmer - Michael L. Raymer  9/20/94                                 */
/* * adn.c and ahphil.c combined by Vishal Thakkar 12/12/96                */
/* * hbond calculation module added by Vishal Thakkar 12/16/96             */
/* * Mobility calculation and formating done by Vishal Thakkar 12/17/96    */
/* * Sum_of_Bvals calcalation module added by Vishal Thakkar 12/26/96      */
/* * Final hbond calculation added by Vishal, April 14, 97                 */
/***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>


/* Here are some global constants, if this grows into a mult-module */
/* application, then these will be moved to a global.h file         */
#define HY_CUTOFF 3.6		/* "Neighbor" distance cutoff */
#define HBOND_CUTOFF 3.5        /* "Neighbor" distance cutoff for Hydrogen Bonds */
#define HBOND_CUTOFF_S 4.0	/* "Neighbor" distance cutoff for Hydrogen Bonds with Sulfur */
#define MAXLN 132		/* Max input line length */
#define BLANK ' '		/* Blank character */
#define TABCHAR '	'	/* Tab character */

/*#define AHP_LIST_FILE "/home/thakkarv/work/ahplist.dat"*/

/* PDB file format for ATOM and HETATM records... */
#define PDB_ANUM_COL 7		/* Atom Serial Number Column */
#define PDB_ANUM_SZ 5		/* Atom Serial Number Size */
#define PDB_RESNUM_COL 23	/* Residue Sequence Number Column */
#define PDB_RESNUM_SZ 5		/* Residue Sequence Number Size */
#define PDB_ANAM_COL 13		/* Atom Name Column */
#define PDB_ANAM_SZ 4		/* Atom Name Size */
#define PDB_RESNAM_COL 18	/* Residue Name Column */
#define PDB_RESNAM_SZ 3		/* Residue Name Size */

#define PDB_XCRD_COL 31         /* X coordinate */
#define PDB_CODE_COL 73		/* The PDB code Column */
#define PDB_CODE_SZ 4		/* The PDB code Size */
#define PDB_BVAL_COL 61		/* The BVal Column */

/* Run Mode types */
#define APP_MODE 1
#define TEST_MODE 0

/* Boolean values */
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

/* H Bond type */
#define WITH_PROTEIN 1
#define WITH_WATER 2

/* Structs needed globally */
typedef struct ahp_rec_type {
  char atom[PDB_ANAM_SZ+1];	     /* The Atom Name */
  char residue[PDB_RESNAM_SZ+1];   /* The Resdiue Name */
  float ahp;			     /* The Atomic Hydrophilicity */
  struct ahp_rec_type *next;	     /* List link field */
} ahp_rec;


/* Debug Variable - may not be fully implemented */
/* 0 = silence is golden                         */
/* 1 = occasional informational messages         */
/* 2 = frequent status reports                   */
/* 3 = spill my guts                             */
/* --------------------------------------------- */
int debug_level = 2;

/***************************************************************************/
/* SUBROUTINES AND FUNCTIONS                                               */
/***************************************************************************/

/***************************************************************************/
/* get_file - get a filename from the user, and open the file              */
/*                                                                         */
/* get_file prompts the user for a filename, when the user enters a        */
/* non-empty string, get_file will attempt to open the file.  It will      */
/* continue to prompt until a valid filename is inputted and the file is   */
/* opened.                                                                 */
/*                                                                         */
/* INPUTS:   prmpt - string to prompt the user with.                       */
/*           omode - the mode to open the file in ("r", "rw", etc...)      */
/*           fdefault - the default file name                              */
/*                      if no default is desired this should point to a    */
/*                      string containing only a null character ('\0')     */
/*                                                                         */
/* OUTPUTS:  get_file returns a pointer to the open file.                  */
/***************************************************************************/
FILE *get_file(prmpt, omode, fdefault,interactive)
     char *prmpt;		/* String to prompt user with */
     char *omode;		/* Mode to open the file in */
     char *fdefault;		/* The default file name, if any */
     int interactive;		/* interactive=1 then ask otherwise don't */
{
  /* Local Variables */
  char input_fn[MAXLN + 1];	   /* User inputted filename */
  FILE *fptr;			   /* Pointer to opened file */

  /* Initializations */
  input_fn[0] = '\0';
  fptr = (FILE *) NULL;

  while (!fptr)
    {
      if(interactive==FALSE) 
	strcpy(input_fn,fdefault);
      else
	/* Get a valid filename */
	while (strlen(input_fn) <= 1){
	    /* prompt user */
	    printf("\n%s\n", prmpt);
	    if (strlen(fdefault) >= 1)
	      printf("(press return to use [%s])\n", fdefault);
	    
	    /* get the filename */
	    fgets(input_fn, MAXLN - 1, stdin);
	    
	    /* If empty, use the default, if any */
	    if ((strlen(input_fn) <= 1) && (strlen(fdefault) >= 1))
	      strcpy(input_fn, fdefault);
	}

      /* Clean up any trailing returns */
      if (input_fn[strlen(input_fn) - 1] == '\n')
	input_fn[strlen(input_fn) - 1] = '\0';

      /* Check, if the file is open for writing are we overwriting it or no*/
      if(*omode=='w')
	{
	  struct stat b;
	  char ans;
	  char bakname[MAXLN+1];

	  if(stat(input_fn,&b)==NULL){/* File exits */
	    if(interactive==TRUE) {/* Go interactive */
	      printf("file: %s exits, Overwrite [y/n] ?",input_fn);
	      scanf("%c",&ans);
	      fflush(stdin);  /*Flush the input buffer of any carriage returns */
	      ans=toupper(ans);
	      if(ans=='N')
		return(NULL); /* Return Null pointer */
	    }  
	    else { /* Non interactive */
	      strcpy(bakname,input_fn);
	      strcat(bakname,".bak");
	      printf("\nWarning: renaming %s to %s\n", input_fn, bakname);
	      rename(input_fn,bakname);
	    }
	  }
	}

      /* Attempt to open the file */
      fptr = fopen(input_fn, omode);

      /* If we fail, reset the filename and print an error */
      if (!fptr)
	{
	  printf("\nUnable to open file: %s\n\n", input_fn);
	  input_fn[0] = '\0';
          if(interactive==FALSE) {
		printf("Exiting...\n");
	        exit(-1);
          }
	}
    }
  return(fptr);
}

/***************************************************************************/
/* load_list - load a list of ahp values from a file into a linked list.   */
/*                                                                         */
/* Synopsis - the file should consist of only a) blank lines, b) lines     */
/*            in which the first character is a # (comments) and c) lines  */
/*            with data of the form ATOM RESIDUE VALUE, with any amount    */
/*            of white space between the three fields.                     */
/*                                                                         */
/* Returns - list_head will point to the start of the list.                */
/***************************************************************************/
void load_list(list_head)
     ahp_rec **list_head;
{
  /* Local Variables */
  ahp_rec *currec, *newrec;	/* Current/New list elements */
  FILE *list_file;		/* Points to the list file */
  char inbuf[MAXLN+1];		/* Input file line buffer */
  char *firstchar;		/* First character in input file line */

  /* Open the list file */
  list_file = fopen(AHP_LIST_FILE, "r");
  
  /* If we fail, reset the filename and print an error */
  if (!list_file) {
    printf("\nUnable to open the AHP values data file.\n");
    printf("File ahplist.dat not found.\n\n");
    exit(-1);
  }

  /* Set up the list */
  *list_head = NULL;
  currec = NULL;
  

  /* Load the data */
  while (!feof(list_file)) {
    /* Get a line */
    if (fgets(inbuf, MAXLN, list_file)) {
      firstchar = inbuf;
      /* Find the first non-whitespace */
      while ((*firstchar == BLANK) || (*firstchar == TABCHAR))
	firstchar++;

      /* If this is a valid line, load it */
      if ((*firstchar != '#') &&
	  (*firstchar != '\n') &&
	  (strlen(firstchar) > 0)) {
	
	/* Allocate the space for the line */
	newrec = (ahp_rec *) malloc (sizeof(ahp_rec));
	newrec->next = NULL;
	if (!currec) {
	  currec = newrec;
	  *list_head = newrec;
	}
	else {
	  currec->next = newrec;
	  currec = newrec;
	}

	/* Load the data */
	sscanf(inbuf, "%s %s %f", currec->atom, currec->residue,
	       &(currec->ahp));
      }
    }
  }
}

/***************************************************************************/
/* free_list - free the linked list of ahp values.                         */
/*                                                                         */
/* WARNING: Does not reset the input pointer to NULL, so only use during   */
/*          final cleanup.                                                 */
/***************************************************************************/
void free_list(lptr)
     ahp_rec *lptr;
{
  /* local variables */
  ahp_rec *tptr;		/* Temp pointer */

  while (lptr) {
    tptr = lptr;
    lptr = lptr->next;
    free(tptr);
  }
}

/*************************************************************************/
/* This function, iswater returns TRUE if the string passed to it is     */
/* WAT,H2O,HOH,D2O,DOD else returns a false                              */
/*************************************************************************/

int iswater(char *str)
{

 if (strcmp(str, "HOH")==NULL)
    return TRUE;
  else if (strcmp(str, "H2O")==NULL)
    return TRUE;
  else if (strcmp(str, "WAT")==NULL)
    return TRUE;
  else if (strcmp(str, "DOD")==NULL)
    return TRUE;
  else if (strcmp(str, "D2O")==NULL)
    return TRUE;
  else 
    return FALSE;
}
/*************************************************************************/
/* This function, hbondpossible works as follows                         */
/* Inputs parameters: aname, resname                                     */
/*                    aname: char *, has the name of the atom            */
/*		      resname: char *, has the name of the residue       */
/* What it does:      it checks if the aname is CL,F,I,BR,N,O,S. If it is*/
/*                    then check that O does not belong to a water.      */
/*                    If the above condition is satisfied, return a 1 else*/
/*                    return a zero                                      */
/*************************************************************************/
int hbondpossible(aname,resname, dist)
      char *aname;
      char *resname;
      float dist;
{
  /* Local Variables */
  char atom[PDB_ANAM_SZ+1];         /* Atom name in CAPS */
  char residue[PDB_RESNAM_SZ+1];    /* Residue name in CAPS */
  int possible=0;                   /* Hbond possible or not */
  int i=0;


 /* Proceed only if HBONDS are possible */
  if( 0 < dist && dist <= HBOND_CUTOFF_S) {
    /*Copy the atom and residue names to local variables, and upcase them */
    strncpy(atom,aname,PDB_ANAM_SZ);
    while(atom[i]!=NULL){
	if(isalpha(atom[i]))
	   atom[i]=toupper(atom[i]);
	i++;
    } 
    i=0;
    strncpy(residue,resname,PDB_RESNAM_SZ);
    while(residue[i]!=NULL) {
      if(isalpha(residue[i])) {
	residue[i]=toupper(residue[i]);
      }
      i++;
    }
    /* Till this point, atom and residue contain uppercase atomname and residue  name */
    if( dist <= HBOND_CUTOFF) {
      if(iswater(residue) == TRUE)
	possible=WITH_WATER;
      else { /* If the residue is not water, proceed to check for Halides, N,O,S */
	   /* No hbonds with Nitrogen of PRO */
        if( (strncmp(atom,"N",1)==NULL) && strncmp(residue,"PRO",3)==NULL)
	  possible=0;
	else if(strncmp(atom,"CL",2)==NULL)
	  possible=WITH_PROTEIN;
	else if(strncmp(atom,"F",1)==NULL)
	  possible=WITH_PROTEIN;
	else if(strncmp(atom,"I",1)==NULL)
	  possible=WITH_PROTEIN;
	else if(strncmp(atom,"BR",2)==NULL)
	  possible=WITH_PROTEIN; 
	else if(strncmp(atom,"N",1)==NULL)
	  possible=WITH_PROTEIN;
	else if(strncmp(atom,"O",1)==NULL)
	  possible=WITH_PROTEIN;
	else if(strncmp(atom,"S",1)==NULL)                     
	  possible=WITH_PROTEIN;
	else
	  /* If none of the above, no bonds possible */
	  possible=0;
      }
    }
    if( HBOND_CUTOFF < dist && dist <= HBOND_CUTOFF_S)    /* Extra checking with this if!*/
      if(strncmp(atom,"S",1)==NULL)                       /* This fraction allows S-H bonds of */
	possible=WITH_PROTEIN;                            /* lengths upto HBOND_CUTOFF_S */


  }
 return possible;
}
/* end of hbondpossible  */

/*************************************************************************/
/* This sub is copied from adn.c                                         */
/* We have rename hphil adn                                              */
/* adn is passed an atom name and a residue name. The atom name is not   */
/* used but the residue name is used.                                    */
/* function: If residue name =  H2O,HOH,WAT,D2O,DOD  return 0            */
/*           else return = 1                                             */
/*************************************************************************/
int adn(anam, resnam)
     char *anam;
     char *resnam;
{
  /* Local Variables */
  char *cpos;			/* A multi-purpose character pointer */
  char atom[PDB_ANAM_SZ+1];	/* The Atom Name of the current atom */
  char residue[PDB_RESNAM_SZ+1]; /* The Resdiue Name - current atom */
  int adnval;			/* The atomic density value */

  /* Copy the atom and residue names to local variables, and upcase them */
  /* ------------------------------------------------------------------- */

  /* Atom name  : Conversion to Upper case*/
  strncpy(atom, anam, PDB_ANAM_SZ);
  cpos=atom;
  while(*cpos)
    {
      if (isalpha(*cpos))
	*cpos = toupper(*cpos);
      cpos++;
    }
  atom[PDB_ANAM_SZ] = '\0';

  /* Residue name :Conversion to Upper case */
  strncpy(residue, resnam, PDB_RESNAM_SZ);
  cpos=residue;
  while(*cpos)
    {
      if (isalpha(*cpos))
	*cpos = toupper(*cpos);
      cpos++;
    }
  residue[PDB_RESNAM_SZ] = '\0';

  /* Compute the adnval for this atom... */
  /* ---------------------------------- */

  if(iswater(residue)==TRUE)
    adnval=0;
  else
    adnval = 1;
	
  /* If we are debugging, print the value */
  if (debug_level >= 3)
    printf("Atomic Density:\tatom=%s\t\tResidue=%s\tValue=%d\n", atom, residue, adnval);

  /* Return the adnval and exit */
  return(adnval);
  
}
/**** End of adn ***/

/***************************************************/
/* This routine takes in the atom name and residue */
/* and corrects them. This is required if that pair*/
/* of atom-residue do not have an entry in the ahp-*/
/* list.dat                                        */
/***************************************************/
void correctName(char *a, char *r){

  int i=0;
  
  if(a[0]=='A') i++;		/* Some atom names start with A for eg, ACN, AN1,etc*/
  if(isalpha(a[i])==0)
    i++;	/* Find the first non 'A' char*/
 
  /* Zinc*/
  if(a[i]=='Z' && a[i+1]=='N') {
    strncpy(a,"ZN",2);a[2]=0;
    strncpy(r,"ZN",2);r[2]=0;
    return;
  }

  /* Sodium*/
  if(a[i]=='N' && a[i+1]=='A') {
    strncpy(a,"NA",2);a[2]=0;
    strncpy(r,"NA",2);r[2]=0;
    return;
  }

  /* Chlorine*/
  if(a[i]=='C' && a[i+1]=='L') {
    strncpy(a,"CL",2);a[2]=0;
    strncpy(r,"AVG",3);r[PDB_RESNAM_SZ]=0;
    return;
  }

  /*Bromine*/
  if(a[i]=='B' && a[i+1]=='R') {
    strncpy(a,"BR",2);a[2]=0;
    strncpy(r,"AVG",3);r[PDB_RESNAM_SZ]=0;
    return;
  }

  /* All other atoms, the first  non-A char is used as atom name
   * and residue of AVG */
  a[0]=a[i];
  a[1]=0;
  
 strncpy(r,"AVG",3);r[PDB_RESNAM_SZ]=0;
 return;
}



/***************************************************************************/
/* hphil - compute the atomic hydrophilicity value for a specific atom     */
/*                                                                         */
/* Brief synopsis - this routine will take an atom and match the atom      */
/*                  name and residue against the ahp_list.  If the atom    */
/*                  name and residue are found, the value from the         */
/*                  ahp_list is returned.  If not, a warning is printed    */
/*                  and a value of 0.0 is returned.                        */
/*                                                                         */
/*            To compute atomic density, just have this routine always     */
/*            return 1.                                                    */
/***************************************************************************/
float hphil(head, anam, resnam)
     ahp_rec *head;		/* Head of ahp_list */
     char *anam;
     char *resnam;
{
  /* Local Variables */
  char *cpos;			/* A multi-purpose character pointer */
  int done;			/* A 'done' flag */
  char atom[PDB_ANAM_SZ+1];	/* The Atom Name of the current atom */
  char residue[PDB_RESNAM_SZ+1]; /* The Resdiue Name - current atom */
  ahp_rec *currec;		/* The current record in the ahplist */

  /* Copy the atom and residue names to local variables, and upcase them */
  /* ------------------------------------------------------------------- */

  /* Atom name */
  strncpy(atom, anam, PDB_ANAM_SZ);
  cpos=atom;
  while(*cpos)
    {
      if (isalpha(*cpos))
	*cpos = toupper(*cpos);
      cpos++;
    }
  atom[PDB_ANAM_SZ] = '\0';

  /* Residue name */
  strncpy(residue, resnam, PDB_RESNAM_SZ);
  cpos=residue;
  while(*cpos)
    {
      if (isalpha(*cpos))
	*cpos = toupper(*cpos);
      cpos++;
    }
  residue[PDB_RESNAM_SZ] = '\0';


  /* Find the hpval for this atom... */
  /* ------------------------------- */
  done = FALSE;
  currec = head;
  while (!done) {
    if (!currec)
      done = TRUE;
    else {
      if ((!strcmp(atom, currec->atom)) || !strcmp("*", currec->atom))
	if ((!strcmp(residue, currec->residue)) ||
	    (!strcmp("*", currec->residue)))
	  done = TRUE;
	else
	  currec = currec->next;
      else
	currec = currec->next;
    }
  }

  /* If none of that caught it, we don't know what the heck it is! */
  if (!currec)
    {
      printf("\nWARNING -- Unknown atom type!\n");
      printf("Can't find AHPhil value for...\n");
      printf("Atom: %s\tResidue: %s\n\n", atom, residue);
      return(-1);
    }

  /* If we are debugging, print the value */
  if (debug_level >= 3)
    printf("HPHIL:\tatom=%s\t\tResidue=%s\tValue=%f\n",
	   atom, residue, currec->ahp);

  /* Return the hpval and exit */
  return(currec->ahp);
  
}

/***************************************************************************/
/* fillActiveTable - This function accepts a file pointer and an integer   */
/* unlocated integer pointer.It then scans the file, find out the water #  */
/* in that file. As the file contains the PDB records of active site waters*/
/* we will have a table of all active site waters. The funtion exits and   */
/* returns the size of the active site table                               */
/***************************************************************************/
int fillActiveTable(FILE *activefile, int **table) {
 int i,numlines;
 char string[MAXLN+1];
 
 /* Init */
 rewind(activefile);
 numlines=0;

 /* Pass 1*/
 fgets(string,MAXLN,activefile);
 while(!feof(activefile)) {
   fgets(string,MAXLN,activefile);
   numlines++;
 }

 /* Pass 2*/
 *table=(int *)malloc(sizeof(int)*numlines);
 rewind(activefile);
 for(i=0;i<numlines;i++) {
   fgets(string,MAXLN,activefile);
   sscanf(string,"%*s%*s%*s%*s%d",(*table+i));
 }
 
 return numlines;
}
 

/***************************************************************************/
/* MAIN ROUTINE                                                            */
/***************************************************************************/
main(int argc, char *argv[])
{
  /* Structures and typedefs */
  typedef struct {
    char anum[PDB_ANUM_SZ+1];	    /* The atom serial number */
    char resnum[PDB_RESNUM_SZ+1];   /* The residue sequence number */
    char anam[PDB_ANAM_SZ+1];       /* The atom name */
    char resnam[PDB_RESNAM_SZ+1];   /* The residue name */
    float x_coord;		    /* The coordinates... */
    float y_coord;
    float z_coord;
    float bval;			    /* The BVal */
  } atom_rec;

  typedef struct {
    char resnum[PDB_RESNUM_SZ+1];   /* The residue sequence number */
    float x_coord;		    /* The coordinates... */
    float y_coord;
    float z_coord;
    float ahphil;		/* The atomic hydrophilicity value */
    int density;                /* Atomic Density of the water molecule */
    int hbondcount_p;           /* Number of hbonds it forms it protein atoms*/
    int hbondcount_w;		/* Number of hbonds with water */
    float occup;                /* Occupancy */
    float bval;                 /* Bond value */
    float c_sum_bval;		/* Sum of bond values of neighbouring nonwater atoms*/
    float avg_bval;		/* Avg of Bval = sum_of_bval/adn or -1 (if adn=0) */
  } water_rec;

  /* Local Variables */
  FILE *protpdb;		/* Pointer to the protein PDB file */
  FILE *watsfile;		/* Pointer to the waters file */
  FILE *consfile;		/* Pointer to the .cons file */
  FILE *activefile;		/* Pointer to the .active.hits file */
  FILE *outfile;		/* Pointer to the output file */
  char prmpt[MAXLN + 1];	/* Buffer for user prompts */
  char fbuf[MAXLN + 1];		/* Buffer for file I/O */
  char ubuf[MAXLN + 1];		/* Buffer for user input */
  char fdefault[MAXLN + 1];	/* Buffer for default filenames */
  char sdescrip[MAXLN + 1];	/* Description string for scanf statements */
  char cfirst_res[PDB_RESNUM_SZ+1]; /* Tag of first residue */
  char clast_res[PDB_RESNUM_SZ+1];  /* Tag of last residue */
  char wats_code[PDB_CODE_SZ+1];    /* PDB code for waters file */
  int *activeSiteTable;		/* Table containing active site water numbers */
  int numActiveSites;		/* Number of active sites */
  int i;			/* A loop counter */
  int num_atoms;		/* Number of atoms in the PDB file */
  int cur_atom;			/* The current atom in the atom_array */
  int ifirst_res;		/* Residue number of the 1st residue */
  int ilast_res;		/* Residue number of last residue */
  int cur_res;			/* The current residue number */
  double r_x;			/* X coordinate difference */
  double r_y;			/* Y coordinate difference */
  double r_z;			/* Z coordinate difference */
  double dist;			/* Distance */
  atom_rec *atom_array;		/* The array of atoms from the PDB file */
  water_rec cur_water;		/* The current water from the wats file */
  ahp_rec *ahp_list;		/* The list of possible AHP values */

  float avg_bval,avg_occup;            /* Variables used for calculating */
  float cur_bval,cur_occup,mobility;   /* mobility */
  int counter=0;                       /* This  variable will be used to count the # of HETATM waters */
  int cons_status;

  int mode;

  /* Initializations */
  ifirst_res = 0;
  ilast_res = 0;


  /* Determine the mode of operation */
  /* if test mode, ask for .cons file and then append it to the end */
  /* if application mode, don't ask for the .cons file */
  switch(argc) {
  case 2: { /* Interactive mode */
    if(argv[1][1]!='t' && argv[1][1]!='a') {
      printf("Usage: %s [-a|-t]\n","all.exe");
      printf("-a    Application mode.\n");
      printf("-t    Test mode.\n");
      exit(-1);
    }
    if(argv[1][1]=='a')
      mode=APP_MODE;
    else
      mode=TEST_MODE;
    /* Print splash banner */
    printf("\nCompute: Atomic Hydrophilicity, Atomic Density, HBond #, Bval and Mobility\n");
    printf("===========================================================================\n\n");
      
    /* Open the two input files */
    /* ======================== */
    
    /* PDB file */
    strcpy(prmpt, "Enter the filename for the protein PDB file...");
    fdefault[0] = '\0';
    protpdb = get_file(prmpt, "r", fdefault,TRUE);
    
    /* Waters file */
    strcpy(prmpt,
	   "\nEnter the filename for the file containing water coordinates...");
    fdefault[0] = '\0';
    watsfile = get_file(prmpt, "r", fdefault,TRUE);
    
    /* .cons file, .active.hits required only in test mode */
    if(mode==TEST_MODE) 
      {/* cons file */
	strcpy(prmpt,"\nEnter the .cons filename");
	fdefault[0] = '\0';
	consfile = get_file(prmpt, "r", fdefault,TRUE);
	
	/* active.hits file */
	strcpy(prmpt,"\nEnter the .active.hits filename");
	fdefault[0] = '\0';
	activefile = get_file(prmpt,"r",fdefault,TRUE);
	numActiveSites=fillActiveTable(activefile,&activeSiteTable);
      }
  }
  break;

  case 5: { /* Non Interactive */
    printf("Running All.EXE...\n");
    if(argv[1][1]!='t' && argv[1][1]!='a') {
      printf("Usage: %s [-a|-t]\n","all.exe");
      printf("-a    Application mode.\n");
      printf("-t    Test mode.\n");
      exit(-1);
    }
    if(argv[1][1]=='a')
      mode=APP_MODE;
    else
      mode=TEST_MODE;
    
    
    protpdb = get_file("","r",argv[2],FALSE);
    watsfile = get_file("","r",argv[3],FALSE);
    if(protpdb==NULL){
        printf("Error opening %s. Exiting...\n",argv[2]);
	exit(-1);
    }
    if(watsfile==NULL){
	printf("Error opening %s. Exiting...\n",argv[3]);
	exit(-1);
    } 
    /* If running in test mode...*/
    if(mode==TEST_MODE) {
      strcpy(fdefault,argv[4]);
      strcat(fdefault,".cons");
      consfile = get_file("", "r", fdefault,FALSE);
	
      /* active.hits file */
      strcpy(fdefault,argv[4]);
      strcat(fdefault,".active.hits");
      activefile = get_file("","r",fdefault,FALSE);
      numActiveSites=fillActiveTable(activefile,&activeSiteTable);
    }
  }
  break;
  
  default: {/* Print Usage */
      printf("Usage: %s [-a|-t] pdbProteinFile pdbWatsFile PDBCode\n",argv[0]);
      printf("-a    Application mode.\n");
      printf("-t    Test mode.\n");
      exit(-1);
  }
  }
    
    /* ================================================================ */
  /* Step 1:  Load the coordinates of every atom in the PDB file into */
  /*          memory for fast access.                                 */
  /* ================================================================ */

  /* Load the list of AHP values */
  load_list(&ahp_list);
    

  /* Pass 1: Determine the amount of memory necessary... */
  /* --------------------------------------------------- */
  printf("\nExamining protein PDB file...\n\n");
  num_atoms = 0;
  while(!feof(protpdb))
    {
      if (fgets(fbuf, MAXLN, protpdb))
	if ((!memcmp(fbuf,"ATOM",4)) || (!memcmp(fbuf,"HETATM",6)))
	  num_atoms++;
    }

  /* Allocate the space for the atom array */
  printf("%d atoms found in PDB file.\n", num_atoms);
  if (debug_level >= 2)
    printf("Allocating %d bytes for atom array...\n",
	   num_atoms*sizeof(atom_rec));
  atom_array = (atom_rec *) malloc(sizeof(atom_rec)*num_atoms);

  /* Pass 2: Load the atom array */
  /* --------------------------- */
  cur_atom = 0;
  avg_bval=0;avg_occup=0;
  rewind(protpdb);
  while((!feof(protpdb)) && cur_atom < num_atoms)
    {

      /* For each line found, check if it is an atom... */
      if (fgets(fbuf, MAXLN, protpdb))
	if ((!memcmp(fbuf,"ATOM",4)) || (!memcmp(fbuf,"HETATM",6)))
	  
	  /* ...if so, then capture it to the atom array */ 
	  {  /****  ****/

	    /* The sprintf's in the following block set up description     */
	    /* strings which allow the sscanf's to capture the appropriate */
	    /* block of data from the PDB file.  Although it looks a bit   */
	    /* confusing, this method allows the PDB column positions      */
	    /* and sizes to be stored in defined constants.  This makes it */
	    /* much easier to fix the source code if the PDB file format   */
	    /* should change or (more likely) a bug is found in the code.  */

	    sprintf(sdescrip, "%%*%dc%%%ds", PDB_ANUM_COL-1, PDB_ANUM_SZ);
	    sscanf(fbuf, sdescrip, &(atom_array[cur_atom].anum));
	    sprintf(sdescrip, "%%*%dc%%%ds", PDB_RESNUM_COL-1, PDB_RESNUM_SZ);
	    sscanf(fbuf, sdescrip, &(atom_array[cur_atom].resnum));
	    sprintf(sdescrip, "%%*%dc%%%ds", PDB_ANAM_COL-1, PDB_ANAM_SZ);
	    sscanf(fbuf, sdescrip, &(atom_array[cur_atom].anam));
	    sprintf(sdescrip, "%%*%dc%%%ds", PDB_RESNAM_COL-1, PDB_RESNAM_SZ);
	    sscanf(fbuf, sdescrip, &(atom_array[cur_atom].resnam));
	    sprintf(sdescrip, "%%*%dc%%f %%f %%f %%f", PDB_XCRD_COL-1);
	    sscanf(fbuf, sdescrip, &(atom_array[cur_atom].x_coord),
		   &(atom_array[cur_atom].y_coord),
		   &(atom_array[cur_atom].z_coord),
		   &cur_occup);
	    sprintf(sdescrip,"%%*%dc%%f",PDB_BVAL_COL-1);
	    sscanf(fbuf,sdescrip,&(atom_array[cur_atom].bval));

	    cur_bval=atom_array[cur_atom].bval;

	     /* checking if the record read is a water. If it is, add the  mobility and bval to */
	     /* the accumulating sum */
	     if( (strncmp(fbuf,"HETATM",6)==NULL) )   /* &&& */
	      {
	       char temp[4];
	     
	       /* copying the current residue name to temp. This helps in simpler, clear coding */
	       strncpy(temp,atom_array[cur_atom].resnam,3);
	       temp[3]=NULL;
	       
	       /* If the HETATM is one of the waters ...*/
	       if( (!strcmp(temp,"HOH")) || (!strcmp(temp,"H2O")) || (!strcmp(temp,"WAT"))
		   || (!strcmp(temp,"DOD")) || (!strcmp(temp,"D2O")) )
		 { avg_bval+=cur_bval;   /* Accumulating the Bval and */
	           avg_occup+=cur_occup; /* occupancy for finding averages */
		   counter++;            /*  Increment  counter which will be used for averaging */
		 } /* end of if*/
	      
	      } /* end of 'if' matching &&& */
	    



	     /* Debug dump */
	    if (debug_level >= 3)
	      {
		printf("Loaded:\n");
		printf("\tanum:\t%s\n", atom_array[cur_atom].anum);
		printf("\tresnum:\t%s\n", atom_array[cur_atom].resnum);
		printf("\tanam:\t%s\n", atom_array[cur_atom].anam);
		printf("\tresnam:\t%s\n", atom_array[cur_atom].resnam);
		printf("\tcoords:\t%f\t%f\t%f\n\n",
		       atom_array[cur_atom].x_coord,
		       atom_array[cur_atom].y_coord,
		       atom_array[cur_atom].z_coord);
	      }

	    /* Note if this is the first or last residue */
	    if (!memcmp(fbuf,"ATOM",4))
	      {
		sscanf(atom_array[cur_atom].resnum, "%d", &cur_res);
		if ((ifirst_res == 0) || (cur_res < ifirst_res))
		  {
		    ifirst_res = cur_res;
		    strncpy(cfirst_res, atom_array[cur_atom].resnum,
			    PDB_RESNUM_SZ);
		  }
		if ((ilast_res == 0) || (cur_res > ilast_res))
		  {
		    ilast_res = cur_res;
		    strncpy(clast_res, atom_array[cur_atom].resnum,
			    PDB_RESNUM_SZ);
		  }
	      }

	    cur_atom++;
	  } /**** ****/
    }/* end of while */

   avg_bval=avg_bval/counter;
   avg_occup=avg_occup/counter;
   printf("\nAvg Bval: %8.2f  Avg Occup: %8.2f\n",avg_bval,avg_occup);

  /* Debug output */
  if (debug_level >= 2)
    {
      printf("\nLowest residue tag found:\t%s\n", cfirst_res);
      printf("(Numeric value:\t%d)\n", ifirst_res);
      printf("Highest residue tag found:\t%s\n", clast_res);
      printf("(Numeric value:\t%d)\n\n", ilast_res);
    }

  /* ================================================================ */
  /* Step 2:  For each water atom in the waters file, check which     */
  /*          protein atoms are within HY_CUTOFF angstroms and        */
  /*          compute the atomic hydrophilicity for each hit.         */
  /*          The sum of the hits is the value for the water.         */
  /* ================================================================ */
  printf("Examining waters...\n");
  
  /* Find the current PDB code for the protein file */
  rewind(protpdb);
  fgets(fbuf, MAXLN, protpdb);
  sprintf(sdescrip, "%%*%dc%%%ds", PDB_CODE_COL-1, PDB_CODE_SZ);
  sscanf(fbuf, sdescrip, wats_code);

  if(argc==2) {
    /* Prompt the user for any changes */
    printf("\nPlease enter the PDB code for the waters file,\n");
    printf("or press Return to use [%s]...\n", wats_code);
    fgets(ubuf, MAXLN-1, stdin);
    
    /* Strip the CR from the response */
    for (i=0; i<strlen(ubuf); i++)
      if (ubuf[i] == '\n')
	ubuf[i] = '\0';
  }
  else {
    /* Non interactive */
    strcpy(ubuf,argv[4]);
  }
    
  /* If there is anything left, then make that the current PDB code */
    if (strlen(ubuf) >= 1)
      for(i=0; ((i<=strlen(ubuf)) && (i < PDB_CODE_SZ)); i++)
	if (isalpha(ubuf[i]))
	  wats_code[i] = toupper(ubuf[i]);
	else
	  wats_code[i] = ubuf[i];
    wats_code[PDB_CODE_SZ] = '\0';
    
  /* Debug printout */
  if (debug_level >= 1)
    printf("Using PDB code \"%s\" for waters.\n\n", wats_code); 

  /* get ready to open the output file */
  strcpy(prmpt, "Enter the filename for the output file...");
  strcpy(fdefault, wats_code);
  strcat(fdefault, ".env");

  /* ... but first put the default filename in lower case... */
  for(i=0; i<strlen(fdefault); i++)
    if (isalpha(fdefault[i]))
      fdefault[i] = tolower(fdefault[i]);

  /* .. Ok, now open it! */
  if(argc==2) { /* Go interactive */
    do{
      outfile = get_file(prmpt, "w", fdefault,TRUE);
    }while(outfile==NULL);
    printf("\n");
  }
  else {
    outfile = get_file("","w",fdefault,FALSE);
    if(outfile==NULL) {
	printf("Error opening %s. Exiting...\n",fdefault);
	exit(-1);
    }
  }
  /* Now examine each water... */
  /* ------------------------- */

  printf("Computing data ...\n");
  rewind(watsfile);
 if(mode==TEST_MODE)
   {
    fprintf(outfile,"#%6s %6s %6s     %4s   %4s      %5s  %5s   %s  %s %s %s %s\n",
    "PBDCode","Resi #","Adn","Ahp","Bval","hbd#P","hbd#W","Mob","Net BVal","Avg Bval","Consrvd","Active");
    fprintf(outfile,"###############################################################################################\n");
   }
 else 
   {
     fprintf(outfile,"#%6s %6s %6s     %4s   %4s      %5s  %5s  %s  %s %s\n",
    "PBDCode","Resi #","Adn","Ahp","Bval","hbd#P","hbd#W","Mob","Net BVal","Avg Bval");
     fprintf(outfile,"##############################################################################\n");
   }
 
  while(!feof(watsfile))
    {
      if(fgets(fbuf, MAXLN, watsfile))
        {
          /* Initialize the structure */
          cur_water.ahphil = 0;

          /* Parse the current water record */
          sprintf(sdescrip, "%%*%dc%%%ds", PDB_RESNUM_COL-1, PDB_RESNUM_SZ);
          sscanf(fbuf, sdescrip, &(cur_water.resnum));
          sprintf(sdescrip, "%%*%dc%%f %%f %%f %%f", PDB_XCRD_COL-1);
          sscanf(fbuf, sdescrip, &(cur_water.x_coord),
                 &(cur_water.y_coord),
                 &(cur_water.z_coord),
		 &(cur_water.occup));
	  sprintf(sdescrip, "%%*%dc%%f",PDB_BVAL_COL-1);
	  sscanf(fbuf,sdescrip,&(cur_water.bval));
          
          /* Debug printout */
          if (debug_level >= 3)
            {
              printf("\nExamining water:\n");
              printf("\tResidue Number:\t%s\n", cur_water.resnum);
              printf("\tCoordinates:\t%f\t%f\t%f\n", cur_water.x_coord,
                     cur_water.y_coord, cur_water.z_coord);
            }

          /* ----------------------------------------------------- */
          /* Check this water against every atom in the atom array */
          /* ----------------------------------------------------- */

	  /* Initialize the values which are computed */
	   cur_water.ahphil       = 0;
           cur_water.density      = 0;
	   cur_water.hbondcount_p = 0;
	   cur_water.hbondcount_w = 0;
	   cur_water.c_sum_bval   = 0;
	   cur_water.avg_bval     = 0;
          for(cur_atom=0; cur_atom < num_atoms; cur_atom++)
            {
              /* Check if this is a neighbor */
              r_x = cur_water.x_coord - atom_array[cur_atom].x_coord;
              r_y = cur_water.y_coord - atom_array[cur_atom].y_coord;
              r_z = cur_water.z_coord - atom_array[cur_atom].z_coord;
              dist = sqrt((r_x*r_x)+(r_y*r_y)+(r_z*r_z));
              
              /* If so, then compute its atomic hydrophilicity and */
              /* add it to the running total.                      */
              
              if ((dist > 0) && (dist <= HY_CUTOFF))
                {
		  float cur_ahp;

		  /* Calculate AHP: START */
 		  cur_ahp=hphil(ahp_list, atom_array[cur_atom].anam,
                                          atom_array[cur_atom].resnam);
		  /* If the atom and residue not found in the Ahplist then */
		  /* Use the first char from the aname as atom name and residue */
		  /* as AVG. This is the final place where all atoms are catched */
		  /* AVG entries are put at the end */
                  if(cur_ahp<0)
		      {
			correctName(atom_array[cur_atom].anam,atom_array[cur_atom].resnam);
   		        cur_ahp = hphil(ahp_list, atom_array[cur_atom].anam,atom_array[cur_atom].resnam);
			printf("Using average AHP value of %f for: %s\n",cur_ahp,
			       atom_array[cur_atom].anam);
		      }
                  cur_water.ahphil += cur_ahp;
		  /* Calculate AHP: END */

		  /* Calculate ADN: START */
                  cur_water.density+=adn(atom_array[cur_atom].anam,
					 atom_array[cur_atom].resnam);
		  /* Calculate ADN: END */
 
	          /* Check for Bval of all nonwaters near the water */
		  if(iswater(atom_array[cur_atom].resnam)==FALSE)
  	            cur_water.c_sum_bval+=atom_array[cur_atom].bval;
                }
	      /* This "switch" for hbond part*/
		switch( hbondpossible(atom_array[cur_atom].anam,
				      atom_array[cur_atom].resnam,
				      dist) ) {
		case WITH_PROTEIN:
		  cur_water.hbondcount_p++;
		  break;
		case WITH_WATER:
		  cur_water.hbondcount_w++;
		default:
		   break;
		}
	     /* End of Hbond calculation */

             /* Calculation of avg_bval */
	      cur_water.avg_bval= (cur_water.density==0)?-1:(cur_water.c_sum_bval/cur_water.density);
            
	  }/**** end of 1 atom in .xfm file per atom in .hits file  and end of "for" loop **/
	  
          /* Debug printout */
          if (debug_level >= 3)
            printf("%4s  %6s   %8.4f %d\n", wats_code, cur_water.resnum, 
                 cur_water.ahphil,cur_water.density);
          
	
	  /* Calculating Mobility */
	  mobility=cur_water.bval*avg_occup/avg_bval/cur_water.occup;


          /* Save the results for this water */
          if(mode==TEST_MODE)
	    { int isactive=0,fini=0,x,n_resnum;
	     /* Load the conserved/displaced info from the .cons file */
	      fgets(fbuf,MAXLN,consfile);
	      sscanf(fbuf,"%*s%*s%d\n",&cons_status);
                  
             /* find whether its active site or not */
	      n_resnum=atoi(cur_water.resnum);
	      for(x=0;x<numActiveSites && !fini;x++)
		if(activeSiteTable[x]==n_resnum)
		  {
		    fini=1;
		    isactive=1;
		  }

	     /* now to print the info */
             fprintf(outfile, "%6s %6s    %4d  %8.4f  %8.2f   %4d %4d  %7.4f %7.2f   %8.2f    %d   %d\n", 
              wats_code, cur_water.resnum,cur_water.density,cur_water.ahphil,cur_water.bval,
	      cur_water.hbondcount_p,cur_water.hbondcount_w,mobility,cur_water.c_sum_bval,
	      cur_water.avg_bval, cons_status,isactive);
	    }
          else
           fprintf(outfile, "%6s %6s    %4d  %8.4f  %8.2f   %4d %4d %7.4f %7.2f  %8.2f\n", 
            wats_code, cur_water.resnum,cur_water.density,cur_water.ahphil,cur_water.bval,
	    cur_water.hbondcount_p,cur_water.hbondcount_w,mobility,cur_water.c_sum_bval,
	    cur_water.avg_bval);
  
    
        }
    }
  printf("...done!\n");
  
  /* Clean up and exit */
  /* ----------------- */
  free(atom_array);
  free_list(ahp_list);
  fclose(protpdb);
  fclose(watsfile);
  if(mode==TEST_MODE) {
    fclose(consfile);
    fclose(activefile);
  }
  fclose(outfile);
  printf("\n");
}






