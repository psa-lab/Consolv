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
/*                                                                         */
/*                ==========================================               */
/*                !!!!!!! READ THIS BEFORE COMPILING !!!!!!!               */
/*                ==========================================               */
/*                                                                         */
/*  This code was origninally written to use only those pdb files where    */
/*  the inhibitor is a protein or peptide.  Thus, it expects that there    */
/*  will be two or more chains in the PDB file, and the protein and        */
/*  inhibitor will each have a distinct chain ID.  Since this code uses    */
/*  the chain ID for loading, it currently cannot handle pdb files where   */
/*  the inhibitor does not have a chain ID.  (For example, where the       */
/*  inhibitor is composed entirely of HETATM's.)                           */
/*                                                                         */
/*  However, it should not be difficult to include this capability         */
/*  by looking at HETATM's during the loading phase.                       */
/*                                                                         */
/*  All in all, this is a quick hack to get some results quickly.  Before  */
/*  further capabilities are added, it would probably be a good idea to    */
/*  go back and break the code into modules, and move some numeric         */
/*  constants into define statements, etc.                                 */
/*                                                                         */
/***************************************************************************/

/***************************************************************************/
/* nearby  -- Find coordinates of water molecules bound to or displaced    */
/*            by the inhibitor.                                            */
/*                                                                         */
/* This simple program is designed to look at two pdb files.  The first    */
/* should contain a protein:inhibitor complex, while the second file       */
/* should have either the protein or the inhibitor alone.  The ATOM        */
/* records of the two PDB files should be in the same coordinate system.   */
/*                                                                         */
/* The program is designed to dump the PDB records for all waters in the   */
/* active site of each protein.  A water is defined to be in the active    */
/* site if it is either:                                                   */
/*  a) within HBOND_CUTOFF angstroms of both the protein and the inhibitor */
/*     (for the "free" protein - the water must be within HBOND_CUTOFF     */
/*      angstroms of the protein and where the inhibitor binds in the      */
/*      the complex form.  This is why the two must be in the same         */
/*      coordinate system.)                                                */
/*  - or -                                                                 */
/*  b) if 1. It is withen HBOND_CUTOFF of a protein AND withen HBOND_CUTOFF*/
/*           of a water which is withen HBOND_CUTOFF of the ligend         */
/*     OR                                                                  */
/*        2. It is withen HBOND_CUTOFF of a ligend AND withen HBOND_CUTOFF */
/*           of a water which is withen HBOND_CUTOFF of the protein.       */
/*                                                                         */
/*     OR 2. It is withen HBOND_CUTOFF of a monowater withen HBOND_CUTOFF  */
/*           of protei or ligand but not both                              */
/* Hacker - Michael Lee Raymer, August, 1994.                              */
/* Hacker II - Vishal C. Thakkar, Dec, 1996.                               */
/*              Modified the code to include the diwater bridges as defined*/
/*              in part b) above.                                          */
/*									   */
/* HBOND_CUTOFF		- revised to 3.5 instead of the original 3.6 (SN)  */
/*									   */
/***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <sys/stat.h>

/* Here are some global constants, if this grows into a mult-module */
/* application, then these will be moved to a global.h file         */
#define HBOND_CUTOFF 3.5           /* "Close" water distance cuttoff */
#define MAXLN 132                  /* Max input line length */
#define MAXCHAINS 10               /* Q&D:  Max number of chains */
#define ID_COL 21                  /* PDB column for chain ID */
#define BLANK ' '                  /* Blank character */
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

/***************************************************************************/
/* SUBROUTINES AND FUNCTIONS                                               */
/***************************************************************************/

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
/*          interactive = if 1, then will prompt for name and overwriting  */
/*                        if 0, will not prompt, and will rename before    */
/*                        overwriting the file                             */
/*                                                                         */
/* OUTPUTS:  get_file returns a pointer to the open file.                  */
/***************************************************************************/
FILE *get_file(prmpt, omode,fdefault, interactive)
char *prmpt;                       /* String to prompt user with */
char *omode;                       /* Mode to open the file in */
char *fdefault;			   /* The default file name, if any */
int interactive;		   /* inteactive=1 means ask, else move file instead of overwriting */
{
  /* Local Variables */
  char input_fn[MAXLN + 1];        /* User inputted filename */
  FILE *fptr;                      /* Pointer to opened file */

  input_fn[0] = '\0';
  fptr = (FILE *) NULL;
  while (!fptr)
    {
      if(interactive==FALSE) {
	strcpy(input_fn,fdefault);
      }
      else
	/* Get a valid filename */
	while (strlen(input_fn) <= 1)
	  {
	    printf("%s\n", prmpt);
	    fgets(input_fn, MAXLN - 1, stdin);
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

	  if(stat(input_fn,&b)==NULL) {/* File exits */
	    if(interactive==TRUE)
	    {/* Ask user */
	      printf("file: %s exits, Overwrite [y/n] ?",input_fn);
	      scanf("%c",&ans);
	      fflush(stdin);  /*Flush the input buffer of any carriage returns */
	      ans=toupper(ans);
	      if(ans=='N')
		return(NULL); /* Return Null pointer */
	    }
	    else {/* don't ask the user, just rename the file */
	      strcpy(bakname,input_fn);
	      strcat(bakname,".bak");
	      printf("\nwarning: moving file %s to %s. \n",input_fn,bakname);
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
	  if(interactive==FALSE) return(0);
        }
    }
  return(fptr);
}


/***************************************************************************/
/* MAIN ROUTINE                                                            */
/***************************************************************************/

/* Type defines */
typedef struct {
  char pdbrec[MAXLN + 1];
  float x_coord;
  float y_coord;
  float z_coord;
  enum bond_type { direct, diwater, none } bond;        /* Type of bond */
  enum isclose_type { yes, no } isclose2p,isclose2i;     /* close to protein or inhibitor*/
  enum isclose2what_type 
  { protein,inhibitor,
    both, neither
  } isclose2watclose2;          /* Close to water which is close to p,i,both,neither */
} water_rec;

typedef struct {
  char struct_type[10];
  water_rec *waters;		/* Pointer to waters array */
  FILE *pdb;			/* Pointer to PDB file */
  FILE *outfile;		/* Pointer to output file. The file is closed inside the sub */
  int w_found;			/* # of waters in that file to process */
  int pchain_len;		/* Length of the Protein chain */
  int ichain_len;		/* Length of the Inhibitor chain */
  float *pchain_coord;		/* Pointer to the protein chain coords */
  float *ichain_coord;		/* Pointer to the inhibitor chain coords */
} process_data;

/**********************************************************************************************/
/* This procedure takes a data structure of type process_data and does the following          */
/* 1. read the HETATOM records from the pdb file (either free or complex and stores them in   */
/*    the waters array. The waters array of an array of water_rec structures.                 */
/* 2. While reading the HETATOM records, sets the enum data types to default values           */
/* 3. Pass 1.                                                                                 */
/*    It checks whether the waters loaded are close the inhibitor chain. If so, sets the      */
/*      isclose2i=yes.                                                                        */
/*    It checks whether the waters loaded are close the protein chain. If so, sets the        */
/*      isclose2p=yes. If the water is close to inhibitor AND close to protein its a monowter */
/*      and hence dumps it.                                                                   */
/* 4. Pass 2.                                                                                 */
/*    It checks whether the waters are close to any other waters. If so, it determines which  */
/*       water is it close to, ie is it close to a water which is close to inhibitor,protein, */
/*       both or neither                                                                      */
/*    Then it checks for diwater sites. A site is diwater if, It is close to protein and close*/
/*      to a water which is close to inhibitor but not the protein OR It is close to inhibitor*/ 
/*      and close  to a water which is close to protein but not the inhibitor.                */
/* 4. Closes the output file                                                                  */
/* ********************************************************************************************/
void process( process_data p )
{
  char prmpt[100];
  char fbuf[MAXLN+1];
  char het_type[MAXLN+1];
  float het_x, het_y, het_z, dist;
  int curwat,w_dumped,found;
  int i;

  if(p.outfile==NULL) {/* if p.outfile==NULL, we have not opened the file, so have to open it*/
                       /* if p.outfile!=NULL, we have opened the file, so have to close later */
    sprintf(prmpt, "Enter filename for %s structure active site output file:",p.struct_type);
    do{
      p.outfile = get_file(prmpt, "w","",TRUE);
    }while(p.outfile==NULL);
  }
    /* Scan the complex file for waters */
  printf("\nScanning %s file for waters:\n",p.struct_type);
  rewind(p.pdb);
  curwat = 0;
  while (!feof(p.pdb))
    {
      /* Find the next HETATM */
      fbuf[0] = '\0';
      while((memcmp(fbuf, "HETATM", 4)) && !feof(p.pdb))
        fgets(fbuf, MAXLN, p.pdb);

      /* Get the coordinates and the type */
      if (!feof(p.pdb))
        {
          sscanf(fbuf, "%*16c %s %*7c %f %f %f",
                 het_type, &het_x, &het_y, &het_z);

          /* If it is a water, store it */
          /*if (!memcmp(het_type, "HOH", 3))*/
	  if(iswater(het_type))
            {
              memcpy(p.waters[curwat].pdbrec, fbuf, MAXLN);
              p.waters[curwat].x_coord = het_x;
              p.waters[curwat].y_coord = het_y;
              p.waters[curwat].z_coord = het_z;
              p.waters[curwat].bond = none;
              p.waters[curwat].isclose2p=no;
              p.waters[curwat].isclose2i=no;
              p.waters[curwat].isclose2watclose2=neither;          
              curwat++;
            }
        }
    }
  
  /* dump pass 1 --                                  */
  /* If a water is within HBOND_CUTOFF angstroms of an atom in the */
  /* inhibitor and an atom in the protein then dump it...          */
  
  printf("%s file:\t\t%d waters found\n",p.struct_type, p.w_found);
  w_dumped = 0;
  for (curwat=0; curwat < p.w_found; curwat++)
    {
      /* check the inhibitor */
      found = FALSE;
      het_x = p.waters[curwat].x_coord;
      het_y = p.waters[curwat].y_coord;
      het_z = p.waters[curwat].z_coord;
      for (i=0; (i<p.ichain_len*3) && (!found); i+=3)
        {
          if(curwat!=i)
          {   dist = hypot(hypot(p.ichain_coord[i] - het_x, p.ichain_coord[i+1] - het_y),
                         p.ichain_coord[i+2] - het_z);
              if (dist < HBOND_CUTOFF)
                {
                  found = TRUE;
                  p.waters[curwat].isclose2i=yes;
                }
          }
        }

      /* check the protein */
       found = FALSE;
       for (i=0; (i<p.pchain_len*3) && (found==FALSE); i+=3)
          {/* 2 */
            dist = hypot(hypot(p.pchain_coord[i] - het_x,
                               p.pchain_coord[i+1] - het_y),
                               p.pchain_coord[i+2] - het_z);
            if (dist < HBOND_CUTOFF)
            {/* 1 */
                p.waters[curwat].isclose2p=yes;
                /* This water is within HBOND_CUTOFF of protein.*/
                /* Is it withen HBOND_CUTOFF of inhibitor? */
                if( p.waters[curwat].isclose2i==yes)
                  { /* If so, it forms a monowater bridge! Hence dump it! */
                     w_dumped++;
                     fprintf(p.outfile, "%s", p.waters[curwat].pdbrec);
                     p.waters[curwat].bond = direct;         /* This bond is direct ie monowater  */
                     found = TRUE;
                  }  
             }/* 1 */
          }/* 2 */
        
      }/* looping through all the water records over!*/ 

  printf("Finished pass 1:\t%d %s active site mono-waters dumped\n", w_dumped,p.struct_type);
  
  /* Dump Pass 2:  Check if any of the remaining waters are */
  /* h-bonded to one of the waters dumped in pass 1.        */
  /* ------------------------------------------------------ */

  /* Checking whether any water is withen HBOND_CUTOFF of any other  water */
  for ( curwat=0; curwat < p.w_found; curwat++ )
    {
      if(p.waters[curwat].bond == none)
        {
          found = FALSE;
          het_x = p.waters[curwat].x_coord;
          het_y = p.waters[curwat].y_coord;
          het_z = p.waters[curwat].z_coord;
          for (i=0; ((i < p.w_found) && (!found)); i++)
            {
              if (i != curwat)
                {
                  dist = hypot( hypot(p.waters[i].x_coord - het_x,p.waters[i].y_coord - het_y),
                               p.waters[i].z_coord - het_z);
                  if (dist < HBOND_CUTOFF)
                    {
                      found = TRUE;
                      if( (p.waters[i].isclose2p==yes) && (p.waters[i].isclose2i==yes) )
                        p.waters[curwat].isclose2watclose2=both;
                      if( (p.waters[i].isclose2p==yes) && (p.waters[i].isclose2i==no) )   
                        p.waters[curwat].isclose2watclose2=protein;
                      if( (p.waters[i].isclose2p==no) && (p.waters[i].isclose2i==yes) )   
                        p.waters[curwat].isclose2watclose2=inhibitor;
                      if( (p.waters[i].isclose2p==no) && (p.waters[i].isclose2i==no) )    
                        p.waters[curwat].isclose2watclose2=neither;
                    }
                }
             }
        }
    }
/* At this point, we have classifed all water atoms to be withen HBOND_CUTOFF another water */
/* Now to process all of the waters in the complex form to check for diwater bridges */

/* Making w_dumped=0 as now we are counting Di-water bridges */
w_dumped=0;
 for ( curwat=0; curwat < p.w_found; curwat++)
   { water_rec w;
   
     w=p.waters[curwat];
     if( 
	 (
	  (w.isclose2p==yes) && (w.isclose2i==no ) && 
	                  ( (w.isclose2watclose2==inhibitor)||(w.isclose2watclose2==both) )
	  ) 
	 || 
         (
	  (w.isclose2p==no ) && (w.isclose2i==yes) && 
	                  ( (w.isclose2watclose2==protein)||(w.isclose2watclose2==both) ) 
	  )
       )
       { /* satisfies diwater bridge conditions */
         w_dumped++;                                       /* increment dump counter */
         fprintf(p.outfile, "%s", p.waters[curwat].pdbrec);/* write record to file */
	 p.waters[curwat].bond=diwater;
       } 
   }
 
 printf("Finished pass 2:\t%d %s active site di-waters dumped\n\n", w_dumped,p.struct_type);
 fclose(p.outfile);
}
/********************************************************************************************/



int main(int argc, char *argv[])
{

  /* Local Variables */
  int i;                           /* Loop counters, etc. */
  char prmpt[MAXLN + 1];           /* Buffer for user prompts */
  char fbuf[MAXLN + 1];            /* Buffer for input file lines */
  char ibuf[MAXLN + 1];            /* Buffer for user input */
  FILE *cplxpdb;                   /* Pointer to complex PDB file */
  FILE *freepdb;                   /* Pointer to "free" PDB file  */
  long first_atom[MAXCHAINS+1];    /* File position of 1st ATOM recrd */
  long curfilepos;                 /* Current file position */
  int cplx_chains;                 /* Number of chains in the complex */
  int labeled_chains;              /* Number of non-blank labeled chains */
  char chain_id[MAXCHAINS];        /* Chain ID of each chain */
  int chain_len[MAXCHAINS];        /* Number of atoms in each chain */
  int ichain_num;                  /* The number of the inhibitor chain */
  int pchain_num;                  /* The number of the protein chain */
  int cur_chain_num;               /* The number of the current chain */
  float *ichain_coord;             /* The coordinate array for the inhibitor */
  float *pchain_coord;             /* The coordinate array for the protein */
  char het_type[MAXLN];            /* Type of the current HETATM */
  int cw_found, fw_found;          /* Number of waters in complex, free */
  water_rec *cplx_waters;          /* Waters in the complex form */
  water_rec *free_waters;          /* Waters in the free form */

  process_data p;


  switch(argc) {
  case 1:{/* interactive mode */
    /* Print splash banner */
    printf("\nnearby - find waters bound to or displaced by the inhibitor\n");
    printf("===========================================================\n\n");
    
    /* Open the two PDB files... */
    /* ========================= */
    
    /* Complex file */
    strcpy(prmpt, "Enter Complex-structure filename:");
    cplxpdb = get_file(prmpt, "r","",TRUE);
    
    /* Free file */
    strcpy(prmpt,"Enter Free-structure filename:");
    freepdb = get_file(prmpt,"r","",TRUE);
  }
  break;
  
  case 7:{ /* non-interactive mode */
    cplxpdb = fopen(argv[1],"r");
    if(cplxpdb==NULL) {printf("Error opening %s. Exiting...\n",argv[1]);exit(-1);}
    freepdb = fopen(argv[2],"r");
    if(freepdb==NULL) {printf("Error opening %s. Exiting...\n",argv[2]);exit(-1);}
  }
  break;
  default: { /* Print usage */
    printf("Usage: %s ComplexPDBfile FreePDBfile ProteinChainID InhibitorChainID CplxPDBCode FreePDBCode\n",argv[0]);
    exit(-1);
  }
  }/* end of switch */

 
    /* ========================================================== */
    /* Step 1: Load the coordinates of each atom in the inhibitor */
    /*         and in the protein.                                */
    /* ========================================================== */

    /* Pass 1 - Scan the file, determine the sizes for the       */
    /*          various arrays, assemble the lists of chain      */
    /*          id's chain lengths, starting points, etc.        */
    /* --------------------------------------------------------- */
    printf("\nExamining complex PDB file...\n\n");
    cplx_chains = 0;
    cw_found = 0;
    fbuf[0] = '\0';
    while (!feof(cplxpdb))
    {
      /* Find the next ATOM or HETATM record */
      while((memcmp(fbuf, "ATOM", 4)) && (memcmp(fbuf, "HETATM", 6))
            && (!feof(cplxpdb)))
        {
          curfilepos = ftell(cplxpdb);
          fgets(fbuf, MAXLN, cplxpdb);
        }
      
      /* If we found a record, check if it is a new chain */
      if ((!memcmp(fbuf, "ATOM", 4)) || (!memcmp(fbuf, "HETATM", 6)))
        {
          /* If we know about this chain, increment the chain */
          /* size counter...                                  */
          if ((cplx_chains > 0) &&
              (memchr(chain_id, fbuf[ID_COL], cplx_chains)))
            {
              /* Find which chain number it is */
              cur_chain_num = 0;
              while ((cur_chain_num < cplx_chains) &&
                     (chain_id[cur_chain_num] != fbuf[ID_COL]))
                cur_chain_num++;

              /* Increment the counter for that chain */
              chain_len[cur_chain_num]++;
            }

          /* If this is a new chain, add its id to the chain list */
          else
            {
              chain_id[cplx_chains] = fbuf[ID_COL];

              /* Initialize the length to one */
              chain_len[cplx_chains] = 1;

              /* Note the location of the first record */
              first_atom[cplx_chains] = curfilepos;

              /* Increment the chains counter */
              cplx_chains++;
            }

          /* If this is a water, increment the water count */
          if (!memcmp(fbuf, "HETATM", 6))
            {
              sscanf(fbuf, "%*16c %s", het_type);
              /*if (!memcmp(het_type, "HOH", 3))*/
	      if(iswater(het_type))
                cw_found++;
            }

          /* Get the next line */
          curfilepos = ftell(cplxpdb);
          fgets(fbuf, MAXLN, cplxpdb);
        }
    }

  /* Pass 2 - Load the coordinate arrays */
  /* ----------------------------------- */

  /* Display the chain IDs found */
  labeled_chains = 0;
  for (i=0; i < cplx_chains; i++)
    if (chain_id[i] != BLANK) labeled_chains++;
  printf("The complex PDB file appears to have %d labeled chains.\n",
         labeled_chains);
  for (i=0; i < cplx_chains; i++)
    if (chain_id[i] != BLANK)
      printf("Chain: %d\tIdentifier: %c\n", i+1, chain_id[i]);
  printf("\n");

  
  /*If interactive ask for chain Id */
  if(argc==1) {
    /* Prompt for protein chain */
    ibuf[0] = '\0';
    while (!memchr(chain_id, ibuf[0], cplx_chains))
      {
	if (ibuf[0] != '\0') printf("\nNot a valid chain identifier!\n");
	printf("Please enter the chain identifier for the protein...\n");
	fgets(ibuf, MAXLN - 1, stdin);
      }
  }
  else{/* Non interactive */
    strcpy(ibuf,argv[3]);
    printf("Using %s as chain ID for the protien chain\n",ibuf);
  }
  
  /* Find which chain number matches the chain id given */
  pchain_num = 0;
  while ((pchain_num < cplx_chains) && (chain_id[pchain_num] != ibuf[0]))
    pchain_num++;

  /* Load the coordinates for all atoms in the protein */
  /* ================================================= */

  /* Move to the first ATOM record */
  fseek(cplxpdb, first_atom[pchain_num], SEEK_SET);

  /* Allocate space for coordinate array */
  pchain_coord = (float *) malloc(sizeof(float)*chain_len[pchain_num]*3);

  /* Load the coordinates */
  i=0;
  while (!feof(cplxpdb))
    {
      /* If we find an ATOM or HETATM with  the correct chain ID... */
      /* then load it...                                            */
      if ((fgets(fbuf, MAXLN, cplxpdb)) &&
          (fbuf[ID_COL] == chain_id[pchain_num]) &&
          ((!memcmp(fbuf, "ATOM", 4)) || (!memcmp(fbuf, "HETATM", 6))))
        {
          sscanf(fbuf, "%*31c %f %f %f", &(pchain_coord[i]),
                 &(pchain_coord[i+1]), &(pchain_coord[i+2]));
          i+=3;
        }
    }

  /* Now for the inhibitor... */
  /* ------------------------ */

  /* Prompt for inhibitor chain */
  printf("\n");
  for (i=0; i < cplx_chains; i++)
    if (chain_id[i] != BLANK)
      printf("Chain: %d\tIdentifier: %c\n", i+1, chain_id[i]);
  printf("\n");

 /*If interactive, ask for chain ID */
  if(argc==1) {
    /* Prompt for inhibitor chain */
    ibuf[0] = '\0';
    while (!memchr(chain_id, ibuf[0], cplx_chains))
      {
      if (ibuf[0] != '\0') printf("\nNot a valid chain identifier!\n");
      printf("Please enter the chain identifier for the inhibitor...\n");
      fgets(ibuf, MAXLN - 1, stdin);
    }
  }
  else { /* Non interactive, chain ID argv[4] */
    strcpy(ibuf,argv[4]);
    printf("Using %s as chain ID for the inhibitor chain\n",ibuf);
  }
  
  
  /* Find which chain number matches the chain id given */
  ichain_num = 0;
  while ((ichain_num < cplx_chains) && (chain_id[ichain_num] != ibuf[0]))
    ichain_num++;

  /* Load the coordinates for all atoms in the inhibitor */
  /* =================================================== */

  /* Move to the first ATOM record */
  fseek(cplxpdb, first_atom[ichain_num], SEEK_SET);

  /* Allocate space for coordinate array */
  ichain_coord = (float *) malloc(sizeof(float)*chain_len[ichain_num]*3);

  /* Load the coordinates */
  i=0;
  while (!feof(cplxpdb))
    {
      /* If we find an ATOM or HETATM with  the correct chain ID... */
      /* then load it...                                            */
      if ((fgets(fbuf, MAXLN, cplxpdb)) &&
          (fbuf[ID_COL] == chain_id[ichain_num]) &&
          ((!memcmp(fbuf, "ATOM", 4)) || (!memcmp(fbuf, "HETATM", 6))))
        {
          sscanf(fbuf, "%*31c %f %f %f", &(ichain_coord[i]),
                 &(ichain_coord[i+1]), &(ichain_coord[i+2]));
          i+=3;
        }
    }

  /* =============================================================== */
  /* Step 2:  Check each water in the complex PDB file for proximity */
  /* to any atom in the inhibitor, dump those which are closer than  */
  /* HBOND_CUTOFF angstroms to both the inhibitor and the protein.   */
  /* Then dump any waters which are closer than HBOND_CUTOFF ang-    */
  /* stroms to any water already dumped.                             */
  /* =============================================================== */

  /* We already know how many waters are in the complex file, but we */
  /* still need to scan the free file...                             */
  fw_found = 0;
  rewind(freepdb);
  while (!feof(freepdb))
    {
      /* Find the next HETATM */
      fbuf[0] = '\0';
      while((memcmp(fbuf, "HETATM", 4)) && !feof(freepdb))
        fgets(fbuf, MAXLN, freepdb);

      /* check the type */
      if (!feof(freepdb))
        {
          sscanf(fbuf, "%*16c %s", het_type);
          /* If it is a water, increment the count */
          /*if (!memcmp(het_type, "HOH", 3))*/
	  if(iswater(het_type))
            fw_found++;
        }
    }

  /* Allocate space for the water arrays */
  cplx_waters = (water_rec *)
                malloc(sizeof(water_rec)*cw_found);
  free_waters = (water_rec *)
                malloc(sizeof(water_rec)*fw_found);

  strcpy(p.struct_type,"Complex\0");
  p.waters=cplx_waters;
  p.pdb=cplxpdb;
  if(argc==1)
    p.outfile=NULL;
  else
    {
      char tmp[MAXLN+1];
      strcpy(tmp,argv[5]);
      strcat(tmp,".active.hits");
      p.outfile=get_file("","w",tmp,FALSE);
    }
  p.w_found=cw_found;
  p.pchain_len=chain_len[pchain_num];
  p.ichain_len=chain_len[ichain_num];
  p.pchain_coord=pchain_coord;
  p.ichain_coord=ichain_coord;

  process(p);
 
  strcpy(p.struct_type,"Free\0");
  p.waters=free_waters;
  p.pdb=freepdb;
  if(argc==1)
    p.outfile=NULL;
  else
    {
      char tmp[MAXLN+1];
      strcpy(tmp,argv[6]);
      strcat(tmp,".active.hits");
      p.outfile=get_file("","w",tmp,FALSE);
    }
  p.w_found=fw_found;
  p.pchain_len=chain_len[pchain_num];
  p.ichain_len=chain_len[ichain_num];
  p.pchain_coord=pchain_coord;
  p.ichain_coord=ichain_coord;

  process(p);

  printf("Done.\n");

  /* ========================== */
  /* Step 4:  Clean up and exit */
  /* ========================== */

  free (ichain_coord);
  free (pchain_coord);
  free (cplx_waters);
  free (free_waters);
  fclose(cplxpdb);
  fclose(freepdb);

  return (0);
}


