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
/* conswat - Find conserved and displaced waters in a pair of PDB files.   */
/*                                                                         */
/* This program will prompt the user for two pdb files.  One should        */
/* contain a protein in "free" form - and another should contain a         */
/* protien with a bound substrate.  The program will compare all the       */
/* in the free form with those in the bound form, and decide whether each  */
/* one is conserved (present in both files) or displaced (present only in  */
/* the free form file.  A water is considered to be conserved if there is  */
/* a corresponding water in the bound form within CNSRV_CUTOFF angstroms   */
/* of the position of the free water.  If more than one water is found,    */
/* the closer of the two is used.  Each water in the bound form can only   */
/* correspond to a single water in the free form.                          */
/***************************************************************************/
/* Programmers: Michael Raymer, Vishal Thakkar                             */
/* Latest Revision: Date June 6, 1997.                                     */
/***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>

/* Here are some global constants, if this grows into a mult-module */
/* application, then these will be moved to a global.h file         */
#define CNSRV_CUTOFF 1.2	   /* "Close" water distance cuttoff */
#define NODIST -1		   /* Initial value for distance */
#define MAXLN 132		   /* Max input line length */
#define MAXPDB 6		   /* Max length of a PDB code */
/*#define PDBLOC "/psa/pdb"*/	   /* Location of PDB files */
#define BLANK ' '		   /* Blank character */
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

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
/*          interactive = if 1, then will prompt for name and overwriting  */
/*                        if 0, will not prompt, and will rename before    */
/*                        overwriting the file                             */
/*                                                                         */
/* OUTPUTS:  get_file returns a pointer to the open file.                  */
/***************************************************************************/
FILE *get_file(prmpt, omode, fdefault,interactive)
     char *prmpt;		/* String to prompt user with */
     char *omode;		/* Mode to open the file in */
     char *fdefault;		/* The default file name, if any */
     int interactive;		/* interactive=1, move the file if overwriting */
{
  /* Local Variables */
  char input_fn[MAXLN + 1];	   /* User inputted filename */
  FILE *fptr;			   /* Pointer to opened file */

  /* Initializations */
  input_fn[0] = '\0';
  fptr = (FILE *) NULL;

    while (!fptr)
      {
	if(interactive==FALSE) {
	  strcpy(input_fn,fdefault);
	}
	else /* Get a valid filename */
	  while (strlen(input_fn) <= 1 && interactive==TRUE)
	    {
	      /* prompt user */
	      printf("\n%s\n", prmpt);
	      if (strlen(fdefault) >= 1)
		printf("(press return to use [%s])\n\n--> ", fdefault);
	      
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
	    
	    if(stat(input_fn,&b)==NULL) {/* File exit */
	      if(interactive==TRUE)
		{/* Ask user */
		  printf("file: %s exits, Overwrite [y/n] ?",input_fn);
		  scanf("%c",&ans);
		  fflush(stdin);  /*Flush the input buffer of any carriage returns */
		  ans=toupper(ans);
		  if(ans=='N')
		    return(NULL); /* Return Null pointer */
		}
	      else{ /* don't ask user, just rename the file*/
		strcpy(bakname,input_fn);
		strcat(bakname,".bak");
		printf("\nwarning:moving file %s to %s.\n",input_fn,bakname);
		rename(input_fn,bakname);
	      }
	    }
	  }/* end of checking...*/
	
	/* Attempt to open the file */
	fptr = fopen(input_fn, omode);
	
	/* If we fail, reset the filename and print an error */
      if (!fptr)
	{
	  printf("\nUnable to open file: %s\n\n", input_fn);
	  input_fn[0] = '\0';
	  if(interactive==FALSE) return 0;
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
  enum stat_type { matched, unmatched } status;
} water_rec;

main(int argc, char *argv[])
{

  /* Local Variables */
  int i;			   /* Loop counters, etc. */
  int found;			   /* A "done" flag */
  char prmpt[MAXLN + 1];	   /* Buffer for user prompts */
  char fdef[MAXLN + 1];		   /* Buffer for default filenames */
  char fbuf[MAXLN + 1];		   /* Buffer for input file lines */
  char ibuf[MAXLN + 1];		   /* Buffer for user input */
  char pdbcode[MAXPDB+1];	   /* PDB code for the free protein */
  FILE *cplxpdb;		   /* Pointer to complex PDB file */
  FILE *freepdb;		   /* Pointer to "free" PDB file  */
  FILE *outfile;		   /* Pointer to output file */
  char cur_type[MAXLN];		   /* Type of the current HETATM */
  float cur_x, cur_y, cur_z;	   /* Coords of current HETATM */
  float dist;			   /* Currrent dist from water to inhibitor */
  water_rec *closest_wat;	   /* Points to the current closest water */
  float cur_dist;		   /* Current closest water distance */
  int cur_num;			   /* Number of current water */
  int wats_found;		   /* Number of waters in complex form */
  int dsplcd_found;                /* Number of displaced waters found */
  int cnsrvd_found;		   /* Number of conserved waters found */
  water_rec *cplx_waters;          /* Waters in the complex form */


  switch(argc){

  case 1:{/* Interactive execution */
    /* Print splash banner */
    printf("conswats - find conserved and displaced waters in a protein\n");
    printf("===========================================================\n\n");
    
    /* Open the two PDB files... */
    /* ========================= */
    
    /* Get the Free file PDB code */
    printf("\nWhat is the PDB code for the free form protein? --> ");
    fgets(pdbcode, MAXPDB, stdin);
    
    /* Get rid of the linefeed on the end */
    pdbcode[strlen(pdbcode)-1]='\0';
    
    /* Free file */
    sprintf(fdef, "%spdb%s.xfm","", pdbcode); /* Replace 3rd parameter with required path */
    strcpy(prmpt,"Enter filename for the Free-structure file:");
    freepdb = get_file(prmpt, "r", fdef,TRUE);
    
    
    /* Get the complex file pdb code */
    printf("\nWhat is the PDB code for the complex form? --> ");
    fgets(ibuf, MAXPDB, stdin);
    
    /* Get rid of the linefeed on the end */
    ibuf[strlen(ibuf)-1]='\0';
    
    /* Complex file */
    sprintf(fdef, "%s/pdb%s.ent", PDBLOC, ibuf);
    strcpy(prmpt, "Enter filename for the protein:inhibitor coordinates...");
    cplxpdb = get_file(prmpt, "r", fdef,TRUE);
    
    /* Opening Output file */
    /* filename fixed to PDBCODE.cons */
    sprintf(fdef, "%s.cons", pdbcode);
    outfile = get_file("","w",fdef,TRUE);
    if(outfile==NULL) {
      printf("Error opening %s. Exiting...\n",fdef);          
      exit(-1);
    }
   
  }
  break;
  
  case 4:
    {
      freepdb=fopen(argv[1],"r");
      cplxpdb=fopen(argv[2],"r");
      if(freepdb==NULL) {
	printf("Error opening %s.Exiting...\n",argv[1]);
	exit(-1);
      }
      if(cplxpdb==NULL) {
	printf("Error opening %s.Exiting...\n",argv[2]);
	exit(-1);
      }
      strcpy(pdbcode,argv[3]);

      /* Opening Output file */
      /* filename fixed to PDBCODE.cons */
      sprintf(fdef, "%s.cons", pdbcode);
      outfile = get_file("","w",fdef,FALSE);
      if(outfile==NULL) {
	printf("Error opening %s. Exiting...\n",fdef);
	exit(-1);
      }
	
  }
  break;

  default:
    printf("Conswat: finds conserved/displaced waters in free and complex forms\n");
    printf("Usage: %s freefile cplxfile freePDBCode\n",argv[0]);
    exit(-1);
    break;
  }/* finished processing the commandline params,if any. */
  


  /* ========================================================== */
  /* Step 1: Load the coordinates of each atom in the inhibitor */
  /*         and in the protein.                                */
  /* ========================================================== */

  /* Pass 1 - Scan the file, determine the sizes for the       */
  /*          various arrays, assemble the lists of chain      */
  /*          id's chain lengths, starting points, etc.        */
  /* --------------------------------------------------------- */
  

  printf("\nExamining complex PDB file...\n\n");
  fbuf[0] = '\0';
  wats_found = 0;
  while (!feof(cplxpdb))
    {
      /* Find the next HETATM record */
      while((memcmp(fbuf, "HETATM", 6)) && (!feof(cplxpdb)))
	fgets(fbuf, MAXLN, cplxpdb);

      /* If we found a HETATM, check if it is a water */
      if (!memcmp(fbuf, "HETATM", 6))
	{
	  /* If this is a water, increment the water count */
	  sscanf(fbuf, "%*16c %s", cur_type);
	  if ((!memcmp(cur_type, "HOH", 3)) ||
	      (!memcmp(cur_type, "WAT", 3)) ||
	      (!memcmp(cur_type, "H2O", 3)) ||
	      (!memcmp(cur_type, "DOD", 3)) ||
	      (!memcmp(cur_type, "D2O", 3)))
	    wats_found++;
	}

      /* Get the next line */
      fgets(fbuf, MAXLN, cplxpdb);
    }

  /* Pass 2 - Load the water array */
  /* ----------------------------------- */

  rewind(cplxpdb);

  /* Allocate space for water array */
  cplx_waters = (water_rec *) malloc(sizeof(water_rec)*wats_found);

  /* Load the coordinates */
  i=0;
  while (!feof(cplxpdb))
    {
      if ((fgets(fbuf, MAXLN, cplxpdb)) &&
	  (!memcmp(fbuf, "HETATM", 6)))       
	{
	  sscanf(fbuf, "%*16c %s", cur_type);
	  if ((!memcmp(cur_type, "HOH", 3)) ||
	      (!memcmp(cur_type, "WAT", 3)) ||
	      (!memcmp(cur_type, "H2O", 3)) ||
	      (!memcmp(cur_type, "DOD", 3)) ||
	      (!memcmp(cur_type, "D2O", 3)))
	    {

	      sscanf(fbuf, "%*31c %f %f %f", &(cplx_waters[i].x_coord),
		     &(cplx_waters[i].y_coord), &(cplx_waters[i].z_coord));
	      cplx_waters[i].status = unmatched;
	      i++;
	    }
	}
    }

  /* =============================================================== */
  /* Step 2:  Check each water in the free PDB file for proximity    */
  /* to any atom in the complex form.  Dump each water ID, and a tag */
  /* identifying it as conserved or dislpaced.  A water is tagged as */
  /* conserved if it is found in BOTH the complex and the free form. */
  /* Since the water may have moved a little in the ligation process */
  /* a tolerance of CNSRV_CUTOFF angstroms is used to identify       */
  /* conserved waters.                                               */
  /* =============================================================== */

  /* Scan the free file for waters */
  printf("\nScanning free file for waters...\n");
  cnsrvd_found = 0;
  dsplcd_found = 0;
  while (!feof(freepdb))
    {
      /* Find the next HETATM */
      fbuf[0] = '\0';
      while((memcmp(fbuf, "HETATM", 6)) && !feof(freepdb))
	fgets(fbuf, MAXLN, freepdb);

      /* Get the coordinates and the type */
      if (!feof(freepdb))
	{
	  sscanf(fbuf, "%*16c %s %i %f %f %f",
		 cur_type, &cur_num, &cur_x, &cur_y, &cur_z);

	  /* Check if it is a water... */
	  if ((!memcmp(cur_type, "HOH", 3)) ||
	      (!memcmp(cur_type, "WAT", 3)) ||
	      (!memcmp(cur_type, "H2O", 3)) ||
	      (!memcmp(cur_type, "DOD", 3)) ||
	      (!memcmp(cur_type, "D2O", 3)))
	    {
	      /* If so, compare its position to unmatched waters in the */
	      /* complex form...                                        */
	      /* ------------------------------------------------------ */

	      /* Find the closest unmatched water... */
	      cur_dist = NODIST;
	      for (i=0; i<wats_found; i++)
		{
		  if (cplx_waters[i].status == unmatched)
		    {
		      dist = hypot(hypot(cplx_waters[i].x_coord - cur_x,
					 cplx_waters[i].y_coord - cur_y),
				   cplx_waters[i].z_coord - cur_z);
		      if ((cur_dist == NODIST) || (dist < cur_dist))
			{
			  cur_dist = dist;
			  closest_wat = &(cplx_waters[i]);
			}
		    }
		}

	      /* If the closest unmatched water is within CNSRV_CUTOFF */
	      /* angstroms to this water - match and dump...           */
	      if ((cur_dist != NODIST) && (cur_dist <= CNSRV_CUTOFF))
		{
		  closest_wat->status = matched;
		  fprintf(outfile, "%s\t%i\t1\n", pdbcode, cur_num);
		  cnsrvd_found++;
		}
	      else
		{
		  fprintf(outfile, "%s\t%i\t0\n", pdbcode, cur_num);
		  dsplcd_found++;
		}
	    }
	}
    }
  
  printf("Done extracting, waters found --\n");
  printf("\tConserved:\t%4i\n", cnsrvd_found);
  printf("\tDisplaced:\t%4i\n", dsplcd_found);
  printf("\tTOTAL:\t\t%4i\n\n", cnsrvd_found + dsplcd_found);


  /* ========================== */
  /* Step 4:  Clean up and exit */
  /* ========================== */

  free (cplx_waters);
  fclose(cplxpdb);
  fclose(freepdb);
  fclose(outfile);

}
