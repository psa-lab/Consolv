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

/**************************************************************************************/
/* newscale.c                                					      */
/*         * This program reads the .env file and normalizes the data for each feature*/
/*           Predetermined max and min values are used as reference. Min value corres */
/*           ding 1 and max value corresponds to 10.                                  */
/*         * The  program also reads data from the .cons file which contains the list */
/*           of conserved/unconserved waters and adds it to the output file           */
/*         * For making it compatible with bbknn runs, if water is                    */
/*               conserved = 2                                                        */
/*             unconserved = 1                                                        */
/*         * If the water is an activesite water then                                 */
/*               activesite = 2                                                       */
/*            nonactivesite = 1                                                       */
/*										      */
/* Programmer: Vishal Thakkar. Date: Dec 18, 1996.				      */
/**************************************************************************************/


/* Revision
 * April 18, 1997.
 * - separated hbd to hbd with protein atoms and hbd with waters
 * - added netbval, avgbval and active sites
 *
 */

/* Revision
 * June 6, 1997.
 *  - Added features so that it can run with commandline, noninteractivly too.
 */

/*
 * Revision
 * Dec 17, 1997.
 * Changed the order of NetAHP and AvgAHP comments in the output scaled file 
 */

#include <stdio.h>
#include <string.h>
#define FNAME_LEN 120
/* New min/max values extracted from the new 1st shell waters database */
/* A cutoff of 99.5% was chosen. */

#define AHP_MAX 3.9
#define AHP_MIN 0.078
#define ADN_MAX 10.0
#define ADN_MIN  1.0
#define HBD_P_MAX 4.0
#define HBD_P_MIN  0.0
#define HBD_W_MAX 7.0 
#define HBD_W_MIN 0.0 
#define BVAL_MAX 80.0
#define BVAL_MIN  2.0
#define MOB_MAX 2.28
#define MOB_MIN 0.054
#define AVGBVAL_MAX 68.0
#define AVGBVAL_MIN 6.25
#define NETBVAL_MAX 320.0
#define NETBVAL_MIN 3.0

#define APP_MODE 1
#define TEST_MODE 0

void process_array(float *,int,float, float);


void process_array(float *table, int size,float max,float min)
{
 int i;
 float tmp;


/* scaling the array which is passed  */
/* formula:
    new = (old-min)/(max-min)*9 + 1;      */
 tmp=max-min;
for(i=0;i<size;i++)
   table[i]= (table[i]-min)/tmp*9 + 1.0;

/* returning */
}

int main(int argc,char *argv[])
{
 int num_of_lines,i,j;
 float *ahptable, *adntable, *hbd_P_table, *hbd_W_table;
 float *bvaltable, *mobtable, *avgbvaltable, *netbvaltable;
 int *constable, *activetable;
 FILE *fp, *fout;
 char str[FNAME_LEN+1];
 int mode;

 /* Determine the mode of operation */
 /* if test mode, ask for .cons file and then append it to the end */
 /* if application mode, don't ask for the .cons file */
 /* Also determine whether to go interactive or noninteractive */

 switch(argc) {
 case 2:{
   /* Interactively*/
   if(argv[1][1]!='a' && argv[1][1]!='t') {
     printf("Invalid mode. Use -a or -t.\n");
     exit(-1);
   }
   /*find execution mode */
   if(argv[1][1]=='t') 
     mode=TEST_MODE;
   else
     mode=APP_MODE;

   /* The program does this */
   printf("\nscale: Produces a output file for bbknn runs.");
   printf("\n=============================================");
   
   /* Read a string from keyborad and strip the \n from it */
   printf("\nGive the environ filename: ");
   fgets(str,FNAME_LEN,stdin);
   i=strlen(str);
   str[i-1]=NULL;
   
   /* Open the file */
   fp=fopen(str,"r");
   if(fp==NULL){ 
     printf("Error Opening file %s\nexiting...",str);
     exit(-1);
   }
   /* Now the output file */
   fflush(stdin);
   printf("Output filename: ");
   fgets(str,FNAME_LEN,stdin);
   j=strlen(str);str[j-1]=NULL;
   fout=fopen(str,"w");
   if(fout==NULL){ 
     printf("Error Opening file %s\nexiting...",str);
     exit(-1);
   }
   /* Finished reading the input/outputfiles*/
 }
 break;
  
 case 4:{
   /*Noninteractively ...*/
   if(argv[1][1]!='a' && argv[1][1]!='t') {
     printf("Invalid mode. Use -a or -t.\n");
     exit(-1);
   }
   /*find execution mode */
   if(argv[1][1]=='t') 
     mode=TEST_MODE;
   else
     mode=APP_MODE;

   fp=fopen(argv[2],"r");
   fout=fopen(argv[3],"w");
   if(fp==NULL||fout==NULL) {
	if(fp==NULL)
	   printf("Error opening %s.",argv[2]);
	else
	   printf("Error opening %s.",argv[3]);
	printf("\n.Exiting...\n");
	exit(-1);
   }
				 
 }
 break;
 
 default:
   { 
     printf("Usage: %s [-a|-t] infile scaled_outfile\n",argv[0]);
     printf("-a    Application mode.\n");
     printf("-t    Test mode.\n");
     exit(-1);
   }
 
 }  

 /* Pass 1. Find the number of lines in the file */
 /* We are ignoring the lines which start with # */
 /* We are using the following construct as feof will be nonzero when fgets is used AFTER */
 /* the last SUCCESSFUL fgets, ie, after the file is read, one more file read is required */
 /* to trigger fgets */
 num_of_lines=0;

 fgets(str,FNAME_LEN,fp);
 while((feof(fp)==NULL))
  {
     i=strlen(str);str[i-1]=NULL;
     if(str[0]!='#')
       num_of_lines++;
     fgets(str,FNAME_LEN,fp);
  } /* end of while */

 /* End of Pass 1. */

 rewind(fp);
 /* allocating memory for each feature array */
 ahptable     =  (float *)malloc(num_of_lines*sizeof(float));
 adntable     =  (float *)malloc(num_of_lines*sizeof(float));
 bvaltable    =  (float *)malloc(num_of_lines*sizeof(float));
 hbd_P_table  =  (float *)malloc(num_of_lines*sizeof(float));
 hbd_W_table  =  (float *)malloc(num_of_lines*sizeof(float));
 mobtable     =  (float *)malloc(num_of_lines*sizeof(float));
 avgbvaltable =  (float *)malloc(num_of_lines*sizeof(float));
 netbvaltable =  (float *)malloc(num_of_lines*sizeof(float));
 
 if(mode==TEST_MODE) {
  constable=(int *)malloc(num_of_lines*sizeof(int));
  activetable=(int *)malloc(num_of_lines*sizeof(int));
 }

 /*Pass 2. Now read each of the data into respective array */
 i=0;
 fgets(str,FNAME_LEN,fp);
 while(!feof(fp))
  {
    j=strlen(str);str[j-1]=NULL;
    if(str[0]!='#')
      {/* Here we are reading the adn,ahp,bval,hbd,mob from */
	/* the string neglecting the sum_of_bvals field */
	  if(mode==TEST_MODE)
	  sscanf(str,"%*s %*s %f %f %f %f %f %f %f %f %d %d\n",
		 &adntable[i],
		 &ahptable[i],
		 &bvaltable[i],
		 &hbd_P_table[i],
		 &hbd_W_table[i],
		 &mobtable[i],
		 &netbvaltable[i],
		 &avgbvaltable[i],
		 &constable[i],
		 &activetable[i]);
	else
	  sscanf(str,"%*s %*s %f %f %f %f %f %f %f %f %*d %*d\n",
		 &adntable[i],
		 &ahptable[i],
		 &bvaltable[i],
		 &hbd_P_table[i],
		 &hbd_W_table[i],
		 &mobtable[i],
		 &netbvaltable[i],
		 &avgbvaltable[i]);
	i++;
         }/* end of if */
    fgets(str,FNAME_LEN,fp);
  } /* end of while */
 /* End of Pass 2 */

 fclose(fp);

 /* Now processing the feature tables */

 process_array(adntable,num_of_lines,ADN_MAX,ADN_MIN);    
 process_array(ahptable,num_of_lines,AHP_MAX,AHP_MIN);
 process_array(bvaltable,num_of_lines,BVAL_MAX,BVAL_MIN);
 process_array(hbd_P_table,num_of_lines,HBD_P_MAX,HBD_P_MIN);
 process_array(hbd_W_table,num_of_lines,HBD_W_MAX,HBD_W_MIN);
 process_array(mobtable,num_of_lines,MOB_MAX,MOB_MIN);
 process_array(avgbvaltable,num_of_lines,AVGBVAL_MAX,AVGBVAL_MIN);
 process_array(netbvaltable,num_of_lines,NETBVAL_MAX,NETBVAL_MIN);

 fprintf(fout,"#Scaled feature file for bbknn.\n");
 fprintf(fout,"#         Min    Max values used for normalizing:\n");
 fprintf(fout,"#Adn:      %.3f      %.3f\n",ADN_MIN,ADN_MAX);
 fprintf(fout,"#Ahp:      %.3f      %.3f\n",AHP_MIN,AHP_MAX);
 fprintf(fout,"#Bval:     %.3f      %.3f\n",BVAL_MIN,BVAL_MAX);
 fprintf(fout,"#hbd_P:    %.3f      %.3f\n",HBD_P_MIN,HBD_P_MAX);
 fprintf(fout,"#hbd_W:    %.3f      %.3f\n",HBD_W_MIN,HBD_W_MAX);
 fprintf(fout,"#Mobility: %.3f      %.3f\n",MOB_MIN,MOB_MAX); 
 fprintf(fout,"#Net Bval: %.3f      %.3f\n",NETBVAL_MIN,NETBVAL_MAX);
 fprintf(fout,"#Avg Bval: %.3f      %.3f\n",AVGBVAL_MIN,AVGBVAL_MAX);
 fprintf(fout,"#\n");
 fprintf(fout,"Features: 8\n");
 if(mode==TEST_MODE)
   for(i=0;i<num_of_lines;i++)
   {
     fprintf(fout,"%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %d %d\n",
           adntable[i],
	   ahptable[i],
	   bvaltable[i],
	   hbd_P_table[i],
	   hbd_W_table[i],
	   mobtable[i],
	   netbvaltable[i],
	   avgbvaltable[i],
	   constable[i]+1,
	   activetable[i]+1);
   }
 else
   for(i=0;i<num_of_lines;i++)
   { 
     fprintf(fout,"%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
	   adntable[i],
	   ahptable[i],
	   bvaltable[i],
	   hbd_P_table[i],
	   hbd_W_table[i],
	   mobtable[i], 
	   netbvaltable[i],
	   avgbvaltable[i] );
   }

 /* Free the memory buffers and close the open files  */
 free(ahptable);
 free(adntable);
 free(bvaltable);
 free(hbd_P_table);
 free(hbd_W_table);
 free(mobtable);
 if(mode==TEST_MODE) { 
   free(constable);
   free(activetable);
 }
 fclose(fout);

 puts("Done!\n");
 exit(0);
}








