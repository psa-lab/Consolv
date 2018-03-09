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
/* util.c - utility functions                                          */
/***********************************************************************/
#include "defines.h"
#include <ctype.h>

/* Upcase simply converts a null-terminated string to upper case */
void upcase (char *buffer) {
  char *current;		/* current character */

  current = buffer;
  while (*current != '\0') {
    if (islower(*current))
      *current = toupper(*current);
    current++;
  }
}

/* Isblank just checks if there are any non-whitespace characters */
/* in a string */
int isblank (char *buffer) {
  char *current;

  current = buffer;
  while (*current != '\0')
    if (*current != ' '
	&& *current != '\t'
	&& *current != '\n')
      return FALSE;
    else
      current++;
  
  return TRUE;
}
