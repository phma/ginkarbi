// Ancestry composition preprocessor

#include <cstdio>
#include "ancescomp.h"

FILE *ancestryfile,*genomefile;

token readtoken()
// Tokens are space, comma, left and right brace, left and right bracket, EOF, string, and number.
// All whitespace characters are returned as space.
// There are unquoted strings, which are returned as single letter tokens. I ignore them.
{
  token tok;
  bool done=false,number=false,quote=false;
  tok.n=0;
  tok.str="";
  while (!done)
  {
    tok.ch=fgetc(ancestryfile);
    if (quote)
      if (tok.ch=='"')
      {
	done=true;
	quote=false;
      }
      else if (tok.ch==EOF)
	done=true;
      else
        tok.str+=tok.ch;
    else // !quote
    {
      if (isspace(tok.ch))
        tok.ch=' ';
      if (isdigit(tok.ch))
      {
        number=true;
        tok.n=tok.n*10+tok.ch-'0';
        tok.ch='0';
      }
      else if (number)
      {
	ungetc(tok.ch,ancestryfile);
	tok.ch='0';
	done=true;
      }
      else if (tok.ch=='"')
	quote=true;
      else
	done=true;
    }
  }
  return tok;
}

void readancestry(char *ancestryname)
/* These occur at the following braceindent levels:
 * 1: trio, genotyped, owned, split, active, segments. Ignored.
 * 2: north_european, eana, and other ethnic categories. They occur in no particular order and must be sorted later.
 * 3: hap1 and hap2, which designate which chromosome of a pair.
 * 4. chr1 through chr22 and chrX-npar. chrY and mt do not occur.
 */
{
  token tok;
  int braceindent=0,bracketindent=0,i;
  ancestryfile=fopen(ancestryname,"r");
  if (ancestryfile)
  {
    while (tok.ch!=EOF)
    {
      tok=readtoken();
      switch (tok.ch)
      {
	case '{':
	  braceindent++;
	  break;
	case '}':
	  braceindent--;
	  break;
	case '[':
	  bracketindent++;
	  break;
	case ']':
	  bracketindent--;
	  break;
	case '"':
	  for (i=0;i<braceindent;i++)
	    printf("{ ");
	  for (i=0;i<bracketindent;i++)
	    printf("[ ");
          printf("%s\n",tok.str.c_str());
	  break;
	case '0':
	  for (i=0;i<braceindent;i++)
	    printf("{ ");
	  for (i=0;i<bracketindent;i++)
	    printf("[ ");
          printf("%d\n",tok.n);
	  break;
      }
    }
    fclose(ancestryfile);
  }
  else
    fprintf(stderr,"%s: could not open file\n",ancestryname);
}

void usage()
{
  printf("Usage: ancescomp ancestrydata genome\nTo get your ancestry data, bring up\n");
  printf("https://www.23andme.com/you/ancestry/composition/fetch/?profile_id=<id>&threshold=0.75\n");
  printf("in a web browser and save it to a file.\n");
}

int main(int argc,char **argv)
{
  if (argc<2 || argc>3)
    usage();
  else
  {
    readancestry(argv[1]);
  }
  return 0;
}
