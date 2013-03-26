// Ancestry composition preprocessor

#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <map>
#include "ancescomp.h"
#include "ethnicity.h"

FILE *ancestryfile,*genomefile;
map<int64_t,interval> haploid[2],diploid;
/* The longest chromosome is the first, which is 250 megabases (ee00000) long.
 * The index is thus ccssssssseeeeeee, where c is the chromosome number (00-18,
 * where 00=mt, 17=X, 18=Y, and 00 and 18 aren't used in the map), s is the start,
 * and e is the end with its bits flipped.
 */

void interval::clear()
{
  chromosome=start=end=0;
  ethnicity[0]=ethnicity[1]=0;
}

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

int chromosomenumber(string str)
{
  str.erase(0,3);
  if (str[0]=='X')
    return 23;
  else
    return atoi(str.c_str());
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
  string ethnicity;
  interval intvl;
  int hap,chr,startend=0;
  int64_t index;
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
	  startend=0;
	  break;
	case ']':
	  bracketindent--;
	  if (startend)
	  {
	    intvl.ethnicity[0]=find_ethnic(ethnicity);
	    index=(((int64_t)chr)<<56)+(((int64_t)intvl.start)<<28)+(intvl.end^0xfffffff);
	    if (haploid[hap][index].ethnicity[0]<=intvl.ethnicity[0])
	      haploid[hap][index]=intvl;
	    printf("%d %016lx [%d,%d] %s\n",hap,index,intvl.start,intvl.end,ethnicity.c_str());
	    intvl.clear();
	    startend=0;
	  }
	  break;
	case ',':
	  startend++;
	  break;
	case '"':
	  switch (braceindent)
	  {
	    case 2:
	      ethnicity=tok.str;
	      break;
	    case 3:
	      hap=tok.str[3]-'1';
	      break;
	    case 4:
	      chr=chromosomenumber(tok.str);
	      break;
	  }
	  /*for (i=0;i<braceindent;i++)
	    printf("{ ");
	  for (i=0;i<bracketindent;i++)
	    printf("[ ");
          printf("%s\n",tok.str.c_str());*/
	  break;
	case '0':
	  if (startend)
	    intvl.end=tok.n;
	  else
	    intvl.start=tok.n;
	  /*for (i=0;i<braceindent;i++)
	    printf("{ ");
	  for (i=0;i<bracketindent;i++)
	    printf("[ ");
          printf("%d\n",tok.n);*/
	  break;
      }
    }
    fclose(ancestryfile);
  }
  else
    fprintf(stderr,"%s: could not open file\n",ancestryname);
}

void sortancestry()
{
  map<int64_t,interval>::iterator i,j;
  int hap;
  for (hap=0;hap<2;hap++)
    for (i=haploid[hap].begin();i!=haploid[hap].end();i++)
      printf("%d %016lx [%d,%d] %s\n",hap,i->first,i->second.start,i->second.end,ethnicities[i->second.ethnicity[0]].c_str());
}
      
void usage()
{
  printf("Usage: ancescomp ancestrydata genome\nTo get your ancestry data, bring up\n");
  printf("https://www.23andme.com/you/ancestry/composition/fetch/?profile_id=<id>&threshold=0.75\n");
  printf("in a web browser and save it to a file.\n");
}

int main(int argc,char **argv)
{
  init_ethnic();
  if (argc<2 || argc>3)
    usage();
  else
  {
    readancestry(argv[1]);
    sortancestry();
  }
  return 0;
}
