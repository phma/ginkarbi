// Ancestry composition preprocessor

#include <cstdio>
#include <cstdlib>
#include <map>
#include "ancescomp.h"
#include "ethnicity.h"
#include "snpfile.h"

#define max(a,b) ((a>b)?(a):(b))
#define min(a,b) ((a<b)?(a):(b))

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

int64_t interval::index()
{
  return (((int64_t)chromosome)<<56)+(((int64_t)start)<<28)+(end^0xfffffff);
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
	    intvl.chromosome=chr;
	    index=intvl.index();
	    if (haploid[hap][index].ethnicity[0]<=intvl.ethnicity[0])
	      haploid[hap][index]=intvl;
	    //printf("%d %016lx [%d,%d] %s\n",hap,index,intvl.start,intvl.end,ethnicity.c_str());
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
  map<int64_t,interval> temploid;
  interval intvl0,intvl1;
  bool cont,split;
  int hap,chr,pos,end0,end1;
  do
  {
    for (cont=false,hap=0;hap<2;hap++)
    {
      /*         (American)
       * (East Asian and North American)
       * becomes
       * ( EANA )(American)(   EANA    )
       */
      temploid.clear();
      for (i=haploid[hap].begin();i!=haploid[hap].end();i++)
      {
	j=i;
	j++;
	split=false;
	if (j!=haploid[hap].end() && j->second.end<=i->second.end && j->second.chromosome==i->second.chromosome && j->second.ethnicity[0]>i->second.ethnicity[0])
	{
	  split=cont=true;
	  intvl0=intvl1=i->second;
	  intvl0.end=j->second.start;
	  intvl1.start=j->second.end;
	  if (intvl0.start<intvl0.end)
	    temploid[intvl0.index()]=intvl0;
	  if (intvl1.start<intvl1.end)
	    temploid[intvl1.index()]=intvl1;
	}
	else
	  temploid[i->second.index()]=i->second;
	//printf("%d %016lx [%d,%d] %s %s\n",hap,i->first,i->second.start,i->second.end,ethnicities[i->second.ethnicity[0]].c_str(),split?"(split)":"");
      }
      haploid[hap]=temploid;
    }
  } while (cont);
  diploid.clear();
  /* hap1: (    French     )(            African                  )
   * hap2: (            American             )(   South Asian     )
   * becomes
   * dip : (French/American)(African/American)(African/South Asian)
   */
  for (chr=pos=0,i=haploid[0].begin(),j=haploid[1].begin();i!=haploid[0].end() || j!=haploid[1].end();)
  {
    if (i==haploid[0].end())
    {
      intvl0=j->second;
      intvl0.ethnicity[1]=intvl0.ethnicity[0];
      intvl0.ethnicity[0]=-1;
      j++;
    }
    else if (j==haploid[1].end())
    {
      intvl0=i->second;
      intvl0.ethnicity[1]=-1;
      i++;
    }
    else
    {
      if (chr<min(i->second.chromosome,j->second.chromosome))
      {
	chr=min(i->second.chromosome,j->second.chromosome);
	pos=0;
      }
      if (pos<min(i->second.start,j->second.start))
	pos=min(i->second.start,j->second.start);
      intvl0.chromosome=chr;
      intvl0.start=pos;
      if (chr<i->second.chromosome || (chr==i->second.chromosome && pos<i->second.start))
      {
	intvl0.ethnicity[0]=-1;
	end0=(chr==i->second.chromosome)?i->second.start:0x7fffffff;
      }
      else
      {
	intvl0.ethnicity[0]=i->second.ethnicity[0];
	end0=i->second.end;
      }
      if (chr<j->second.chromosome || (chr==j->second.chromosome && pos<j->second.start))
      {
	intvl0.ethnicity[1]=-1;
	end1=(chr==j->second.chromosome)?j->second.start:0x7fffffff;
      }
      else
      {
	intvl0.ethnicity[1]=j->second.ethnicity[0];
	end1=j->second.end;
      }
      intvl0.end=min(end0,end1);
      if (intvl0.ethnicity[0]>=0 || intvl0.ethnicity[1]>=0)
      {
        diploid[intvl0.index()]=intvl0;
	printf("%d [%d,%d] %s:%s\n",intvl0.chromosome,intvl0.start,intvl0.end,
	       (intvl0.ethnicity[0]<0)?"":ethnicities[intvl0.ethnicity[0]].c_str(),
	       (intvl0.ethnicity[1]<0)?"":ethnicities[intvl0.ethnicity[1]].c_str());
      }
      pos=intvl0.end;
      if (j->second.end<=pos || j->second.chromosome<chr)
	j++;
      if (i->second.end<=pos || i->second.chromosome<chr)
	i++;
    }
  }
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
