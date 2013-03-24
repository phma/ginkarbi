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
{
  token tok;
  ancestryfile=fopen(ancestryname,"r");
  if (ancestryfile)
  {
    while (tok.ch!=EOF)
    {
      tok=readtoken();
      printf("%c\t%d\t%s\n",tok.ch,tok.n,tok.str.c_str());
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
