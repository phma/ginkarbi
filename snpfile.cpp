#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include "snpfile.h"

const string calls="--D I DDDIIDIIA C G T AAACAGATCACCCGCTGAGCGGGTTATCTGTT  ";
const string whitespace=" \n\r\t\f\v";

string linebuffer;
char bytebuffer[64];

snp::snp()
{
  snpname=chromosome=position=allele=0;
  ethnicity[0]=ethnicity[1]=-1;
}

int snpnumber(string snpstr)
{
  int rsi,num,i;
  if (snpstr.substr(0,1)=="i")
  {
    rsi=0x80000000;
    snpstr.erase(0,1);
  }
  else if (snpstr.substr(0,2)=="rs")
  {
    rsi=0;
    snpstr.erase(0,2);
  }
  else
    rsi=-1;
  for (i=num=0;i<snpstr.length();i++)
    num=10*num+snpstr[i]-'0';
  return rsi|num;
}

string snpstring(int snpnum)
{
  char buf[24];
  memset(buf,0,24);
  sprintf(buf,"%s%d",(snpnum&0x80000000)?"i":"rs",snpnum&0x7fffffff);
  return string(buf);
}

int chromnumber(string chromstr)
{
  int num,i;
  if (chromstr=="X")
    num=23;
  else if (chromstr=="Y")
    num=24;
  else if (chromstr=="MT")
    num=0;
  else
    for (i=num=0;i<chromstr.length();i++)
      num=10*num+chromstr[i]-'0';
  return num;
}

string chromstring(int chromnum)
{
  char buf[24];
  memset(buf,0,24);
  sprintf(buf,"%d",chromnum);
  if (chromnum==0)
  {
    buf[0]='M';
    buf[1]='T';
  }
  if (chromnum>22)
  {
    buf[0]=chromnum+'X'-23;
    buf[1]=0;
  }
  return string(buf);
}

int allelenumber(string allelestr)
{
  int i,num;
  while (allelestr.length()<2)
    allelestr+=' ';
  for (i=num=0;i<calls.length()/2;i++)
    if (allelestr==calls.substr(2*i,2))
      num=i;
  return num%27;
}

string allelestring(int allelenum)
{
  return calls.substr(2*allelenum,2);
}

int posnumber(string posstr)
{
  int num,i;
  for (i=num=0;i<posstr.length();i++)
    num=10*num+posstr[i]-'0';
  return num;
}

string posstring(int posnum)
{
  char buf[24];
  memset(buf,0,24);
  sprintf(buf,"%d",posnum);
  return string(buf);
}

vector<string> words(string line)
{
  string word;
  vector<string> wordlist;
  size_t whpos;
  while (line.length())
  {
    whpos=line.find_first_of(whitespace);
    if (whpos!=string::npos)
      whpos++;
    word=line.substr(0,whpos-1);
    if (word.length())
      wordlist.push_back(word);
    line.erase(0,whpos);
  }
  return wordlist;
}

snp readgenometextline(FILE *genomefile)
{
  string line;
  size_t nread,eolpos;
  snp readsnp;
  vector<string> wordlist;
  while (linebuffer.find_first_of("\r\n")==string::npos && !feof(genomefile) && !ferror(genomefile))
  {
    nread=fread(bytebuffer,1,64,genomefile);
    linebuffer.append(bytebuffer,nread);
  }
  eolpos=linebuffer.find_first_of("\r\n");
  if (eolpos!=string::npos)
    eolpos++;
  line=linebuffer.substr(0,eolpos-1);
  linebuffer.erase(0,eolpos);
  eolpos=line.find_first_of("#");
  if (eolpos!=string::npos)
    line.erase(eolpos); //ignore comments
  wordlist=words(line);
  if (wordlist.size()>=4)
  {
    readsnp.snpname=snpnumber(wordlist[0]);
    readsnp.chromosome=chromnumber(wordlist[1]);
    readsnp.position=posnumber(wordlist[2]);
    readsnp.allele=allelenumber(wordlist[3]);
  }
  return readsnp;
}

vector<snp> readgenometextfile(char *genomename)
{
  FILE *genomefile;
  snp readsnp;
  vector<snp> genome;
  genomefile=fopen(genomename,"r");
  while ((!feof(genomefile) && !ferror(genomefile)) || linebuffer.length())
  {
    readsnp=readgenometextline(genomefile);
    if (readsnp.snpname)
    {
      //cout<<snpstring(readsnp.snpname)<<"\t"<<chromstring(readsnp.chromosome)<<"\t"<<posstring(readsnp.position)<<"\t"<<allelestring(readsnp.allele)<<endl;
      genome.push_back(readsnp);
    }
  }
  return genome;
}
