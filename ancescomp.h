#include <string>
using namespace std;

/* Bitmap for alleles:
 * A 000001
 * C 000010
 * G 000100
 * T 001000
 * I 010000
 * D 100000
 * Stored as 18 bits, the top six being the exclusive-or of the middle and bottom.
 */
struct snp
{
  string snpname;
  int chromosome;
  int position;
  int allele;
  string ethnicity[2];
};

struct interval
{
  int chromosome;
  int start,end;
  string ethnicity[2];
};
