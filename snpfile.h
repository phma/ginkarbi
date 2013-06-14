/* Binary file format for snp files:
 * One file for each person contains all that person's SNPs (just the AT GG DI ... parts of the line);
 * one file shared by all people contains the names and locations of all SNPs.
 * Header common to both files:
 * 6a 67 69 6e 61 30 : Magic string. SNP name file ends in 30; individual SNP file ends in 31.
 * 02 : Domain Eukarya. Bacteria are 00, Archaea are 01.
 * 00 : Kingdom Animalia. Plants are 80. Bikonts are 00-7f; unikonts are 80-ff. Division of bacteria to be defined.
 * 48 6f 6d 6f 00 : Genus (or genus and species) of organism. Kingdom is required to disambiguate e.g. Aotus.
 * 00 25 00 68 : Build 37, annotation 104.
 * 23 20 54 68 ... 00 : Null-terminated comment.
 * 
 * Name and location file:
 * 00 44 51 1c 00 01 00 01 40 ea : rs4477212, chromosome 1, position 82154
 * Chromosome takes 2 bytes because there are plants with over 1000 chromosomes.
 * Position is unsigned and can be as much as 249210707 (0xedaa753) for humans
 * or about 3.75G (0xdf800000) for some organism.
 * 
 * Personal SNP file:
 * 0b 0d 0d 0d 15 10 0c 10 ... : AA AG AG AG GG CC AC CC ...
 * Each SNP is stored in one byte, 00-1a, representing these calls in order:
 * -- D  I  DD DI ID II A  C 
 * G  T  AA AC AG AT CA CC CG
 * CT GA GC GG GT TA TC TG TT
 * The upper three bits are used internally to store uncertainty during phasing.
 */
#include <string>
#include <vector>
using namespace std;

extern const string calls;
struct snp
{
  int snpname; // one bit tells whether it's rs or i, and the rest is a number
  int chromosome;
  int position;
  int allele;
  int ethnicity[2];
  snp();
};

vector<snp> readgenometextfile(char *genomename);
string snpstring(int snpnum);
string chromstring(int chromnum);
string posstring(int posnum);
string allelestring(int allelenum);
