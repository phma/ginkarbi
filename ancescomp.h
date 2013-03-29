#include <string>
#include <stdint.h>
using namespace std;

struct token
{
  int ch;
  int n;
  string str;
};

struct interval
{
  int chromosome;
  int start,end;
  int ethnicity[2];
  int64_t index();
  void clear();
};
