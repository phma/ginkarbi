#include <cstdio>
#include <vector>

/* List of ethnicities. root includes all humans; every subcategory
 * is after the category it is part of.
 */
#include "ethnicity.h"

vector<string> ethnicities;

void init_ethnic()
{
  // Level 0
  ethnicities.push_back("root");
  // Level 1
  ethnicities.push_back("eana"); // East Asian Native American
  ethnicities.push_back("mena"); // Middle Eastern North African
  ethnicities.push_back("european");
  ethnicities.push_back("african"); // south of the Sahara
  ethnicities.push_back("south_asian");
  ethnicities.push_back("oceanian");
  // Level 2
  ethnicities.push_back("east_asian");
  ethnicities.push_back("american");
  ethnicities.push_back("middle_eastern");
  ethnicities.push_back("north_african");
  ethnicities.push_back("north_european");
  ethnicities.push_back("south_european");
  ethnicities.push_back("east_european");
  ethnicities.push_back("ashkenazi");
  // Level 3
  ethnicities.push_back("british_irish");
  ethnicities.push_back("french_german");
  ethnicities.push_back("scandinavian");
  ethnicities.push_back("finnish");
  ethnicities.push_back("iberian");
  ethnicities.push_back("italian");
  ethnicities.push_back("balkan");
  ethnicities.push_back("sardinian");
  //ethnicities.push_back("");
}

int find_ethnic(string eth)
{
  int i;
  for (i=0;i<ethnicities.size();i++)
    if (ethnicities[i]==eth)
      return i;
  fprintf(stderr,"unknown ethnicity: %s\n",eth.c_str());
  return 0;
}

string ethstr(int eth)
{
  if (eth<0)
    return "";
  else
    return ethnicities[eth];
}
