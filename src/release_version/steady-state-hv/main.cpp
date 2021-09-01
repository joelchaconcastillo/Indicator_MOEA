#include <bits/stdc++.h>
#include "steady_state_moea.hpp"
using namespace std;
int main()
{
  srand(time(0));
  max_nfes=2500000;
  nvar=24;
  nobj=2;
  npop=100;
  Pm=1.0/nvar;
  param_k=4;
  run();
  return 0;
}
