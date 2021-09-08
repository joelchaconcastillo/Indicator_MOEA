#include <bits/stdc++.h>
#include "steady_state_moea.hpp"
using namespace std;
int main()
{
  srand(1);
  max_nfes=10000;
  nvar=2;
  nobj=2;
  npop=100;
  Pm=1.0/nvar;
  param_k=1;
  run();
  return 0;
}
