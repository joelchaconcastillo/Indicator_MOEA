#include <bits/stdc++.h>
#include "HypervolumeIndicator.h"
using namespace std;
using namespace shark;
int main()
{
//  IndicatorBasedSelection<HypervolumeIndicator> m_selection; ///< Selection operator relying on the (contributing) hypervolume indicator.
  HypervolumeIndicator m_indicator;
  vector<vector<double> > p;
  p.push_back({0.0, 0.0, 0.5});
  p.push_back({0.0, 0.5, 0.0});
  p.push_back({0.5, 0.0, 0.0});
  vector<double> ref={1.0, 1.0, 1.0};
  m_indicator.setReference(ref);
  auto res=m_indicator.leastContributors(p,1);
  cout << res[0].first<<endl;
 
  return 0;
}
