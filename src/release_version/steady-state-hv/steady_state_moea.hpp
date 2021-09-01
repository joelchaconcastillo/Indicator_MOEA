#ifndef STEADY_STATE_MOEA_HPP
#define STEADY_STATE_MOEA_HPP
#include "HypervolumeIndicator.h"
#include "problem.h"
#define EPS 1.0e-14
using namespace std;
using namespace shark;
double rnd_uni(){
  return (double)rand() / (double)((unsigned)RAND_MAX + 1) ;
}
//configuration params...
char strpath[800];
double lb[1000], ub[1000];
int etax=2, etam=50, seed=1, nvar, nobj, npop, param_k;
double Di=0.4, Df=0.5, Px=0.4, Pm, Dt;
long long max_nfes, nfes=0;
HypervolumeIndicator m_indicator;
vector<vector<double> > parent_xvar, parent_yobj;
vector<vector<int> > fronts;
vector<int> Np, Rp;
vector<unordered_set<int> > Sp;
bool operator<(const vector<double> &ind1, const vector<double> &ind2)
{
    bool dominated = true;
    for(int n=0; n<nobj; n++)
        if(ind2[n]<ind1[n]) return false;
    if(ind2==ind1) return false;
    return dominated;
}
void dominance_info(){
   Sp.assign(npop+1, unordered_set<int>());
   Np.assign(npop+1, 0);
   Rp.assign(npop+1, 0);
   fronts.assign(1, vector<int>());
   for(int pidx1=0; pidx1< npop; pidx1++){
      for(int pidx2=0; pidx2< npop; pidx2++){
	if(pidx1 == pidx2) continue;
        if( parent_yobj[pidx1] < parent_yobj[pidx2]) Sp[pidx1].insert(pidx2);
 	else if( parent_yobj[pidx2] < parent_yobj[pidx1]) Np[pidx1]++;
      }
      if( Np[pidx1]==0)fronts[0].push_back(pidx1);
   }    
   int rank = 0;
   vector<int> innerNp=Np;
   while(true){
     vector<int> nxt;
     for(auto idx:fronts[rank]){
       for(auto idx_dom:Sp[idx]){
	  innerNp[idx_dom]--;
	  if(innerNp[idx_dom]==0)
	     nxt.push_back(idx_dom), Rp[idx_dom]=rank+1;
       }
     }
     if(nxt.empty())break;
     fronts.push_back(nxt);
     rank++;
   }
}
void classic_hv_selection()
{
//   for(int i = 0; i < npop; i++){
//      if( parent_yobj.back() < parent_yobj[i]) Sp[npop].insert(i);
//      else if(parent_yobj[i] < parent_yobj.back()) Np[npop]++, Sp[i].insert(npop); 
//   }
  npop++;
  dominance_info();
  npop--;
  vector<int> new_idx_parent;
  int rank=0;
  while(new_idx_parent.size()+fronts[rank].size() < npop)
  {
     for(auto idx:fronts[rank])new_idx_parent.push_back(idx);
     rank++;
  }
  vector<vector<double> > lastfront;
  for(auto idx:fronts[rank]) lastfront.push_back(parent_yobj[idx]);
  vector<pair<double, size_t> > to_remove=m_indicator.leastContributors(lastfront, 1);
  iter_swap(parent_yobj.begin()+fronts[rank][to_remove[0].second], parent_yobj.end()-1);
  iter_swap(parent_xvar.begin()+fronts[rank][to_remove[0].second], parent_xvar.end()-1);
}
void updateDt(){
      Dt = Di - Di * (nfes/ (double)(max_nfes*Df));
}
double distance_var( vector<double> &xa, vector<double> &xb)
{
   double dist = 0;
   for(int i = 0; i < nvar; i++)
      dist += ((xa[i]-xb[i])*(xa[i]-xb[i])) / ( (ub[i]-lb[i])*(ub[i]-lb[i]));
   return sqrt(dist);
}
bool realmutation(vector<double> &x){
    double rnd, delta1, delta2, mut_pow, deltaq, y, yl, yu, val, xy;
    bool chaged=false;
    for (int j=0; j<nvar; j++)
    {
        if (rnd_uni()<=Pm)
        {
	    chaged=true;
            y  = x[j];
            yl = lb[j];
            yu = ub[j];
            delta1 = (y-yl)/(yu-yl);
            delta2 = (yu-y)/(yu-yl);
            rnd = rnd_uni();
            mut_pow = 1.0/(etam+1.0);
            if (rnd <= 0.5)
            {
                xy = 1.0-delta1;
                val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(etam+1.0)));
                deltaq =  pow(val,mut_pow) - 1.0;
            }
            else
            {
                xy = 1.0-delta2;
                val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(etam+1.0)));
                deltaq = 1.0 - (pow(val,mut_pow));
            }
            y = y + deltaq*(yu-yl);
	    y = max(y, yl);
	    y = min(y, yu);
            x[j] = y;
        }
    }
    return chaged;
}
bool real_sbx_xoverA(vector<double> &parent1, vector<double> &parent2, vector<double> &child1)
{
    vector<double> child2=child1;
    bool changed=false;
    double rand, y1, y2, yl, yu, c1, c2, alpha, beta, betaq;
    if (rnd_uni() <= Px){
	changed=true;
        for (int i=0; i<nvar; i++){
            if (rnd_uni()<=0.5 ){
                if (fabs(parent1[i]-parent2[i]) > EPS){
                    if (parent1[i] < parent2[i]){
                        y1 = parent1[i];
                        y2 = parent2[i];
                    }
                    else{
                        y1 = parent2[i];
                        y2 = parent1[i];
                    }
                    yl = lb[i];
                    yu = ub[i];
                    rand = rnd_uni();
                    beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(etax+1.0));
                    if (rand <= (1.0/alpha))
                        betaq = pow ((rand*alpha),(1.0/(etax+1.0)));
                    else
                        betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(etax+1.0)));
                    c1 = 0.5*((y1+y2)-betaq*(y2-y1));
                    beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(etax+1.0));
                    if (rand <= (1.0/alpha))
                        betaq = pow ((rand*alpha),(1.0/(etax+1.0)));
                    else
                        betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(etax+1.0)));
                    c2 = 0.5*((y1+y2)+betaq*(y2-y1));
		    c1 = max(c1, yl);
		    c2 = max(c2, yl);
                    c1 = min(c1, yu);
                    c2 = min(c2, yu);
                    if (rnd_uni()<=0.5){
                        child1[i] = c2;
                        child2[i] = c1;
                    }
                    else{
                        child1[i] = c1;
                        child2[i] = c2;
                    }
                }
                else{
                    child1[i] = parent1[i];
                    child2[i] = parent2[i];
                }
            }
            else{
                child1[i] = parent1[i];
                child2[i] = parent2[i];
            }
        }
    }
    else child1 = parent1, child2 = parent2;
    return changed;
}

void eval(vector<double> &yobj, vector<double> &xvar){
   wfg8(yobj, xvar, param_k);
   //dtlz1(yobj, xvar);
}
pair<int, int> binary_torunament(){
  int ps[2];
  for(int i = 0; i < 2; i++){
    int idx1=rand()%npop, idx2=rand()%npop;
     while(idx1==idx2) idx2=rand()%npop;
     if(Rp[idx1]<Rp[idx2]) ps[i]=idx1;
     else if(Rp[idx1]>Rp[idx2]) ps[i]=idx2;
     else ps[i]=(rand()%2)?idx1:idx2;
  }
  return make_pair(ps[0], ps[1]); 
}
void replacement(){
   
}
void evol(){
   updateDt();
   pair<int, int> idxs=binary_torunament();
   real_sbx_xoverA(parent_xvar[idxs.first], parent_xvar[idxs.second], parent_xvar.back());
   realmutation(parent_xvar.back());
   eval(parent_yobj.back(), parent_xvar.back()); 
   //replacement();
   classic_hv_selection();
}
void init(){
   for(int i = 0; i < nvar; i++) lb[i]=0.0, ub[i]=(i+1.0)*2.0;
   parent_xvar.assign(npop+1, vector<double>(nvar));
   parent_yobj.assign(npop+1, vector<double>(nobj));
   nfes=0; 
   for(int i = 0; i < npop; i++, nfes++){
     for(int j = 0; j < nvar; j++)
	parent_xvar[i][j] = lb[j] + rnd_uni()*(ub[j]-lb[j]); 
     eval(parent_yobj[i], parent_xvar[i]);
   }
   dominance_info();
   updateDt();
}
void end(){
  for(auto i1:parent_yobj){
    for(auto i2:i1)  cout << i2<<" ";
     cout <<endl;
  }
}
void run(){
   init();
   while(nfes<max_nfes){
     evol();
     nfes++;
   }
   end();
}
#endif
