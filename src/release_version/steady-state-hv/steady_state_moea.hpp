#ifndef STEADY_STATE_MOEA_HPP
#define STEADY_STATE_MOEA_HPP
#include "HypervolumeIndicator.h"
#include "problem.h"
#define EPS 1.0e-240
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
vector<unordered_set<int> > fronts;
vector<int> Np, Rp;
vector<unordered_set<int> > Sp;
vector<bool> prev_survivor;
bool operator<(const vector<double> &ind1, const vector<double> &ind2)
{
    bool dominated = true;
    for(int n=0; n<nobj; n++)
        if(ind2[n]<ind1[n]) return false;
    if(ind2==ind1) return false;
    return dominated;
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
//    vector<double> child2=child1;
    bool changed=false;
	changed=true;
    double rand, y1, y2, yl, yu, c1, c2, alpha, beta, betaq;
    if (rnd_uni() <= Px){
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
                    //    child2[i] = c1;
                    }
                    else{
                        child1[i] = c1;
                      //  child2[i] = c2;
                    }
                }
                else{
                    child1[i] = parent1[i];
                    //child2[i] = parent2[i];
                }
            }
            else{
                child1[i] = parent1[i];
               // child2[i] = parent2[i];
            }
        }
    }
    else child1 = parent1;//, child2 = parent2;

    return changed;
}


void updateDt(){
      Dt = Di - Di * (nfes/ (double)(max_nfes*Df));
}
void dominance_info(){
   Sp.assign(npop, unordered_set<int>());
   Np.assign(npop, 0);
   Rp.assign(npop+1, 0);
   fronts.assign(npop+1, unordered_set<int>());
   for(int pidx1=0; pidx1< npop; pidx1++){
      for(int pidx2=0; pidx2< npop; pidx2++){
	if(pidx1 == pidx2) continue;
        if( parent_yobj[pidx1] < parent_yobj[pidx2]) Sp[pidx1].insert(pidx2);
 	else if( parent_yobj[pidx2] < parent_yobj[pidx1]) Np[pidx1]++;
      }
      if( Np[pidx1]==0)fronts[0].insert(pidx1), Rp[pidx1]=0;
   }    
   int rank = 0;
   vector<int> innerNp=Np;
   while(true){
     for(auto idx:fronts[rank]){
       for(auto idx_dom:Sp[idx]){
	  innerNp[idx_dom]--;
	  if(innerNp[idx_dom]==0)
	     fronts[rank+1].insert(idx_dom), Rp[idx_dom]=rank+1;
       }
     }
     if(fronts[rank+1].empty())break;
     rank++;
   }
}
int update_fronts(){
  Rp.back()=0;
  int lastRank=0;
  for(int i = 0; i < npop; i++){
     if( parent_yobj[i]<parent_yobj.back()){
       Rp.back()=max(Rp[i]+1, Rp.back());
     }
     if(parent_yobj.back()< parent_yobj[i]){
	 fronts[Rp[i]].erase(i);
         Rp[i]++;
	 fronts[Rp[i]].insert(i);
     }
    lastRank=max(lastRank, Rp[i]);
  }
  fronts[Rp.back()].insert(npop);
  lastRank=max(lastRank, Rp.back());
  return lastRank;
}
void classic_hv_selection(){
  int lastRank = update_fronts();
  vector<vector<double> > lastfront;
  vector<int> idxs;
  for(auto idx:fronts[lastRank]) lastfront.push_back(parent_yobj[idx]), idxs.push_back(idx);
  vector<pair<double, size_t> > to_remove=m_indicator.leastContributors(lastfront, 1);
  int id=idxs[to_remove[0].second];
  iter_swap(parent_yobj.begin()+id, parent_yobj.end()-1);
  iter_swap(parent_xvar.begin()+id, parent_xvar.end()-1);
  fronts[Rp[id]].erase(id);
  fronts[Rp[npop]].erase(npop);
  fronts[Rp[npop]].insert(id);
  iter_swap(Rp.begin()+id, Rp.end()-1);
}
void replacement(){
   int lastRank = update_fronts();
   vector<vector<double> > survivors, candidates, penalized;
   vector<int> idx_survivors, idx_candidates, idx_penalized;  
   for(int i = 0; i < parent_yobj.size(); i++){ //first optimization trick
       if(Rp[i] < Rp.back() && prev_survivor[i])
	   survivors.push_back(parent_yobj[i]), idx_survivors.push_back(i);
       else
	   candidates.push_back(parent_yobj[i]), idx_candidates.push_back(i);
       prev_survivor[i]=false;
   }
   for(int i = 0; i < survivors.size(); i++){
      for(int j = 0; j < candidates.size(); j++){
       double dist=distance_var(survivors[i], candidates[j]);
         if(dist < Dt){
	     penalized.push_back(survivors[i]);
	     idx_penalized.push_back(i);
	 }
      }
   } 


//   while(survivors.size() < npop && !candidates.empty()){
//
//   }
   
   for(auto idx:idx_survivors) prev_survivor[idx]=true;
   vector<double> mdist(npop+1, DBL_MAX);
      for(int i = 0; i < penalized.size(); i++){
          for(int j = 0; j < survivors.size(); j++){
	    mdist[i] = min( mdist[i], distance_var(penalized[i], survivors[j]));
      }
   }
   while(survivors.size() < npop){
      pair<double, int> maxdcn(-1, -1);
      for(int i = 0; i < penalized.size(); i++){
	      maxdcn=max(maxdcn, make_pair(mdist[i], i) );
      } 
      survivors.push_back(penalized[maxdcn.second]); 
      idx_survivors.push_back(idx_penalized[maxdcn.second]);
      iter_swap(penalized.begin()+maxdcn.second, penalized.end()-1);
      iter_swap(idx_penalized.begin()+maxdcn.second, idx_penalized.end()-1);
      penalized.pop_back();
      idx_penalized.pop_back();
   }
}

void eval(vector<double> &yobj, vector<double> &xvar){
   //wfg8(yobj, xvar, param_k);
   //wfg1(yobj, xvar, param_k);
   dtlz1(yobj, xvar);
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

void evol(){
   updateDt();
   pair<int, int> idxs=binary_torunament();
   real_sbx_xoverA(parent_xvar[idxs.first], parent_xvar[idxs.second], parent_xvar.back());
   realmutation(parent_xvar.back());
   eval(parent_yobj.back(), parent_xvar.back()); 
   replacement();
   //classic_hv_selection();
}
void init(){
   prev_survivor.assign(npop+1, false);
   //for(int i = 0; i < nvar; i++) lb[i]=0.0, ub[i]=(i+1.0)*2.0;
   for(int i = 0; i < nvar; i++) lb[i]=0.0, ub[i]=1.0;
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
  for(int i1=0; i1<npop+1; i1++){
    for(int i2=0; i2<nobj; i2++)  cout << parent_yobj[i1][i2]<<" "; 
     cout << Rp[i1];
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
