#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include <cfloat>
#include <set>
#include <queue>
#include <map>
#include <unordered_set>
#include <iomanip>
#include "global.h"
#include "recomb.h"
#include "individual.h"
#include "HypervolumeIndicator.h"
using namespace shark;
class CMOEAD
{

   public:
	CMOEAD();
	virtual ~CMOEAD();

	void init_population();                
	void replacement_phase();
	void evol_population();                                    
	// execute MOEAD
	void exec_emo(int run);
	void save_front(char savefilename[4024]); 
	void save_pos(char savefilename[4024]);
	void update_parameterD();
	double distance_var( vector<double> &a, vector<double> &b);
        void dominance_information(); 
	void full_dominance_information();
	void pick_penalized(unordered_set<int> &penalized, unordered_set<int> &survivors, vector<double> &distances);
	int max_HV_contribution(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors, unordered_set<int> &survivors_front);
	void dominance_information_remove(int rm_idx);
	void update_lowest_front(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors_front);
	void diversity_information(unordered_set<int> &survivors, unordered_set<int> &penalized, vector<double> &distances);
	void penalize_nearest(int idx_survivor, unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &penalized);
	void classic_hv_selection();
   private:
	vector <CIndividual> pool;
	vector<int> child_idx, parent_idx, inv_parent_idx;
        vector<int> Np,Rp;//rank
	vector<unordered_set<int> > Sp;//dominated indexes and inverse
        vector<vector<int> > fronts ;
  	HypervolumeIndicator m_indicator;

	// algorithm parameters
	long long nfes;          //  the number of function evluations
	double	D;	//Current minimum distance

};
CMOEAD::CMOEAD()
{

}
CMOEAD::~CMOEAD()
{

}
void CMOEAD::update_parameterD()
{
      double TElapsed = nfes, TEnd = max_nfes;
      D = Di - Di * (TElapsed / (TEnd*Df));
}
double CMOEAD::distance_var( vector<double> &a, vector<double> &b)
{
   double dist = 0 ;
   for(int i = 0; i < a.size(); i++)
   {
      double factor = (a[i]-b[i])/(vuppBound[i]-vlowBound[i]);
      dist += factor*factor;
   }
   return sqrt(dist);
}
void CMOEAD::init_population()
{
    // Load weight vectors
    for(int i=0; i< nPop+nOffspring; i++)
    {
	       CIndividual ind;
		// Randomize and evaluate solution
		ind.rnd_init();
		ind.obj_eval();
		// Initialize the reference point
		pool.push_back(ind); 
		if( i < nPop)
		{
		   parent_idx.push_back(i);
		}
		else
		   child_idx.push_back(i);
		nfes++;
     }
}
void CMOEAD::evol_population()
{
   dominance_information(); 
   for(int i = 0; i < nOffspring; )
   {
      int idx1=rand()%nPop, idx2=rand()%nPop, idx3=rand()%nPop, idx4=rand()%nPop;
      while(idx2==idx1) idx2=rand()%nPop;
      while(idx3==idx4) idx4=rand()%nPop;
      if(Rp[parent_idx[idx1]] > Rp[parent_idx[idx2]]) idx1=idx2;
      if(Rp[parent_idx[idx3]] > Rp[parent_idx[idx4]]) idx3=idx4;
      // produce a child solution
      CIndividual child1 = pool[parent_idx[idx1]], child2 = pool[parent_idx[idx2]];
      bool crosed = real_sbx_xoverA(pool[parent_idx[idx1]], pool[parent_idx[idx3]], child1, child2);
      // apply polynomial mutation
      bool mut1 = realmutation(child1, 1.0/nvar);
      bool mut2 = realmutation(child2, 1.0/nvar);
      if(crosed || mut1) nfes++;
      if(crosed || mut2) nfes++;
      child1.obj_eval();
      child2.obj_eval();
      swap(child1, child2);
      if(i < nOffspring){ //this modification allows to modify its behaviour from stady-state to traditional MOO
	pool[child_idx[i++]]=child1;
      }
      if(i < nOffspring){
	pool[child_idx[i++]]=child2;
      }
   }
   if(D>0)
     replacement_phase();
   else
     classic_hv_selection();
}
void CMOEAD::exec_emo(int run)
{
        char filename1[5024];
        char filename2[5024];
	seed = run;
	srand(run);
	seed = (seed + 23)%1377;
	rnd_uni_init = -(long)seed;

	//initialization
	nfes      = 0;
	init_population();

	sprintf(filename1,"%s/POS/POS_R2_EMOA_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_Px_%lf",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar), Df, CR);
	sprintf(filename2,"%s/POF/POF_R2_EMOA_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_Px_%lf",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar), Df, CR);
        long long current = nfes;
	long long accumulator = 0, bef = nfes;
	save_pos(filename1);
        save_front(filename2);
	while(nfes<max_nfes)
	{
		update_parameterD();
		evol_population();
		accumulator += nfes - bef ;
                if(accumulator > 0.1*(max_nfes)  )
		{
	           accumulator -= 0.1*(max_nfes);
		   save_pos(filename1);
		   save_front(filename2);
		}
		bef=nfes;
	        //nfes += nOffspring;
	}
	save_pos(filename1);
	save_front(filename2);
}
void CMOEAD::save_front(char saveFilename[4024])
{

    std::fstream fout;
    fout.open(saveFilename,fstream::app|fstream::out );
    for(int n=0; n < nPop; n++)
    {
       for(int k=0;k<nobj;k++)
          fout<<pool[parent_idx[n]].y_obj[k]<<"  ";
       fout<<"\n";
    }
    fout.close();
}

void CMOEAD::save_pos(char saveFilename[4024])
{
   std::fstream fout; //fout.open(saveFilename,std::ios::out);
   fout.open(saveFilename, fstream::app|fstream::out);
   for(int n=0; n<nPop; n++)
   {
      for(int k=0;k<nvar;k++)
         fout<<pool[parent_idx[n]].x_var[k] << "  ";
      for(int k=0;k<nvar;k++)
   	 fout<<pool[parent_idx[n]].x_var[k]/(vuppBound[k]-vlowBound[k]) << "  "; //fout<<population[n].indiv.x_var[k]<< fixed << setprecision(30) << "  ";
   	fout<<"\n";
   }
   fout.close();
}

void CMOEAD::diversity_information(unordered_set<int> &survivors, unordered_set<int> &candidates, vector<double> &distances)
{
   for(auto s_idx:survivors)
   {
      for(auto c_idx:candidates)
      {
	distances[c_idx] = min(distances[c_idx], distance_var(pool[s_idx].x_var, pool[c_idx].x_var));
      }
   } 
}
void CMOEAD::pick_penalized(unordered_set<int> &penalized, unordered_set<int> &survivors, vector<double> &distances)
{
   pair<double, int> max_spread = make_pair(-DBL_MAX, -1);
   for( auto p_idx : penalized)
   {
      if( distances[p_idx] > max_spread.first) max_spread = make_pair(distances[p_idx], p_idx);
   }
   survivors.insert(max_spread.second);
   penalized.erase(max_spread.second);
   for( auto p_idx : penalized)
   {
      distances[p_idx] = min(distances[p_idx], distance_var( pool[p_idx].x_var, pool[max_spread.second].x_var));
   }
}
void CMOEAD::update_lowest_front(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors_front)
{
  for(auto s_idx:survivors_front)
  {
	for(auto idx:Sp[s_idx])
        {
           Np[idx]--;
           if(Np[idx] == 0) candidates_front.insert(idx);
        }
        Np[s_idx]--;
  }
  survivors_front.clear();
}
void CMOEAD::replacement_phase()
{
  unordered_set<int> survivors, candidates, penalized, survivors_front, candidates_front;
  vector<double> distances((int)pool.size(), DBL_MAX);
  for(int idx = 0; idx < pool.size(); idx++) candidates.insert(idx); 
  dominance_information(); 
  for(auto c_idx:candidates) if( Np[c_idx] == 0) candidates_front.insert(c_idx);

  while( survivors.size() < nPop && !candidates.empty())
  {  
     //move to survivors..
     int new_survivor_idx = max_HV_contribution(candidates, candidates_front, survivors, survivors_front);
     //penalize nearest ind..
     penalize_nearest(new_survivor_idx, candidates, candidates_front, penalized);
     if(candidates.empty())break;
     //update ranks
     if(candidates_front.empty())
       update_lowest_front(candidates, candidates_front, survivors_front);
  }
  diversity_information(survivors, candidates, distances);
  while( survivors.size() < nPop ) pick_penalized(penalized, survivors, distances);
  int idx = 0; 
  vector<bool> setted( pool.size(), false);
  for(auto ite_s = survivors.begin() ; ite_s != survivors.end(); ite_s++, idx++)
  {
    parent_idx[idx] = *ite_s;
    setted[*ite_s] = true;
  }
  idx = 0;
  for(int i = 0; i < pool.size();  i++) if(!setted[i]) child_idx[idx++]=i;
}
void CMOEAD::penalize_nearest(int idx_survivor, unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &penalized)
{
     unordered_set<int> tmp=candidates;
    for(auto c_idx:tmp)
     {
	double dist = distance_var(pool[idx_survivor].x_var, pool[c_idx].x_var);
	if( dist < D ) 
        {
	  penalized.insert(c_idx);
          for(auto idx:Sp[c_idx]) Np[idx]--; 
          for(int idx=0; idx<Sp.size(); idx++) Sp[idx].erase(c_idx);
          candidates_front.erase(c_idx);
	  candidates.erase(c_idx);
       }		
     }
     for(auto c_idx:candidates) if(Np[c_idx]==0) candidates_front.insert(c_idx);
}
int CMOEAD::max_HV_contribution(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors, unordered_set<int> &survivors_front)
{
  
  vector<vector<double> > all(survivors_front.size()+1);
  int size=0;
  vector<double> ref(nobj, -DBL_MAX);
  for(auto idx:survivors_front){
	all[size++]=pool[idx].y_obj;
	for(int m = 0; m < nobj; m++)ref[m]=max(ref[m], pool[idx].y_obj[m]);
  }
  for(auto idx:candidates_front)
	for(int m = 0; m < nobj; m++)ref[m]=max(ref[m], pool[idx].y_obj[m]);
  for(int m = 0; m < nobj; m++)ref[m] +=1.0;
  m_indicator.setReference(ref);
  pair<double, int> max_ctr(-DBL_MAX, -1);

  for(auto idx:candidates_front)
  {
     all[size]=pool[idx].y_obj; 
     vector<pair<double, size_t> > res=m_indicator.lastIdxContributor(all, 1);
     double cur_ctr=res[0].first;
//     for(auto pp:res){
//	if( pp.second==size){
//	   cur_ctr=pp.first; break;
//	}
//     } 
     max_ctr = max(max_ctr, make_pair(cur_ctr, idx));
  }
  survivors_front.insert(max_ctr.second);
  survivors.insert(max_ctr.second);
  candidates.erase(max_ctr.second);
  candidates_front.erase(max_ctr.second);
  return max_ctr.second;
}
void CMOEAD::classic_hv_selection()
{
  full_dominance_information();
  vector<double>r;
  m_indicator.setReference(r);
  vector<int> new_idx_parent;
  int rank=0;
  while(new_idx_parent.size()+fronts[rank].size() < nPop)
  {
     for(auto idx:fronts[rank])new_idx_parent.push_back(idx);
     rank++;
  }
  vector<vector<double> > lastfront;
  for(auto idx:fronts[rank]) lastfront.push_back(pool[idx].y_obj);
  while(new_idx_parent.size()+fronts[rank].size() > nPop)
  {
    vector<pair<double, size_t> > to_remove=m_indicator.leastContributors(lastfront, 1);
    iter_swap(lastfront.begin()+to_remove[0].second, lastfront.end()-1);   
    iter_swap(fronts[rank].begin()+to_remove[0].second, fronts[rank].end()-1);   
    lastfront.pop_back();
    fronts[rank].pop_back();
  }
  for(auto idx:fronts[rank])new_idx_parent.push_back(idx);
  vector<bool> isparent(pool.size(), false);
  int i=0;
  for(auto idx:new_idx_parent) isparent[idx]=true, parent_idx[i++]=idx;
  for(int i = 0,j=0; i < pool.size(); i++) if(!isparent[i]) child_idx[j++]=i;
}
void CMOEAD::full_dominance_information()
{

   Sp.assign(nPop+nOffspring, unordered_set<int>());
   Np.assign(nPop+nOffspring, 0);
   Rp.assign(nPop+nOffspring, 0);
   fronts.assign(1, vector<int>());
   int rank = 0;
   for(int pidx1=0; pidx1< nPop+nOffspring; pidx1++)
   {
      for(int pidx2=0; pidx2< nPop+nOffspring; pidx2++)
      {
	if(pidx1 == pidx2) continue;
        if( pool[pidx1] < pool[pidx2]) Sp[pidx1].insert(pidx2);
 	else if( pool[pidx2] < pool[pidx1]) Np[pidx1]++;
      }
      if( Np[pidx1] == 0)
      {
         fronts[rank].push_back(pidx1);
	 Rp[pidx1]=rank;
      }
   }
   while(true)
   {
      vector<int> next_front;
      for(auto idx:fronts[rank])
      {
	for(auto idx_dominated:Sp[idx])
        {
	  Np[idx_dominated]--;
          if(Np[idx_dominated]  == 0) 
          {
	     next_front.push_back(idx_dominated);
	     Rp[idx_dominated] = rank+1;
          }
        }
      }
      if(next_front.empty()) break;
      fronts.push_back(next_front);
      rank++;
   }
}
void CMOEAD::dominance_information()
{
   Sp.assign(nPop+nOffspring, unordered_set<int>());
   Np.assign(nPop+nOffspring, 0);
   Rp.assign(nPop+nOffspring, 0);
   fronts.assign(1, vector<int>());
   int rank = 0;
   for(int pidx1=0; pidx1< nPop+nOffspring; pidx1++)
   {
      for(int pidx2=0; pidx2< nPop+nOffspring; pidx2++)
      {
	if(pidx1 == pidx2) continue;
        if( pool[pidx1] < pool[pidx2]) Sp[pidx1].insert(pidx2);
 	else if( pool[pidx2] < pool[pidx1]) Np[pidx1]++;
      }
      if( Np[pidx1] == 0)
      {
         fronts[rank].push_back(pidx1);
	 Rp[pidx1]=rank;
      }
   }
}
#endif
