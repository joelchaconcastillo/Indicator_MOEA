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
#include "common.h"
#include "individual.h"

//#define EXPAND(n)((long long)(n*1e6))
#define EXPAND(n)(n)
class CMOEAD
{

   public:
	CMOEAD();
	virtual ~CMOEAD();

	void init_population();                
	bool update_reference(CIndividual &ind); 
	void replacement_phase();
	void evol_population();                                    
	// execute MOEAD
	void exec_emo(int run);
	void save_front(char savefilename[4024]); 
	void save_pos(char savefilename[4024]);
	void update_parameterD();
	double distance_var( vector<double> &a, vector<double> &b);
        void dominance_information(); 
	void pick_penalized(unordered_set<int> &penalized, unordered_set<int> &survivors, vector<double> &distances);
	int max_R2_contribution(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors, unordered_set<int> &survivors_front, unordered_set<int> &penalized, vector<set<pair<double, int> > > &survivors_weight);
	void dominance_information_remove(int rm_idx);
	void update_lowest_front(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors_front);
	void diversity_information(unordered_set<int> &survivors, unordered_set<int> &penalized, vector<double> &distances);
	void penalize_nearest(int idx_survivor, unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &penalized);
        void table_fitness_information();
   private:
	vector <CIndividual> pool;
	vector<int> child_idx, parent_idx, inv_parent_idx;
	vector<vector<double> > namda;     // weight vector
        vector<int> Np,Rp;//rank
	vector<unordered_set<int> > Sp;//dominated indexes and inverse
        vector<vector<int> > fronts ;
        vector<vector<double > > table_fitness;
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
//      D = max(D, 0.0);
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
    idealpoint = vector<double>(nobj, 1.0e+30);
    char filename[1024];
    // Read weight vectors from a data file
    sprintf(filename,"%s/ParameterSetting/Weight/W%dD_%d.dat", strpath, nobj, nWeight);
    std::ifstream readf(filename);
    namda.resize(nWeight, vector<double> (nobj, 0.0));
    table_fitness.assign(nWeight, vector<double> (nPop+nOffspring));
    // Load weight vectors
    for(int i=0; i< nWeight; i++)
	for(int j=0; j<nobj; j++)
	 readf>>namda[i][j];

    for(int i=0; i< nPop+nOffspring; i++)
    {
	       CIndividual ind;
		// Randomize and evaluate solution
		ind.rnd_init();
		ind.obj_eval();
		// Initialize the reference point
		update_reference(ind);
		// Save in the population
		pool.push_back(ind); 
		if( i < nPop)
		{
		   parent_idx.push_back(i);
		}
		else
		   child_idx.push_back(i);
		nfes++;
     }
     readf.close( );
}
bool CMOEAD::update_reference(CIndividual &ind)
{
   bool changed = false;
   for(int n=0; n<nobj; n++)
   {
      if(ind.y_obj[n]<idealpoint[n])
      {
         idealpoint[n] = ind.y_obj[n];
         changed = true;
      }
   }
   return changed;
}
void CMOEAD::evol_population()
{
   dominance_information(); 
   for(int i = 0; i < nOffspring; i +=2)
   {
      int idx1=rand()%nPop, idx2=rand()%nPop, idx3=rand()%nPop, idx4=rand()%nPop;
      while(idx2==idx1) idx2=rand()%nPop;
      while(idx3==idx4) idx4=rand()%nPop;
      if(Rp[parent_idx[idx1]] > Rp[parent_idx[idx2]]) idx1=idx2;
      if(Rp[parent_idx[idx3]] > Rp[parent_idx[idx4]]) idx3=idx4;
      // produce a child solution
      CIndividual &child1 = pool[child_idx[i]], &child2 = pool[child_idx[i+1]];
      
      bool crosed = real_sbx_xoverA(pool[parent_idx[idx1]], pool[parent_idx[idx3]], child1, child2);
      // apply polynomial mutation
      bool mut1 = realmutation(child1, 1.0/nvar);
      bool mut2 = realmutation(child2, 1.0/nvar);
      if(crosed || mut1) nfes++;
      if(crosed || mut2) nfes++;
      child1.obj_eval();
      child2.obj_eval();
      // update the reference points and other solutions in the neighborhood or the whole population
      update_reference(child1); //O(M)
      update_reference(child2); //O(M)
   }
   replacement_phase();
}
void CMOEAD::exec_emo(int run)
{
        char filename1[5024];
        char filename2[5024];
	seed = run;
	seed = (seed + 23)%1377;
	rnd_uni_init = -(long)seed;

	//initialization
	nfes      = 0;
	init_population();

	sprintf(filename1,"%s/POS/POS_R2_EMOA_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_CR_%lf_F_%lf",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar), Df, CR, F);
	sprintf(filename2,"%s/POF/POF_R2_EMOA_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_CR_%lf_F_%lf",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar), Df, CR, F);
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
void CMOEAD::table_fitness_information()
{
  //getting R2-contribution of each individual...
  for(int w_id = 0; w_id < nWeight; w_id++)
  {
     for(int idx = 0; idx < pool.size(); idx++)
     {
        table_fitness[w_id][idx] = fitnessfunction(pool[idx].y_obj, namda[w_id]);
     }
  }
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
  survivors_front.clear(); //reset survivors front, this is important and is related to the R2-Contribution, each front requires a different computation of R2-indicator
}
int CMOEAD::max_R2_contribution(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors, unordered_set<int> &survivors_front, unordered_set<int> &penalized, vector<set<pair<double, int> > > &survivors_weight)
{
  pair<double, int> max_contribution(-DBL_MAX, -1);
  if(survivors_front.empty())
  {
     survivors_weight.assign(nWeight, set<pair<double, int>>());
     for(auto c_idx:candidates_front)
     {
	double sum = 0;
        for(int w_id=0; w_id < nWeight; w_id++)
        {
           sum -= table_fitness[w_id][c_idx];
        }
	if(max_contribution.first < sum) max_contribution = make_pair(sum, c_idx);
     }
    max_contribution.first *=-1; //this is not necessary..
  }
  else
  {
     for(auto c_f:candidates_front)
     {
       double Total_contribution = 0.0;
       for(int w_id = 0; w_id < nWeight; w_id++)
          Total_contribution += max(0.0, survivors_weight[w_id].begin()->first - table_fitness[w_id][c_f]);
       if(Total_contribution > max_contribution.first) max_contribution = make_pair(Total_contribution, c_f);
     }
  }
  survivors_front.insert(max_contribution.second);
  survivors.insert(max_contribution.second);
  candidates.erase(max_contribution.second);
  candidates_front.erase(max_contribution.second);
  for(int w_id=0; w_id < nWeight; w_id++)  survivors_weight[w_id].insert(make_pair(table_fitness[w_id][max_contribution.second], max_contribution.second));
  return max_contribution.second;
}
void CMOEAD::replacement_phase()
{
  unordered_set<int> survivors, candidates, penalized, survivors_front, candidates_front;
  vector<double> distances((int)pool.size(), DBL_MAX);
  vector<set<pair<double, int> > > survivors_weight;
  for(int idx = 0; idx < pool.size(); idx++) candidates.insert(idx); 
  table_fitness_information();
  dominance_information(); 
  for(auto c_idx:candidates) if( Np[c_idx] == 0) candidates_front.insert(c_idx);

  while( survivors.size() < nPop && !candidates.empty())
  {  
     //move to survivors..
     int new_survivor_idx = max_R2_contribution(candidates, candidates_front, survivors, survivors_front, penalized, survivors_weight);
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
#endif
