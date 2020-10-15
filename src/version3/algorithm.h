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
	void penalization(unordered_set<int> &candidates, unordered_set<int> &penalized, vector<double> &distances, vector<set<pair<double, int> > > &ranking_weights);
	void pick_best_R2(unordered_set<int> &survivors, unordered_set<int> &candidates, unordered_set<int> &penalized, vector<double> &distances, vector<set<pair<double, int> > > &ranking_weights);
	void pick_penalized(unordered_set<int> &penalized, unordered_set<int> &survivors, vector<double> &distances);
	void diversity_information(unordered_set<int> &survivors, unordered_set<int> &penalized, vector<double> &distances);

        void table_fitness_information();
   private:
	vector <CIndividual> pool;
	vector<int> child_idx, parent_idx, inv_parent_idx;
	vector<vector<double> > namda;     // weight vector
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
      D = max(D, 0.0);
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
   for(int i = 0; i < nOffspring; i++)
   {
      int idx_target = rand()%nPop;
      int idx1=rand()% nPop, idx2=rand()%nPop, idx3=rand()%nPop;
      while(idx1 == idx_target) idx1=rand()%nPop;
      //while(idx2 == idx1 ) idx2=rand()%nPop;
      while(idx2 == idx1 || idx2 == idx_target) idx2=rand()%nPop;
      //while(idx3 == idx2 || idx3 == idx1 ) idx3=rand()%nPop;
      while(idx3 == idx2 || idx3 == idx1 || idx3 == idx_target) idx3=rand()%nPop;
      // produce a child solution
      CIndividual &child = pool[child_idx[i]];
      diff_evo_xoverA(pool[parent_idx[idx_target]], pool[parent_idx[idx1]], pool[parent_idx[idx2]], pool[parent_idx[idx3]], child, CR, F);
      // apply polynomial mutation
      realmutation(child, 1.0/nvar);
      child.obj_eval();
      // update the reference points and other solutions in the neighborhood or the whole population
      update_reference(child); //O(M)
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
	//	   save_pos(filename1);
		   save_front(filename2);
		}
		bef=nfes;
	        nfes += nOffspring;
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
void CMOEAD::penalization(unordered_set<int> &candidates, unordered_set<int> &penalized, vector<double> &distances, vector<set<pair<double, int> > > &ranking_weights)
{
  if(candidates.empty()) return;
  unordered_set<int> tmp = candidates;
  for(auto c_idx:tmp)
  {
     if( distances[c_idx] <= D ) 
     {
	penalized.insert(c_idx);
	candidates.erase(c_idx);
	for(int w = 0; w < nWeight; w++) ranking_weights[w].erase(make_pair(table_fitness[w][c_idx], c_idx));
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
void CMOEAD::replacement_phase()
{
  unordered_set<int> survivors, candidates, penalized;
  vector<double> distances((int)pool.size(), DBL_MAX);
  vector<set<pair<double, int> > > ranking_weights(nWeight, set<pair<double, int>>());
  for(int idx = 0; idx < pool.size(); idx++) candidates.insert(idx); 
  table_fitness_information();
  for(int w = 0; w < nWeight; w++)
  {
     for(auto c_idx:candidates)
	{
	   ranking_weights[w].insert(make_pair(table_fitness[w][c_idx], c_idx));
	}
  }

  pick_best_R2(survivors, candidates, penalized, distances, ranking_weights);
  diversity_information(survivors, candidates, distances);

  while( survivors.size() < nPop)
  {  
     penalization(candidates, penalized, distances, ranking_weights);
     if( candidates.empty()) break;
     pick_best_R2(survivors, candidates, penalized, distances, ranking_weights);
  }
  while( survivors.size() < nPop )   pick_penalized(penalized, survivors, distances);

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
void CMOEAD::pick_best_R2(unordered_set<int> &survivors, unordered_set<int> &candidates, unordered_set<int> &penalized, vector<double> &distances, vector<set<pair<double, int> > > &ranking_weights)
{
  pair<double, int> best(0, -1); 
  vector<double> contribution(pool.size(), 0);
  while(best.second ==-1)
  {
   for(int w = 0; w < nWeight; w++)   
   {
	if(ranking_weights[w].empty()) continue;
	if(survivors.find(ranking_weights[w].begin()->second) != survivors.end()) continue;
        if(ranking_weights[w].size() >1)
        contribution[ranking_weights[w].begin()->second] += max(next(ranking_weights[w].begin(), 2)->first - ranking_weights[w].begin()->first, 0.0001);
//	if(best.first > ranking_weights[w].begin()->first) best = make_pair(ranking_weights[w].begin()->first, ranking_weights[w].begin()->second);
   }
   for(int i =0; i < contribution.size(); i++) if( best.first <  contribution[i]) best = make_pair(contribution[i], i);
   if( best.second ==-1) //remove first current rank
   {
      for(int w = 0; w < nWeight; w++)  if(!ranking_weights[w].empty())ranking_weights[w].erase(ranking_weights[w].begin());
   }
  }
  survivors.insert(best.second);
  candidates.erase(best.second);
  for(auto idx_c:candidates) distances[idx_c] = min(distances[idx_c], distance_var(pool[idx_c].x_var, pool[best.second].x_var));
  for(auto idx_p:penalized) distances[idx_p] = min(distances[idx_p], distance_var(pool[idx_p].x_var, pool[best.second].x_var));
}
#endif
