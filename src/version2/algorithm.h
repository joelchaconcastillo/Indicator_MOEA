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
	void select_extremes(unordered_set<int> &candidates, unordered_set<int> &survivors, unordered_set<int> &survivors_front);
	void penalization(unordered_set<int> &candidates, unordered_set<int> &penalized, vector<double> &distances, unordered_set<int> &candidates_front);
        void dominance_information(); 
        void update_fronts(); 
	void pick_penalized(unordered_set<int> &penalized, unordered_set<int> &survivors, vector<double> &distances);
	void R2_contribution_subset(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors, unordered_set<int> &survivors_front, unordered_set<int> &penalized, vector<set<pair<double, int> > > &survivors_weight, vector<double> &distances);
	void dominance_information_new(int idx_new);
	void dominance_information_remove(int rm_idx, unordered_set<int> &current_front);
	void lowest_front(unordered_set<int> &candidates, unordered_set<int> &candidates_front);
	void update_lowest_front(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors_front);
	void diversity_information(unordered_set<int> &survivors, unordered_set<int> &penalized, vector<double> &distances);

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
void CMOEAD::replacement_phase()
{
  unordered_set<int> survivors, candidates, penalized, survivors_front, candidates_front;
  vector<double> distances((int)pool.size(), DBL_MAX);
  vector<set<pair<double, int> > > survivors_weight;
  for(int idx = 0; idx < pool.size(); idx++) candidates.insert(idx); 
  table_fitness_information();

  dominance_information(); 
  lowest_front(candidates, candidates_front);

  R2_contribution_subset(candidates, candidates_front, survivors, survivors_front, penalized, survivors_weight, distances);
  diversity_information(survivors, candidates, distances);
//  cout<<candidates.size()<<" "<<candidates_front.size()<<endl;
//  cout<<survivors.size()<< " "<<survivors_front.size()<<endl;
  while( survivors.size() < nPop)
  {  
     penalization(candidates, penalized, distances, candidates_front);
     if( candidates.empty())
     {
       //get the farthest penalized individual and move to candidates
       pick_penalized(penalized, survivors, distances);
       continue;
     }
     if( candidates_front.empty())
     update_lowest_front(candidates, candidates_front, survivors_front);
     //get list of R2 contribution
     R2_contribution_subset(candidates, candidates_front, survivors, survivors_front, penalized, survivors_weight, distances);
  }
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
   for(int i = 0; i < nOffspring; i++)
   {
      int idx_target = rand()%nPop;
      int idx1=rand()% nPop, idx2=rand()%nPop, idx3=rand()%nPop;
      while(idx1 == idx_target) idx1=rand()%nPop;
      while(idx2 == idx1 || idx2 == idx_target) idx2=rand()%nPop;
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
                if(accumulator > 0.01*(max_nfes)  )
		{
	           accumulator -= 0.01*(max_nfes);
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

/*
   Dominance information....
*/
void CMOEAD::dominance_information_new(int idx_new)
{
   Sp[idx_new].clear();
   Np[idx_new] = 0;
   for(auto idx_p:parent_idx)
   {
      if( pool[idx_new] < pool[idx_p] )
      {
         Np[idx_p]++;
         Sp[idx_new].insert(idx_p);
      }
      else if( pool[idx_p] < pool[idx_new])
      {
        Np[idx_new]++;
        Sp[idx_p].insert(idx_new);
      }
   }
}
void CMOEAD::dominance_information_remove(int rm_idx, unordered_set<int> &current_front)
{
  for(auto idx:Sp[rm_idx])
  {
     Np[idx]--;
    if(Np[idx] == 0) current_front.insert(idx);
  }
  for(auto idx:parent_idx) Sp[idx].erase(rm_idx);
}
void CMOEAD::update_fronts()
{
   vector<int> current_Np = Np;
   fronts.assign(1, vector<int>());
   int rank = 0;
   for(auto idx:parent_idx)
   {
      if(current_Np[idx]==0)
      {
        fronts[rank].push_back(idx);
        Rp[idx] = rank;
      }
   }
   while(true)
   {
      vector<int> next_front;
      for(auto idx:fronts[rank])
      {
	for(auto idx_dominated:Sp[idx])
        {
	  current_Np[idx_dominated]--;
          if(current_Np[idx_dominated]  == 0) 
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
   for(auto pidx1:parent_idx)
   {
      for(auto pidx2:parent_idx)
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
//   vector<int> current_Np = Np;
//   while(true)
//   {
//      vector<int> next_front;
//      for(auto idx:fronts[rank])
//      {
//	for(auto idx_dominated:Sp[idx])
//        {
//	  current_Np[idx_dominated]--;
//          if(current_Np[idx_dominated]  == 0) 
//          {
//	     next_front.push_back(idx_dominated);
//	     Rp[idx_dominated] = rank+1;
//          }
//        }
//      }
//      if(next_front.empty()) break;
//      fronts.push_back(next_front);
//      rank++;
//   }
}
void CMOEAD::select_extremes(unordered_set<int> &candidates, unordered_set<int> &survivors, unordered_set<int> &survivors_front)
{
 for(int m = 0; m < nobj; m++)
  {
     pair<double, unordered_set<int>::iterator> b_best =  make_pair(DBL_MAX, candidates.begin());
     for(auto idx = candidates.begin(); idx != candidates.end(); idx++)
        if( pool[*idx].y_obj[m] < b_best.first ) b_best =  make_pair(pool[*idx].y_obj[m], idx);

     if( survivors.find(*b_best.second)  == survivors.end() )
     {
        survivors.insert(*b_best.second);
        survivors_front.insert(*b_best.second);
	candidates.erase(b_best.second);
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
void CMOEAD::penalization(unordered_set<int> &candidates, unordered_set<int> &penalized, vector<double> &distances, unordered_set<int> &candidates_front)
{
  if(candidates.empty()) return;
  unordered_set<int> tmp = candidates;
  for(auto c_idx = tmp.begin(); c_idx != tmp.end(); c_idx++)
  {
     if( distances[*c_idx] < D) 
     {
	penalized.insert(*c_idx);
	dominance_information_remove(*c_idx, candidates_front);
        candidates_front.erase(*c_idx);
	candidates.erase(*c_idx);
     }
  } 
}
void CMOEAD::pick_penalized(unordered_set<int> &penalized, unordered_set<int> &survivors, vector<double> &distances)
{
   pair<double, int> max_spread = make_pair(-100000, -1);
   for( auto p_idx : penalized)
   {
      if( distances[p_idx] > max_spread.first) max_spread = make_pair(distances[p_idx], p_idx);
   }
   survivors.insert(max_spread.second);
   for( auto p_idx : penalized)
   {
	if(p_idx == max_spread.second) continue;
      distances[p_idx] = min(distances[p_idx], distance_var( pool[p_idx].x_var, pool[max_spread.second].x_var));
   }
   penalized.erase(max_spread.second);
}
void CMOEAD::lowest_front(unordered_set<int> &candidates, unordered_set<int> &candidates_front)
{
   for(auto c_idx:candidates) if( Np[c_idx] == 0) candidates_front.insert(c_idx);
}
void CMOEAD::update_lowest_front(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors_front)
{
  for(auto s_idx:survivors_front)
  {
     if( Np[s_idx] == 0 ) 
     {
	for(auto idx:Sp[s_idx])
        {
           Np[idx]--;
           if(Np[idx] == 0) candidates_front.insert(idx);
        }
        Np[s_idx]--;
     }
  }
  survivors_front.clear();
}
void CMOEAD::R2_contribution_subset(unordered_set<int> &candidates, unordered_set<int> &candidates_front, unordered_set<int> &survivors, unordered_set<int> &survivors_front, unordered_set<int> &penalized, vector<set<pair<double, int> > > &survivors_weight, vector<double> &distances)
{
   //select_extremes(candidates, survivors, survivors_front);
  if( survivors_front.empty())
  {
     survivors_weight.assign((int)pool.size(), set<pair<double, int>>());
     pair<double, int> max_contribution(DBL_MAX, -1);
     for(auto c_idx:candidates_front)
     {
	double sum = 0.0;
        for(int w_id=0; w_id < nWeight; w_id++)
        {
           sum += table_fitness[w_id][c_idx];
        }
	if( max_contribution.first > sum) max_contribution = make_pair(sum, c_idx);
     }
     survivors_front.insert(max_contribution.second);
     survivors.insert(max_contribution.second);
     candidates.erase(max_contribution.second);
     candidates_front.erase(max_contribution.second);
     for(int w_id=0; w_id < nWeight; w_id++)  survivors_weight[w_id].insert(make_pair(table_fitness[w_id][max_contribution.second], max_contribution.second));
     return;
  }
  pair<double, int> max_contribution(-10000, -1);
  for(auto c_f:candidates_front)
  {
    double Total_contribution = 0.0;
    for(int w_id = 0; w_id < nWeight; w_id++)
       Total_contribution += max(0.0, survivors_weight[w_id].begin()->first - table_fitness[w_id][c_f]);
    if(Total_contribution > max_contribution.first) max_contribution = make_pair(Total_contribution, c_f);
  }
  survivors_front.insert(max_contribution.second);
  survivors.insert(max_contribution.second);
  candidates.erase(max_contribution.second);
  candidates_front.erase(max_contribution.second);
  for(int w_id=0; w_id < nWeight; w_id++)  survivors_weight[w_id].insert(make_pair(table_fitness[w_id][max_contribution.second], max_contribution.second));
  for(auto idx_c:candidates) distances[idx_c] = min(distances[idx_c], distance_var(pool[idx_c].x_var, pool[max_contribution.second].x_var));
  for(auto idx_p:penalized) distances[idx_p] = min(distances[idx_p], distance_var(pool[idx_p].x_var, pool[max_contribution.second].x_var));

}
#endif
