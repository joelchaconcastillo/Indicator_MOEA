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
        void update_fronts(); 
	void dominance_information_new(int idx_new);
	void dominance_information_remove(int new_child_idx);
   private:
	vector <CIndividual> pool;
	vector<int> child_idx, parent_idx, inv_parent_idx;
	vector<vector<double> > namda;     // weight vector
        vector<int> Np,Rp;//rank
	vector<unordered_set<int> > Sp;//dominated indexes and inverse
        vector<vector<int> > fronts ;

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
	double TElapsed = nfes;
        double TEnd = max_nfes;
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
void CMOEAD::replacement_phase()
{
   vector<int> candidates = parent_idx;
   vector< set<>> 
   parent_idx.clear();
   candidates.insert(candidates.end(), child_idx.begin(), child_idx.end());
   
   pair<double, int> idx_R2(DBL_MIN, -1);
   for(auto idx:candidates)
   {
      
      idx_R2 = max(idx_R2, fitnessfunction(pool[idx]).y_obj,  )	
   }

   
}
void CMOEAD::init_population()
{
    idealpoint = vector<double>(nobj, 1.0e+30);
    char filename[1024];
    // Read weight vectors from a data file
    sprintf(filename,"%s/ParameterSetting/Weight/W%dD_%d.dat", strpath, nobj, nWeight);
    std::ifstream readf(filename);
    namda.resize(nWeight, vector<double> (nobj, 0.0));
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
void CMOEAD::dominance_information_remove(int ridx)
{
  for(auto idx:Sp[ridx]) Np[idx]--;
  for(auto idx:parent_idx) Sp[idx].erase(ridx);//it can be faster...
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
   vector<int> current_Np = Np;
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
#endif
