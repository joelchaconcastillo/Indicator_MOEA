#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include <queue>
#include <iomanip>
#include "global.h"
#include "recomb.h"
#include "common.h"
#include "individual.h"

class CMOEAD
{

public:
	CMOEAD();
	virtual ~CMOEAD();

	void init_population();                  // initialize the population
	void init_neighbourhood();
	void update_reference(CIndividual &ind);                 // update ideal point which is used in Tchebycheff or NBI method
	double get_distance_near_point( vector<int> &SetA, int index,  vector<CIndividual> &candidates);
	void mate_selection(vector<int> &list, int cid, int size, int type);
	void mate_selection_decision(vector<int> &list, int cid, int type);
	void replacement_phase();

	void evol_population();                                      // DE-based recombination
	

	// execute MOEAD
	void exec_emo(int run);

	void save_front(char savefilename[4024]);       // save the pareto front into files
	void save_pos(char savefilename[4024]);


	void update_parameterD();

	double distance( vector<double> &a, vector<double> &b);
	double distance_obj( vector<double> &a, vector<double> &b);
	vector <CSubproblem> population;
	vector<CIndividual> child_pop;//, best;	// memory solutions

	void operator=(const CMOEAD &moea);

public:

	// algorithm parameters
	long long max_gen, curren_gen;       //  the maximal number of generations and current gen
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
//	D = max(0.0, D);
}
double CMOEAD::distance_obj( vector<double> &a, vector<double> &b)
{
	double dist = 0 ;
   for(int i = 0; i < a.size(); i++)
	{
	   double factor = (a[i]-b[i]);
	   dist += factor*factor;
	}
   return sqrt(dist);
}
double CMOEAD::distance( vector<double> &a, vector<double> &b)
{
	double TElapsed = nfes;
        double TEnd = max_nfes;

	double dist = 0 ;
        double ff = max(1.0 -  (TElapsed / (TEnd*Df)),-0.01);
   for(int i = 0; i < a.size(); i++)
	{
	   double factor = (a[i]-b[i])/(vuppBound[i]-vlowBound[i]);
	   //if(fabs(factor)< 0.0001*ff) return 0.0;
	   dist += factor*factor;
	}
   return sqrt(dist);
}
double CMOEAD::get_distance_near_point( vector<int> &SetA, int index,  vector<CIndividual> &candidates)
{

   double min_distance = INFINITY;

   if(SetA.empty()) return min_distance;
   for(int i = 0 ; i < SetA.size(); i++)
	min_distance = min(min_distance, distance( candidates[SetA[i]].x_var, candidates[index].x_var) );

   return min_distance;
}

void CMOEAD::replacement_phase()
{
       vector<int> selected_pop;
       vector<CIndividual> Candidates;
   
        priority_queue< pair<double, pair<int, int> > > pq;
        double f1;
        for(int i = 0 ; i < pops; i++)
	{
		Candidates.push_back(population[i].indiv);
		Candidates.push_back(child_pop[i]);
	}

        for(int i = 0 ; i < Candidates.size(); i++)
	{
		for(int k = 0; k < pops; k++)
		{
			f1 = fitnessfunction( Candidates[i].y_obj , population[k].namda);
			pq.push(make_pair(-f1, make_pair(i, k)));
		}
	}

        vector<int> penalized;//( Candidates.size(), 1);
	vector<bool> active_subproblem(pops, true), idxpenalized(Candidates.size(), true);

	while(!pq.empty())
	{
	   pair<double, pair<int, int> > data = pq.top();
	   int idxindividual = data.second.first;
	   int idxsubproblem = data.second.second;
	   pq.pop(); 

	   if(!active_subproblem[idxsubproblem]) continue;
	   if( !idxpenalized[idxindividual]) continue;
	   
	   double dist_near = INFINITY; 
	   for(int i = 0 ; i < selected_pop.size(); i++)
	   {
	   	dist_near = min(dist_near, distance(Candidates[selected_pop[i]].x_var, Candidates[idxindividual].x_var));
	   }
	   if( dist_near < D)
	   {
		  penalized.push_back(idxindividual);
		  idxpenalized[idxindividual] = false;
	   }
	   else
	   {
	        selected_pop.push_back(idxindividual);
	 	if( D > 0)
		idxpenalized[idxindividual] = false;
		population[idxsubproblem].indiv = Candidates[idxindividual];
		active_subproblem[idxsubproblem] = false;
	   }
	}	
     vector<double> v_distances(penalized.size(), INFINITY);
     for(int i = 0 ;  i < penalized.size(); i++)
        {
           for(int j = 0; j< selected_pop.size(); j++)
           {
              v_distances[i] = min(v_distances[i], distance( Candidates[penalized[i]].x_var, Candidates[selected_pop[j]].x_var));
           }
        }
       vector<int> unset_subproblem;
        for(int i = 0; i < active_subproblem.size(); i++)
	{
	   if(active_subproblem[i]) unset_subproblem.push_back(i);
	}
	//while(selected_pop.size() < pops)
	while(!unset_subproblem.empty())
	{
	    double maxd = -INFINITY;
            int idx_individual = -1;

            for(int i = 0 ; i < penalized.size(); i++)
            {
                    if( v_distances[i] > maxd)
                    {
                            maxd = v_distances[i];
                            idx_individual = i;
                    }
            }
	    int idx_subproblem = -1;
	    double minfit = INFINITY;
	    //find the best subproblem...
	    for(int i = 0; i < unset_subproblem.size(); i++)
	    {
		double fitness = fitnessfunction( Candidates[penalized[idx_individual]].y_obj , population[unset_subproblem[i]].namda);
		if(fitness < minfit )
		{
		   minfit = fitness;
		   idx_subproblem = i;
		}
	    }

	  for(int i = 0 ; i < penalized.size(); i++)
          {
	     if( i==idx_individual) continue;
             v_distances[i] = min(v_distances[i] , distance(Candidates[penalized[idx_individual]].x_var, Candidates[penalized[i]].x_var ));
          }

	   population[unset_subproblem[idx_subproblem]].indiv = Candidates[penalized[idx_individual]];

//	   selected_pop.push_back(penalized[idx_individual]);
	   iter_swap(penalized.begin() + idx_individual, penalized.end()-1);
	   penalized.pop_back();
	   iter_swap(v_distances.begin() + idx_individual, v_distances.end()-1);
	   v_distances.pop_back();

	   iter_swap(unset_subproblem.begin() + idx_subproblem, unset_subproblem.end()-1);
	   unset_subproblem.pop_back();

	}
}
void CMOEAD::init_population()
{

	idealpoint = vector<double>(nobj, 1.0e+30);

	char filename[1024];
	// Read weight vectors from a data file
//	sprintf(filename,"/home/joel.chacon/Current/MyResearchTopics/MOEA-D-Diversity/MOEAD-DE/vsd-moead-opt/ParameterSetting/Weight/W%dD_%d.dat", nobj, pops);
	sprintf(filename,"%s/ParameterSetting/Weight/W%dD_%d.dat", strpath, nobj, pops);
	std::ifstream readf(filename);

    for(int i=0; i<pops; i++)
	{

		CSubproblem sub;

		// Randomize and evaluate solution
		sub.indiv.rnd_init();
		sub.indiv.obj_eval();

		sub.saved = sub.indiv;

		// Initialize the reference point
		update_reference(sub.indiv);

		// Load weight vectors
		for(int j=0; j<nobj; j++)
		{
		    readf>>sub.namda[j];
		}
	//	printf("\n"); getchar();

		// Save in the population
		population.push_back(sub);
		child_pop.push_back(sub.indiv);
		nfes++;
	}
	readf.close( );
}

void CMOEAD::operator=(const CMOEAD &alg)
{
	//population = alg.population;
}
void CMOEAD::init_neighbourhood()
{
    vector<double> dist   = vector<double>(pops, 0);
        vector<int>    indx   = vector<int>(pops, 0);

        for(int i=0; i<pops; i++)
        {
                // calculate the distances based on weight vectors
                for(int j=0; j<pops; j++)
                {
                    dist[j]    = dist_vector(population[i].namda,population[j].namda);
                        indx[j]  = j;
                }

                // find 'niche' nearest neighboring subproblems
                minfastsort(dist,indx,population.size(),niche);

                // save the indexes of the nearest 'niche' neighboring weight vectors
                for(int k=0; k<niche; k++)
                {
                        population[i].table.push_back(indx[k]);
                }

        }
    dist.clear();
        indx.clear();
}
void CMOEAD::update_reference(CIndividual &ind)
{
	//ind: child solution
	for(int n=0; n<nobj; n++)
	{
		if(ind.y_obj[n]<idealpoint[n])
		{
			idealpoint[n] = ind.y_obj[n];
		}
	}
}
void CMOEAD::mate_selection(vector<int> &list, int cid, int size, int type){
        // list : the set of the indexes of selected mating parents
        // cid  : the id of current subproblem
        // size : the number of selected mating parents
        // type : 1 - neighborhood; otherwise - whole population
        int ss   = population[cid].table.size(), id, parent;
    while(list.size()<size)
        {
                if(type==1){
                    id      = int(ss*rnd_uni(&rnd_uni_init));
                        parent  = population[cid].table[id];
                }
                else
                        parent  = int(population.size()*rnd_uni(&rnd_uni_init));

                // avoid the repeated selection
                bool flag = true;
                for(int i=0; i<list.size(); i++)
                {
                        if(list[i]==parent) // parent is in the list
                        {
                                flag = false;
                                break;
                        }
                }

                if(flag) list.push_back(parent);
        }
}
void CMOEAD::mate_selection_decision(vector<int> &list, int cid, int type){
        // list : the set of the indexes of selected mating parents
        // cid  : the id of current subproblem
        // size : the number of selected mating parents
        // type : 1 - neighborhood; otherwise - whole population
	if(type==1)
	{
        	vector<pair<double, int>> dists;
                for(int i = 0; i < population.size(); i++)
        	{
        		if(i == cid) continue;
        		dists.push_back( make_pair(distance_obj(population[cid].indiv.y_obj,population[i].indiv.y_obj), i));
        	}
        	sort(dists.begin(), dists.end());
		random_shuffle(dists.begin(), dists.begin()+niche+1);//random_shuffle is exclusive
		list.push_back(dists[0].second);
		list.push_back(dists[1].second);
		list.push_back(dists[2].second);
	}
	else
	{
		int idx1 = int(population.size()*rnd_uni(&rnd_uni_init));
	        while( idx1 == cid)
		{
		   idx1 = int(population.size()*rnd_uni(&rnd_uni_init));
		}
	      int idx2 = int(population.size()*rnd_uni(&rnd_uni_init));
	        while( idx2 == cid || idx2 == idx1)
		{
		   idx2 = int(population.size()*rnd_uni(&rnd_uni_init));
		}
	      int idx3 = int(population.size()*rnd_uni(&rnd_uni_init));
		while( idx3 == cid || idx3 == idx1 ||  idx3==idx2)
		{
		   idx3 = int(population.size()*rnd_uni(&rnd_uni_init));
		}
		list.push_back(idx1);
		list.push_back(idx2);
		list.push_back(idx3);
	}
}
void CMOEAD::evol_population()
{

    for(int sub=0; sub<pops; sub++)
	{

                // select the indexes of mating parents
                vector<int> plist;
                mate_selection(plist, sub, 3, 2);  // neighborhood selection

		// produce a child solution
		CIndividual child;

		diff_evo_xoverA(population[sub].indiv,population[plist[0]].indiv,population[plist[1]].indiv, population[plist[2]].indiv, child, CR, F);

		// apply polynomial mutation
		realmutation(child, 1.0/nvar);
		child.obj_eval();
		// update the reference points and other solutions in the neighborhood or the whole population
		update_reference(child);

		child_pop[sub] = child;
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

	// initialization
	nfes      = 0;
	init_population();
	init_neighbourhood();

	sprintf(filename1,"%s/POS/v2_POS_MOEAD_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_CR_%lf_F_%lf",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar), Df, CR, F);
	sprintf(filename2,"%s/POF/v2_POF_MOEAD_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_CR_%lf_F_%lf",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar), Df, CR, F);
        long long current = nfes;
	long long accumulator = 0, bef = nfes;
	save_pos(filename1);
        save_front(filename2);

	while(nfes<max_nfes)
	{
		update_parameterD();
		evol_population();
		accumulator += nfes - bef ;
        //        if(accumulator > 0.1*(max_nfes)  )
	//	{
	//           accumulator -= 0.1*(max_nfes);
	//	   save_pos(filename1);
	//	   save_front(filename2);
	//	}
		bef=nfes;
	        nfes += pops;
	}
		save_pos(filename1);
		save_front(filename2);
	population.clear();
	idealpoint.clear();

}


void CMOEAD::save_front(char saveFilename[4024])
{

    std::fstream fout;
	//fout.open(saveFilename,std::ios::out);
	fout.open(saveFilename,fstream::app|fstream::out );
	for(int n=0; n<pops; n++)
	{
		//for(int k=0;k<nobj;k++)
		//	fout<<best[n].y_obj[k]<<"  ";
		for(int k=0;k<nobj;k++)
			fout<<population[n].indiv.y_obj[k]<<"  ";
	for(int k=0;k<nobj;k++)
			fout<<child_pop[n].y_obj[k]<<"  ";

	fout<< "|||| " << fitnessfunction( population[n].indiv.y_obj , population[n].namda) << " ";
	fout << "|||| " << fitnessfunction( child_pop[n].y_obj , population[n].namda) << " ";
		fout<<"\n";
	}
	fout.close();
}

void CMOEAD::save_pos(char saveFilename[4024])
{
    std::fstream fout;
	//fout.open(saveFilename,std::ios::out);
	fout.open(saveFilename, fstream::app|fstream::out);
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nvar;k++)
			fout<<population[n].indiv.x_var[k] << "  ";
		for(int k=0;k<nvar;k++)
			fout<<population[n].indiv.x_var[k]/(vuppBound[k]-vlowBound[k]) << "  ";

			//fout<<population[n].indiv.x_var[k]<< fixed << setprecision(30) << "  ";
//	  for(int k=0;k<nvar;k++)
//			fout<<best[n].x_var[k]<<"  ";
//	  for(int k=0;k<nvar;k++)
//			fout<<child_pop[n].x_var[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}



#endif
