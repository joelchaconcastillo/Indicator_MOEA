PATH1=/home/joel.chacon/Final_Experiment_Dominance/StoppingCriterion4_3/MOEA_D/maxnfes_2.5pow6/POF
for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF1 UF2 UF3 UF4 UF5 UF6 UF7
do
   for sed in {1..35}   
   do
   tail -100 $PATH1/POF_MOEAD_${instance}_RUN${sed}_seed*_2_nvar_* > MOEAD/${instance}_2_${sed}
   done
done
for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF8 UF9 UF10
do
   for sed in {1..35}   
   do
   tail -100 $PATH1/POF_MOEAD_${instance}_RUN${sed}_seed*_3_nvar_* > MOEAD/${instance}_3_${sed}
   done
done
PATH1=/home/joel.chacon/Final_Experiment_Dominance/StoppingCriterion4_3/R2-EMOA/maxnfes_2.5pow6/POF
for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF1 UF2 UF3 UF4 UF5 UF6 UF7
do
   for sed in {1..35}
   do
   cat $PATH1/${instance}_2_${sed}_nvar_*[!R] | tail -100 > R2EMOA/${instance}_2_${sed}
   done
done
for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF8 UF9 UF10
do
   for sed in {1..35}
   do
   cat $PATH1/${instance}_3_${sed}_nvar_*[!R] | tail -100 > R2EMOA/${instance}_3_${sed}
   done
done
PATH1=/home/joel.chacon/Final_Experiment_Dominance/StoppingCriterion4_3/VSD_MOEA/released-vsd-moea_maxnfes_2.5pow6/POF
DF=0.5
T=LINEAL
for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF1 UF2 UF3 UF4 UF5 UF6 UF7
do
   for sed in {1..35}
   do
   tail -100 $PATH1/POF_VSD-MOEA_${instance}_RUN_${sed}_*_nvar_*_nobj_2_Di_0.400000_nfes_250*_Df_${DF}*pc_${pc}*_$T > VSDMOEA/${instance}_2_${sed}
   done
done
for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF8 UF9 UF10
do
   for sed in {1..35}
   do
   tail -100 $PATH1/POF_VSD-MOEA_${instance}_RUN_${sed}_*_nvar_*_nobj_3_Di_0.400000_nfes_250*_Df_${DF}*pc_${pc}*_$T > VSDMOEA/${instance}_3_${sed}
   done
done
maxfes=2500000
T=CPDEA
PATH1=/home/joel.chacon/Final_Experiment_Dominance/StoppingCriterion4_3/PlatEMO_xover_probability/PlatEMO/Data/CPDEA
for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF1 UF2 UF3 UF4 UF5 UF6 UF7
do
   for sed in {1..35}
   do
     tail -100 $PATH1/OBJ_*_${instance}_M2_D*_maxFE_${maxfes}_pc_1.0*_${sed}.txt > ${T}/${instance}_2_${sed}
   done
done
for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF8 UF9 UF10
do
   for sed in {1..35}
   do
     tail -100 $PATH1/OBJ_*_${instance}_M3_D*_maxFE_${maxfes}_pc_0.6_${sed}.txt > ${T}/${instance}_3_${sed}
   done
done
