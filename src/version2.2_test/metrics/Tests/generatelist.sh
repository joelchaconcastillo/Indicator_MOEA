PATHA=../HV/  #/home/joel.chacon/Current/MyResearchTopics/Data_Decomposition/StateOfTheArt/HV

   for i in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF1 UF2 UF3 UF4 UF5 UF6 UF7;
   do
     echo "--"${i}_2; 
     #for Di in 0.0 0.2 0.4 0.6 0.8 1.0
     for Di in 0.4
     do
       #echo ${PATHA}/POF/${i}_2_${Di}
       echo ${PATHA}/POF/${i}_2
     done
   done
   for i in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF8 UF9 UF10;
   do
     echo "--"${i}_3; 
     #for Di in 0.0 0.2 0.4 0.6 0.8 1.0
     for Di in 0.4
     do
       #echo ${PATHA}/POF/${i}_3_${Di}
       echo ${PATHA}/POF/${i}_3
     done
   done


