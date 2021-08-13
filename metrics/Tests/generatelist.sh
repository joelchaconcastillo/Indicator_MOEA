#PATHA=/home/joel.chacon/Current/MyResearchTopics/Data_Dominance/StateOfTheArt_2/HV
PATHA=../HV/  #/home/joel.chacon/Current/MyResearchTopics/Data_Decomposition/StateOfTheArt/HV
 for i in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF1 UF2 UF3 UF4 UF5 UF6 UF7;
do
    echo "--"${i}_2; 
      for j in VSDMOEA CPDEA MOEAD R2EMOA 
      do
        echo ${PATHA}/${j}/${i}_2
      done
done

for i in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF8 UF9 UF10;
do
   echo "--"${i}_3; 
    for j in VSDMOEA CPDEA MOEAD R2EMOA 
    do
      echo ${PATHA}/${j}/${i}_3
    done
done
