 algname=(CPDEA MOEAD R2EMOA VSDMOEA)
 ##Calculating the HV considering plus 10% for each objective
    cont=0
 for alg in ${algname[@]}
 do
 #   alg=($alg"/POF/")
    for i in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF1 UF2 UF3 UF4 UF5 UF6 UF7;
    do
 	   rm ${algname[$cont]}/${i}_2
 	   for file in ../$alg/${i}_2_*
 	   do
              ./igdp --r ../Optimals/${i}_2.txt --d $file  >> ${algname[${cont}]}/${i}_2
 	   done
    done
 #################3obj
    for i in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF8 UF9 UF10;
    #for i in UF8 UF9 UF10;
    do
 	   rm ${algname[$cont]}/${i}_3
 	   for file in ../$alg/${i}_3_*
 	   do
              ./igdp --r ../Optimals/${i}_3.txt --d $file  >> ${algname[${cont}]}/${i}_3
 	   done
    done
    cont=$((cont+1))
 done
