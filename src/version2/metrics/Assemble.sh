
PATH1=../POF

DI=0.4;
  for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF1 UF2 UF3 UF4 UF5 UF6 UF7
  do
     for sed in {1..35}
     do
     tail -100 $PATH1/POF_*_${instance}_RUN${sed}_seed_*nobj_2_*_DI_${DI}*CR_0.0*_F_*0.75* | cut -f1,3 -d' '  > POF/${instance}_2_${sed}_${DI}
     done
  done
  for instance in WFG1 WFG2 WFG3 WFG4 WFG5 WFG6 WFG7 WFG8 WFG9 DTLZ1 DTLZ2 DTLZ3 DTLZ4 DTLZ5 DTLZ6 DTLZ7 UF8 UF9 UF10
  do
     for sed in {1..35}
     do
     tail -100 $PATH1/POF_*_${instance}_RUN${sed}_seed_*nobj_3*_DI_${DI}*CR_0.0*_F_*0.75* | cut -f1,3,5 -d' '  > POF/${instance}_3_${sed}_${DI}
     done
  done
