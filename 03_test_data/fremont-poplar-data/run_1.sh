# SA of leafAreaIndex and RootBeta

echo "EXPERIMENT 1 POPULATION $1" >> exp-1.out
nohup python sa.py -i exp-1_sa-problem.json -o ./exp_1 -m ../../02_program_code -w 30 -s 64 -p $1 2>&1 >> exp-1.out &