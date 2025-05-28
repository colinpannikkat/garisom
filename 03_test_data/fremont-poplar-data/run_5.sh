# SA of baperga, height, leafAreaIndex, and RootBeta

echo "EXPERIMENT 5 POPULATION $1" >> exp-5.out
nohup python sa.py -i exp-5_sa-problem.json -o ./exp_5 -m ../../02_program_code -w 30 -s 64 -p $1 2>&1 >> exp-5.out &