# SA of leafAreaIndex independently

echo "EXPERIMENT 2 POPULATION $1" >> exp-2.out
nohup python sa.py -i exp-2_sa-problem.json -o ./exp_2 -m ../../02_program_code -w 30 -s 64 -p $1 2>&1 >> exp-2.out &