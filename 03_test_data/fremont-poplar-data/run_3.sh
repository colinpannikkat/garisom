# SA of rootBeta independently

echo "EXPERIMENT 3 POPULATION $1" >> exp-3.out
nohup python sa.py -i exp-3_sa-problem.json -o ./exp_3 -m ../../02_program_code -w 30 -s 64 -p $1 2>&1 >> exp-3.out &