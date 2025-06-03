# SA of baperga only

echo "EXPERIMENT 6 POPULATION $1" >> exp-6.out
nohup python sa.py -i exp-6_sa-problem.json -o ./exp_6 -m ../../02_program_code -w 30 -s 64 -p $1 2>&1 >> exp-6.out &