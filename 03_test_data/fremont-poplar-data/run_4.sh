# SA of baperga and height

echo "EXPERIMENT 4 POPULATION $1" >> exp-4.out
nohup python sa.py -i exp-4_sa-problem.json -o ./exp_4 -m ../../02_program_code -w 30 -s 64 -p $1 2>&1 >> exp-4.out &