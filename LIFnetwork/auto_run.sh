#!/bin/sh

seed_num=$1
echo $seed_num

i=1
while [ $i -le $2 ]
    do 
        python run_w_various_g.py --seed $seed_num
        i=$(($i+1))
        seed_num=$(($seed_num+10))
        # echo $i
        # i=$(($i*10))
done
# for i in 1 2 3
# i=0
# for ((i=0; i<10; i++))
# while [${i} -le 5]
#     do 
#         echo $i
#         i=$(($i+1))
#         # seed_num=$(($sed_num+330))
#         # python run_w_various_g.py --seed $seed_num
# done