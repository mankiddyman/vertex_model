#!/bin/bash
result=$(ls simulations/*pkl | wc -l)
for ((i=1;i<=result;i++))
do
    echo "running experiment ${i} out of ${result}"
    a=$((i-1))
    python make_df.py $a
done

