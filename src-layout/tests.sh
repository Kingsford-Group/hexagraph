#!/bin/sh

graphs=(circuit7,call3,flor15,tree50,dolphin,flowchart,karate)

for g in $graphs ; do
    for i in 1 2 3 4 5 6 7 8 9 10 ; do
        ./tiling_layout.py "$g-t$i" ../paper-gd/$g-in.graph > "$g-ti.out"
    done
done
