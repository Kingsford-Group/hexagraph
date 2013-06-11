#!/bin/sh

x=circuit7
graphs="call3 flor15 tree50 dolphin flowchart karate chrom10"

mkdir tests
cd tests
for g in $graphs ; do
    echo "S Graph: $g"
    for i in 1 2 3 4 5 6 7 8 9 10 ; do
        ../tiling_graph.py "$g-t$i" ../../paper-gd/$g-in.graph > "$g-t$i.out"
    done
done

