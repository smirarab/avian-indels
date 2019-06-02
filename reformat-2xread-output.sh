#!/bin/sh

grep "sequence_indel_" $1 > $1.indel.list.txt
sed "s/{//g" $1.indel.list.txt > t.0
mv t.0 $1.indel.list.txt
sed "s/sequence_indel_/sequence_indel /g" $1.indel.list.txt > t.0
mv t.0 $1.indel.list.txt
sed "s/-/ /g" $1.indel.list.txt > t.0
mv t.0 $1.indel.list.txt

sed "s/\[01\]/?/g" $1 > t.0
mv t.0 $1
sed "s/\[02\]/?/g" $1 > t.0
mv t.0 $1
sed "s/\[03\]/?/g" $1 > t.0
mv t.0 $1
sed "s/\[12\]/?/g" $1 > t.0
mv t.0 $1
sed "s/\[13\]/?/g" $1 > t.0
mv t.0 $1
sed "s/\[23\]/?/g" $1 > t.0
mv t.0 $1
sed "s/\[012\]/?/g" $1 > t.0
mv t.0 $1
sed "s/\[013\]/?/g" $1 > t.0
mv t.0 $1
sed "s/\[023\]/?/g" $1 > t.0
mv t.0 $1
sed "s/\[123\]/?/g" $1 > t.0
mv t.0 $1
sed "s/\[0123\]/?/g" $1 > t.0
mv t.0 $1

sed "1d" $1 > t.0
mv t.0 $1
sed "1d" $1 > t.0
mv t.0 $1
sed "1d" $1 > t.0
mv t.0 $1
