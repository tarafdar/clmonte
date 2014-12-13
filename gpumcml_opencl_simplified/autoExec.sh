#!/bin/bash

echo -e "A_rz\tRd_ra\tTt_ra" > result.txt

#For 1M photons
echo "1M" >> result.txt
mkdir Output1M
for i in {1..5}
do
  cmd="./gpumcml_opencl test.mci Output1M/output$i.txt 1000000"
  eval $cmd
done

for f in Output1M/output*
do
  lineresult=""
  while read line
  do
    lineresult="$lineresult$line\t"
  done < $f
  echo -e $lineresult >> result.txt
done
echo -e "\n\n" >> result.txt

#For 2M photons
echo "2M" >> result.txt
mkdir Output2M
for i in {1..5}
do
  cmd="./gpumcml_opencl test.mci Output2M/output$i.txt 2000000"
  eval $cmd
done

for f in Output2M/output*
do
  lineresult=""
  while read line
  do
    lineresult="$lineresult$line\t"
  done < $f
  echo -e $lineresult >> result.txt
done
echo -e "\n\n" >> result.txt

#For 4M photons
echo "4M" >> result.txt
mkdir Output4M
for i in {1..5}
do
  cmd="./gpumcml_opencl test.mci Output4M/output$i.txt 4000000"
  eval $cmd
done

for f in Output4M/output*
do
  lineresult=""
  while read line
  do
    lineresult="$lineresult$line\t"
  done < $f
  echo -e $lineresult >> result.txt
done
echo -e "\n\n" >> result.txt

#For 8M photons
echo "8M" >> result.txt
mkdir Output8M
for i in {1..5}
do
  cmd="./gpumcml_opencl test.mci Output8M/output$i.txt 8000000"
  eval $cmd
done

for f in Output8M/output*
do
  lineresult=""
  while read line
  do
    lineresult="$lineresult$line\t"
  done < $f
  echo -e $lineresult >> result.txt
done
echo -e "\n\n" >> result.txt

#For 16M photons
echo "16M" >> result.txt
mkdir Output16M
for i in {1..5}
do
  cmd="./gpumcml_opencl test.mci Output16M/output$i.txt 16000000"
  eval $cmd
done

for f in Output16M/output*
do
  lineresult=""
  while read line
  do
    lineresult="$lineresult$line\t"
  done < $f
  echo -e $lineresult >> result.txt
done
echo -e "\n\n" >> result.txt

#clear the raw data, comment out this line if want to keep them
rm -rf Output*M
