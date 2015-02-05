#!/bin/sh

mkdir outputToy
echo "max events from CRAB: $MaxEvents"
n="$MaxEvents"
./combine 110GeVmodel.root -M MarkovChainMC -D data_mass -m 110 -t $n -s -1 -S 1 --rMin=0 --rMax=60 -b 3000 --iteration 50000 --tries 2 --optimizeSim=1
rm CMS-HGG.root
rm 110GeVmodel.root
mv *.root outputToy/
rm *.root
tar cvfz outputToy.tgz outputToy/
