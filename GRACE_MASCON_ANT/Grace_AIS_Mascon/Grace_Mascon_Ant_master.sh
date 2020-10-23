#!/bin/bash
#export LD_LIBRARY_PATH=/opt/python3/lib:$LD_LIBRARY_PATH
echo "Copy MODELS content .nc files into Maindir"
cp -r MODELS/*.nc .

echo "Make the IMAGES directory"
mkdir OUT
#mkdir OUT
#mkdir WORK
echo "First I call in the first part of the Grace script."
python3 ./AIS_Mass_Change.py
echo "We copy all the jpg results that are not already in the IMAGES folder to the IMAGES folder."
cp *.png OUT

#echo "Copy the content of the WORK directory to the OUT folder."
#cp -r IMAGES/* OUT

