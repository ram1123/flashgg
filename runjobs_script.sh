#!/bin/sh
# The purpose of this bash script is to run fggrunjobs one file at a time since a direct output on eos is not allowed. 
# Setup cms environment and access grid 
cmsenv
voms-proxy-init --voms cms --valid 168:00
# Params
label=$1 
num_events=$2
max_df=0 # max data paths 
max_bf=-1 # max background paths
data_output=$label
data_output+="_Data"
bkg_output=$label
bkg_output+="_Bkg"
# Make output directories if they don't exist 
mkdir -p $data_output;
mkdir -p $bkg_output;
mkdir -p /eos/cms/store/user/atishelm/$data_output;
mkdir -p /eos/cms/store/user/atishelm/$bkg_output;
# Get all files 
data_direc="/afs/cern.ch/work/a/atishelm/2FebFlashgglxplus7/CMSSW_10_2_1/src/flashgg/Taggers/test/Data_Jsons" # where singular data configs are 
bkg_direc="/afs/cern.ch/work/a/atishelm/2FebFlashgglxplus7/CMSSW_10_2_1/src/flashgg/Taggers/test/Bkg_Jsons" # where singular bkg configs are 
unset data_paths # Make sure array name is free in memory 
declare -a data_paths # unassociative array 
unset bkg_paths # Make sure array name is free in memory 
declare -a bkg_paths # unassociative array 
data_path=()
bkg_paths=() 
#Data 
for path in `find $data_direc -name '*.json'`
do
    # Use if statement below if you want to ignore files that contain a specified string 
    #if [[ $path != *"inLHE"* ]]; then # double brackets allows you to use * outside quotes 
    echo "found path: $path"
    data_paths+=("$path"); 
    #fi 
    done 
# Bkg 
for path in `find $bkg_direc -name '*.json'`
do
    # Use if statement below if you want to ignore files that contain a specified string 
    #if [[ $path != *"inLHE"* ]]; then # double brackets allows you to use * outside quotes 
    echo "found path: $path"
    bkg_paths+=("$path"); 
    #fi 
    done 
# Data 
dp_i=0 # data path i 
for path in "${data_paths[@]}"
do 
    if [ "$dp_i" -eq $max_df ] 
    then
        echo "Reached max desired number of data files: $dp_i"
        break
    fi 
    echo "Running on path number $dp_i"
    echo "Running on path: $path"
    command='fggRunJobs.py --load '
    command+=$path
    command+=' -D -P -n 100 -d ' # Might need to be careful not to have too many output files. EOS has a limit. 
    command+=$data_output
    command+=' -x cmsRun Taggers/test/HHWWggTest_cfg.py maxEvents='
    command+=$num_events
    command+=' -q microcentury --no-use-tarball'
    echo "command: $command"
    eval "$command" 
    echo "Finished job for file: $path"
    mv $data_output/*.root /eos/cms/store/user/atishelm/$data_output/
    echo "Finished moving files for: $path"
    echo " "
    dp_i=$(($dp_i+1))
done 
# Background
bp_i=0 # bkg path i 
for path in "${bkg_paths[@]}"
do 
    if [ "$bp_i" -eq $max_bf ] 
    then
        echo "Reached max desired number of background files: $bp_i"
        break
    fi 
    echo "Running on path number $bp_i"
    echo "Running on path: $path"
    command='fggRunJobs.py --load '
    command+=$path
    command+=' -D -P -n 100 -d ' # Might need to be careful not to have too many output files. EOS has a limit. 
    command+=$bkg_output
    command+=' -x cmsRun Taggers/test/HHWWggTest_cfg.py maxEvents='
    command+=$num_events
    command+=' -q microcentury --no-use-tarball'
    echo "command: $command"
    eval "$command" 
    echo "Finished job for file: $path"
    mv $bkg_output/*.root /eos/cms/store/user/atishelm/$bkg_output/
    echo "Finished moving files for: $path"
    echo " "
    bp_i=$(($bp_i+1))
done 
