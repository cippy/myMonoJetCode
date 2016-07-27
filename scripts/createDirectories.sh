#!/bin/bash                                             

#launch as --> source <script_name>.sh                                                                                                                                 

wwwDir="/afs/cern.ch/user/m/mciprian/www/"
basePath=$CMSSW_BASE/src

# this macro should be invoked passing the name of config file, where directories'names are provided
# the name passed can be the absolute path or simply it's name, or any relative path, because we will get the <name>.txt part and build again the absolute path
# we are assuming the file is in myMonoJetCode/config directory, so please keep to this convention  

if [ $# -lt 2 ]; 
then 
    echo "ERROR: too few arguments passed to $0. Launch as --> prompt> source $0 <configFile> <last_folder_name>"
    echo "Exiting"
    return 1
fi

fileToGrep=$1
absPathToConfig=${basePath}/myMonoJetCode/config/
fileToGrep=${fileToGrep##*/}  # remove everything before <name>.txt (remove nothing if no '/' is present)
#fileToGrep="monojet_signalRegion_config.txt"
fileToGrep=${absPathToConfig}${fileToGrep}

#check existance of file
if [ ! -f "$fileToGrep" ]; 
then
    echo "ERROR: file $fileToGrep"
    echo "not found in $absPathToConfig !"
    echo "Exiting"
    return 1
fi

echo "changing directory to --> ${basePath}"
echo "setting environment with cmsenv"
cd $basePath
#cmsenv  # it's an alias, so it might not be already known
eval `scramv1 runtime -sh`
echo "Back to previous directory"
cd -
#echo                                                                                                                                                                  
#echo "changing directory to $MY_PATH/myMonoJetCode/src"                                                                     
#echo
#cd ${basePath}/myMonoJetCode/src

dirTocreate=`cat ${fileToGrep} | grep DIRECTORY_PATH | awk '{print $3}'`  # name is 3rd object in line
#finalDir=`cat ${fileToGrep} | grep DIRECTORY_NAME | awk '{print $3}'`
finalDir=$2 # second argument is passed from main, because there might be some suffixes in the name given in config file
finalDir=${dirTocreate}${finalDir}
mkdir -p ${basePath}${finalDir}

# ex: finalDir=/myMonoJetCode/output/monojet/80X/lumi_7p65fb/SignalRegion
#echo "\$finalDir = $finalDir"
#echo $finalDir
#echo ${finalDir#/myMonoJetCode/output/}  # remove shortest pattern matching '/myMonoJetCode/output/' from the beginning
dirToCreateInWWW=${finalDir#/myMonoJetCode/output/}
echo "mkdir -p ${wwwDir}${dirToCreateInWWW}"
mkdir -p ${wwwDir}${dirToCreateInWWW}
cd ${wwwDir}
./copyphp.sh
cd -

