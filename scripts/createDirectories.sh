#!/bin/bash                                             

#launch as --> source <script_name>.sh                                                                                                                                 

wwwDir="/afs/cern.ch/user/m/mciprian/www/"
basePath=$CMSSW_BASE/src

# this macro should be invoked passing the name of config file, where directories'names are provided
# the name passed can be the absolute path or simply it's name, or any relatve path, because we will get the <name>.txt part and build again the absolute path
# we are assuming the file is in myMonoJetCode/config directory, so please keep to this convention  

if [ $# -eq 0 ]; 
then 
    echo "ERROR: no argument passed to $0. Launch as --> prompt> source $0 <configFile>"
    echo "Exiting and returning 5"
    return 5
fi

fileToGrep=$1
absPathToConfig=${basePath}/myMonoJetCode/config/
fileToGrep=${fileToGrep#*/}  # remove everything before <name>.txt (remove nothing if no '/' is present
#fileToGrep="monojet_signalRegion_config.txt"
fileToGrep=${absPathToConfig=}${fileToGrep}

#check existance of file
if [ ! -f "$fileToGrep" ]; then
    echo "ERROR: file $fileToGrep"
    echo "not found in $absPathToConfig !"
    echo "Exiting and returning 5"
    return 5
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
finalDir=`cat ${fileToGrep} | grep DIRECTORY_NAME | awk '{print $3}'`
finalDir=${dirTocreate}${finalDir}
mkdir -p ${basePath}${finalDir}

# ex: finalDir=/myMonoJetCode/output/monojet/80X/lumi_7p65fb/SignalRegion
#echo "\$finalDir = $finalDir"
#echo $finalDir
echo ${finalDir#/myMonoJetCode/output/}  # remove shortestpattern matching '/myMonoJetCode/output/' from the beginning
dirToCreateInWWW=${finalDir#/myMonoJetCode/output/}
echo "mkdir -p ${wwwDir}${dirToCreateInWWW}"
mkdir -p ${wwwDir}${dirToCreateInWWW}