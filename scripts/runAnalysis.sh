#!/bin/bash
# ok, it's working
#launch as --> source <script_name>.sh
MY_PATH=$CMSSW_BASE/src
echo "changing directory to --> $MY_PATH and setting environment with cmsenv"
cd $MY_PATH
cmsenv
#echo
#echo "changing directory to $MY_PATH/myMonoJetCode/src"
echo
MY_PATH=$MY_PATH/myMonoJetCode/src
cd $MY_PATH
echo "Compiling code with make command"
echo
make
echo
cd ..
MACRO_FOLDER="macro"
MACRO_NAME="distribution.C"
EXE_NAME="distribution"
echo "Compiling $MACRO_FOLDER/$MACRO_NAME"
g++ -Wall -pedantic -lm -o $MACRO_FOLDER/$EXE_NAME $MACRO_FOLDER/$MACRO_NAME `rootlib`
echo
echo "Calling --> bsub -q cmsan -o std_output.txt -J runAnalysis \"source $MY_PATH/scripts/AnaLaunch.sh\""
#bsub -q cmsan -Is /bin/bash < AnaLaunch.sh
bsub -q cmsan -o std_output.txt -J runAnalysis "source $MY_PATH/scripts/AnaLaunch.sh"
echo
echo "THE END!"

