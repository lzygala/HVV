#!/bin/sh -e

export X509_USER_PROXY=/afs/cern.ch/user/l/lzygala/.x509up_u115293
voms-proxy-info

source /cvmfs/sft.cern.ch/lcg/views/LCG_97/x86_64-centos7-gcc8-opt/setup.sh

OUTPUTFILE=$2;
INPUTFILE=$3;
NAME=$1;

root -l -b -q  "GenLevelAnalyzer_TTreeMaker.C++(\"$INPUTFILE\")"

xrdcp -f --nopbar  TTree_tmp.root $OUTPUTFILE 

echo $NAME

rm ./*
