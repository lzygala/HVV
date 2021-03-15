#!/bin/sh -e

export X509_USER_PROXY=/afs/cern.ch/user/l/lzygala/x509up_u115293
voms-proxy-info

xrdcp root://eoscms.cern.ch//store/user/lzygala/HVV/MG5_aMC_v2_9_2.tar.gz ./
tar -xzf MG5_aMC_v2_9_2.tar.gz
#rm MG5_aMC_v2_9_2.tar
cd MG5_aMC_v2_9_2/pp_hvvjj/

echo -e "evaluate"
python2.7 ./bin/generate_events 

JOBID=$1;  


echo -e "Copying result";
xrdcp -f --nopbar  Events/run_01/unweighted_events.lhe.gz root://eoscms.cern.ch//store/user/lzygala/HVV/pp_hvvjj/unweighted_events_${JOBID}.lhe.gz;

echo -e "DONE";

