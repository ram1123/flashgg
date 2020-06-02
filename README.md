
Cloning the flashgg Repository
------

``````
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_5_0 
cd CMSSW_10_5_0/src
cmsenv
git cms-init
cd $CMSSW_BASE/src 
git clone -b HHWWgg_dev https://github.com/atishelmanch/flashgg 
source flashgg/setup_flashgg.sh
cd $CMSSW_BASE/src
scram b -j 9
```````
Setting up a voms Proxy
-------
``````
cmsenv
voms-proxy-init --voms cms --valid 168:00
. proxy.sh X509_USER_PROXY
````````
HHWWggTagProducer
------------------
Running Locally
---------
```````
sh testSig.sh
```````
Running on Condor
--------
Before you run those cmds,you should change the output directory to where you want to Save.
you can change it in MetaData/python/parallel.py line:253
filesOutPath = '/afs/cern.ch/user/c/chuw/work/HHWWgg/FullLep/CMSSW_10_5_0/src/flashgg/Samples/'

------------
`````
Signal:
-------

sh HHWWgg_Run_Jobs.sh --labelName HHWWggTaggerTest_Signal --nEvents all --json Taggers/test/HHWWgg_Full/HHWWgg_Signal_2017.json -g -w

Data:
------
sh HHWWgg_Run_Jobs.sh --labelName HHWWggTaggerTest_Data --nEvents all --json Taggers/test/HHWWgg_Full/HHWWgg_Data_All_2017.json -g -w
`````




