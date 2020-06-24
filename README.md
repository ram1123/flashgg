
# Cloning the flashgg Repository

```bash
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_5_0 
cd CMSSW_10_5_0/src
cmsenv
git cms-init
cd $CMSSW_BASE/src 
git clone -b FullLeptonic_V1.1 git@github.com:chuwang1/flashgg.git
source flashgg/setup_flashgg.sh
cd $CMSSW_BASE/src
scram b -j 9
```

# Setting up a voms Proxy

```bash
cmsenv
cd $CMSSW_BASE/src/flashgg
voms-proxy-init --voms cms --valid 168:00
sh proxy.sh X509_USER_PROXY
```

# HHWWggTagProducer

## Running Locally

```bash
sh testSig.sh
```

## Running on Condor

Before you run those commands,you should change the output directory to where you want to Save.
you can change it in [MetaData/python/parallel.py#L253](https://github.com/IHEP-CMS-HHWWgg/flashgg/blob/310ab91cfb51c2e490148e3836cec54c627b5fa2/MetaData/python/parallel.py#L253)
filesOutPath = '/afs/cern.ch/user/c/chuw/work/HHWWgg/FullLep/CMSSW_10_5_0/src/flashgg/Samples/'

```bash
# Signal:

sh HHWWgg_Run_Jobs.sh --labelName HHWWggTaggerTest_Signal --nEvents all --json Taggers/test/HHWWgg_Full/HHWWgg_Signal_2017.json -g -w

#Data:

sh HHWWgg_Run_Jobs.sh --labelName HHWWggTaggerTest_Data --nEvents all --json Taggers/test/HHWWgg_Full/HHWWgg_Data_All_2017.json -g -w
```

# Dump Variables

You can dump the variables you want to dump by edit:Systematics/python/HHWWggCustomize.py line:16
Systematics/test/workspaceStd.py line:289
DataFormats/src/HHWWggTag.cc
DataFormats/interface/HHWWggTag.h


