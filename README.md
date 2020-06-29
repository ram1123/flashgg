# General Introduction

# Code setup

```bash
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_5_0 
cd CMSSW_10_5_0/src
cmsenv
git cms-init
cd $CMSSW_BASE/src 
# git clone -b HHWWgg_dev https://github.com/atishelmanch/flashgg 
git clone -b abe_HHWWgg_dev_fullyHad git@github.com:ram1123/flashgg.git
source flashgg/setup_flashgg.sh
```

If everything now looks reasonable, you can build:

```bash
cd $CMSSW_BASE/src
scram b -j 4
```

To access grid files to run the tagger on, you must run the following commands:

```bash
cmsenv
voms-proxy-init --voms cms --valid 168:00
```

after the `voms` command, you should receive an output similar to:

```bash
Created proxy in /tmp/x509up_u95168
```

to set this proxy to your `X509_USER_PROXY` environment variable for the example above, simply use the command:

```bash
. proxy.sh x509up_u95168
```

where `x590up_u95168` would be replaced by whatever your proxy name is. 


# General Workflow/Information of Package

## HHWWgg Tagger

The HHWWgg Tagger is developed to tag events as coming from the HH->WWgg process, and is compatible with workspaceStd in order to include the standard systematics workflow, 
and if desired to include tagging of other flashgg tags on the same events. 

### Running Locally 

The HHWWgg Tagger can be run locally on signal (with 2017 metaConditions) with:

```bash
cmsRun Systematics/test/workspaceStd.py metaConditions=MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1.json campaign=HHWWgg_v2-6 dataset=ggF_X600_HHWWgg_qqlnu doHHWWggTag=1 HHWWggTagsOnly=1 maxEvents=500 doSystematics=0 dumpWorkspace=0 dumpTrees=1 useAAA=1 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1
```

