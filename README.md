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

The HHWWgg Tagger is developed to tag events as coming from the `HH->WWgg` process, and is compatible with workspaceStd in order to include the standard systematics workflow, and if desired to include tagging of other flashgg tags on the same events. 

### Running Locally 

The HHWWgg Tagger can be run locally on signal (with 2017 metaConditions) with:

#### On Data

```bash
cmsRun Systematics/test/workspaceStd.py metaConditions=MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1.json campaign=Era2017_RR-31Mar2018_v2 dataset=/DoubleEG/spigazzi-Era2017_RR-31Mar2018_v2-legacyRun2FullV1-v0-Run2017B-31Mar2018-v1-d9c0c6cde5cc4a64343ae06f842e5085/USER doHHWWggTag=1 HHWWggTagsOnly=1 maxEvents=500 doSystematics=0 dumpWorkspace=1 dumpTrees=1 useAAA=1 processId=Data processType=Data doHHWWggTagCutFlow=0 saveHHWWggFinalStateVars=0
```

1. With default conditions:
   ```bash
   cmsRun Systematics/test/workspaceStd.py metaConditions=MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1.json campaign=HHWWgg_v2-6 dataset=ggF_X600_HHWWgg_qqlnu doHHWWggTag=1 HHWWggTagsOnly=1 maxEvents=500 doSystematics=0 dumpWorkspace=0 dumpTrees=1 useAAA=1 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1
   ```

1. With FullyHad benchmark11
   2. Without Systematics:
      ```bash
      cmsRun Systematics/test/workspaceStd.py metaConditions=MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1.json campaign=HHWWgg_v2-6 dataset=ggF_X600_HHWWgg_qqlnu doHHWWggFullyHadTag=1 HHWWggTagsOnly=1 maxEvents=500 doSystematics=0 dumpWorkspace=1 dumpTrees=0 useAAA=1 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1
      ```
   2. With Systematics:
      ```bash
      cmsRun Systematics/test/workspaceStd.py metaConditions=MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1.json campaign=HHWWgg_v2-6 dataset=ggF_X600_HHWWgg_qqlnu doHHWWggFullyHadTag=1 HHWWggTagsOnly=1 maxEvents=500 doSystematics=1 dumpWorkspace=1 dumpTrees=1 useAAA=1 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1
      ```

      ```bash
      cmsRun Systematics/test/workspaceStd.py metaConditions=MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1.json campaign=rasharma-HHWWgg_v2-2_Test_94X_mc2017 dataset=ggF_X250_WWgg_qqlnugg doHHWWggFullyHadTag=1 HHWWggTagsOnly=1 maxEvents=500 doSystematics=0 dumpWorkspace=0 dumpTrees=1 useAAA=1 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1
      ```

      1. **Signal**: 
      ```bash
        cmsRun Systematics/test/workspaceStd.py metaConditions=MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1.json campaign=rasharma-HHWWgg_v2-2_Test_94X_mc2017 dataset=ggF_HHWWgg_qqqq_node11 doHHWWggTag=1 HHWWggTagsOnly=1 maxEvents=500 doSystematics=0 dumpWorkspace=1 dumpTrees=1 useAAA=1 doHHWWggTagCutFlow=0 saveHHWWggFinalStateVars=0
      ```


## Condor Job

### On Data

```bash
. HHWWgg_Run_Jobs.sh --labelName HHWWgg_2017_Data_Trees --nEvents all --output /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/HHWWgg_5July_v3/ --json Taggers/test/HHWWgg_2017_Data_All/HHWWgg_Data_All_2017.json --condorQueue longlunch --year 2017 -g -c -t -w -s
```

### On Signal
```bash
. HHWWgg_Run_Jobs.sh --labelName GluGluToHHTo_WWgg_qqqq_node11 --nEvents all --output /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/HHWWgg_5July_v3/ --json Taggers/test/HHWWgg_v2-6/HHWWggFullyHad.json  --condorQueue longlunch --year 2017 -g -c -t -w -s
```

### Hadd script
```bash
. HHWWgg_Process_Files.sh --inFolder  GluGluToHHTo_WWgg_qqqq_node11 --outFolder GluGluToHHTo_WWgg_qqqq_node11_Hadded -s --signalType EFT
```

#### On data

Need to run both commands one by one:
- Hadd-ed different RunEra

   ```bash
   . HHWWgg_Process_Files.sh --inFolder HHWWgg_2017_Data_Trees --outFolder HHWWgg_2017_Data_Trees_Hadded -d
   ```

- Hadd-ed everything

   ```bash
   . HHWWgg_Process_Files.sh --inFolder HHWWgg_2017_Data_Trees_Hadded --outFolder HHWWgg_2017_Data_Trees_Hadded_Combined -d -c
   ```

### Important point

1. First point

```bash
$cat Taggers/test/HHWWgg_v2-6/HHWWggFullyHad.json 
{
   "processes" : {
        "GluGluToHHTo_WWgg_qqqq_node11" : [ "/ggF_node11_HHWWgg_qqqq"  ]

},
    "cmdLine" : "campaign=rasharma-HHWWgg_v2-2_Test_94X_mc2017 targetLumi=1e+3 useAAA=1 useEOS=0 puTarget=6.245e-06,...,8.814e-12"
```

The name `ggF_node11_HHWWgg_qqqq` should be defined in `MetaData/data/cross_sections.json`. Something like:

```json
   "ggF_node11_HHWWgg_qqqq": {
  "xs": 0.001
 },
```

## Check the json file status

```bash
fggManageSamples.py -C rasharma-HHWWgg_v2-2_Test_94X_mc2017 
```
This command will give the summary status of the DAS sample present with campaign name `rasharma-HHWWgg_v2-2_Test_94X_mc2017`.