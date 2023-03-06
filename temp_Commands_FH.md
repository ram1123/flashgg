# Commands summary

## Nov 2021




## bbgg all three years



### 2016

```bash
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2016/Signal/bbgg_2016.json -D -P -n 500 -d FHWW_bbgg_NLO_2016 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Signal/FHWW_bbgg_NLO_2016/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=0 useParentDataset=0 recalculatePDFWeights=0
```

```bash
python Systematics/scripts/resubmit_jobs.py --dir FHWW_bbgg_NLO_2016 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Signal/FHWW_bbgg_NLO_2016/
python Systematics/scripts/resubmit_jobs.py --dir FHWW_bbgg_NLO_2017 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/Signal/FHWW_bbgg_NLO_2017/
python Systematics/scripts/resubmit_jobs.py --dir FHWW_bbgg_NLO_2018 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Signal/FHWW_bbgg_NLO_2018/
```

### 2017
```bash
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2017/Signal/bbgg_2017.json -D -P -n 500 -d FHWW_bbgg_NLO_2017 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/Signal/FHWW_bbgg_NLO_2017/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q workday --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=0 useParentDataset=0 recalculatePDFWeights=0
```

### 2018

```bash
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2018/Signal/bbgg_2018.json -D -P -n 500 -d FHWW_bbgg_NLO_2018 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Signal/FHWW_bbgg_NLO_2018/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2018_RR-17Sep2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=0 useParentDataset=0 recalculatePDFWeights=0
```


## 2016

## Single

```bash
#NLO
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2016/Signal/FHZZ_NLO_2016.json -D -P -n 500 -d FHZZ_NLO_2016 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Signal/FHZZ_NLO_2016/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=1 useParentDataset=1 recalculatePDFWeights=1
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2016/Signal/FH_NLO_2016.json   -D -P -n 500 -d FHWW_NLO_2016 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Signal/FHWW_NLO_2016/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=1 useParentDataset=1 recalculatePDFWeights=1
#LO
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2016/Signal/FHZZ_LO_2016.json -D -P -n 500 -d FHZZ_LO_2016_noPdfWeight --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Signal/FHZZ_LO_2016_noPdfWeight/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=0 useParentDataset=0 recalculatePDFWeights=0
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2016/Signal/FH_LO_2016.json   -D -P -n 500 -d FHWW_LO_2016_noPdfWeight --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Signal/FHWW_LO_2016_noPdfWeight/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=0 useParentDataset=0 recalculatePDFWeights=0
```

**Resubmit commands**

```bash
python Systematics/scripts/resubmit_jobs.py --dir FHWW_NLO_2016 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Signal/FHWW_NLO_2016/
python Systematics/scripts/resubmit_jobs.py --dir FHZZ_NLO_2016 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Signal/FHZZ_NLO_2016/

python Systematics/scripts/resubmit_jobs.py --dir FHWW_LO_2016_noPdfWeight -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Signal/FHWW_LO_2016_noPdfWeight/
python Systematics/scripts/resubmit_jobs.py --dir FHZZ_LO_2016_noPdfWeight -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Signal/FHZZ_LO_2016_noPdfWeight/
```

**Hadd commands**

```bash
. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Signal/ --inFolder  FHWW_NLO_2016 --outFolder  FHWW_NLO_2016_Hadded -s --signalType NORES -t

. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Signal/ --inFolder  FHZZ_NLO_2016 --outFolder  FHZZ_NLO_2016_Hadded -s --signalType NORES -t
```

## Data

```bash
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2016/Data/HHWWgg_Data_All_2016.json   -D -P -n 500 -d Data_Trees_2016 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Data_Trees_2016/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2018/Data/HHWWgg_Data_All_2018.json   -D -P -n 500 -d Data_Trees_2018 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Data_Trees_2018/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2018_RR-17Sep2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1
```

**Resubmit commands**
```bash
python Systematics/scripts/resubmit_jobs.py --dir  Data_Trees_2016 -s  /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Data_Trees_2016
```

**Hadd Commands**

```bash
. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/ --inFolder Data_Trees_2016 --outFolder Data_Trees_2016_Hadded -d -t

. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/ --inFolder Data_Trees_2016_Hadded --outFolder  Data_Trees_2016_Hadded_Combined -d -c -t
```

## Single Higgs
```bash
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2016/Single_H/Higgs_bkg_2016_120125130.json -D -P -n 500 -d Single_H_2016 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Single_H_2016/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=1 useParentDataset=1 recalculatePDFWeights=1
```

**Resubmit commands**
```bash
python Systematics/scripts/resubmit_jobs.py --dir Single_H_2016 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Single_H_2016
```


**Hadd Commands**

```bash
. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/ --inFolder  Single_H_2016 --outFolder  FSingle_H_2016_Hadded -b -t
```

## Backgrounds


```bash
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2016/Backgrounds/Flashgg_bkg.json -D -P -n 500 -d Backgrounds_2016v2 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Backgrounds_2016v2/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=1 useParentDataset=1 recalculatePDFWeights=1
```

**Resubmit commands**
```bash
python Systematics/scripts/resubmit_jobs.py --dir Backgrounds_2016v2 -s  /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2016/Backgrounds_2016v2/
```

# 2017

## Data

```bash
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2017/Data/HHWWgg_Data_All_2017.json  -D -P -n 500 -d Data_Trees_2017 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/Data_Trees_2017/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1
```

**Resubmit commands**
```bash
python Systematics/scripts/resubmit_jobs.py --dir  Data_Trees_2017 -s  /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/Data_Trees_2017
```

**Hadd Commands**

```bash
. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/ --inFolder Data_Trees_2017 --outFolder Data_Trees_2017_Hadded -d -t

. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/ --inFolder Data_Trees_2017_Hadded --outFolder  Data_Trees_2017_Hadded_Combined -d -c -t
```



## Single Higgs
```bash
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2017/Single_H/Higgs_bkg_2017_120125130.json -D -P -n 500 -d Single_H_2017 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/Single_H_2017/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=1 useParentDataset=1 recalculatePDFWeights=1
```


## QCD Jobs

```bash
#fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2017/HHWWgg_QCD_HTBinned_50to100.json  -D -P -n 500 -d  QCD_HT50to100_2017 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/QCD_HT50to100_2017/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=500 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1-HHWWgg.json   doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=1 useParentDataset=1 recalculatePDFWeights=1
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2017/Backgrounds/HHWWgg_QCD_HTBinned.json  -D -P -n 500 -d  QCD_HTBinned_2017 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/QCD_HTBinned_2017/  -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q microcentury --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1-HHWWgg.json   doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=0 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1
```

**Resubmit Commands**

```bash
#python Systematics/scripts/resubmit_jobs.py --dir QCD_HT50to100_2017 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/QCD_HT50to100_2017/
python Systematics/scripts/resubmit_jobs.py --dir QCD_HTBinned_2017 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/QCD_HTBinned_2017/
```

. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/ --inFolder  QCD_HTBinned_2017 --outFolder  QCD_HTBinned_2017_Hadded -b -t


## Backgrounds


```bash
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2017/Backgrounds/Flashgg_bkg.json   -D -P -n 500 -d Background_2017     --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/Background_2017/     -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=0 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2017/Backgrounds/HHWWgg_bkg_v2.json -D -P -n 500 -d Background_2017ext1 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/Background_2017ext1/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=0 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2017/Backgrounds/HHWWgg_bkg_v3.json -D -P -n 500 -d Background_2017ext2 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/Background_2017ext2/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=0 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2017/Backgrounds/HHWWgg_bkg_v4.json -D -P -n 500 -d Background_2017ext3 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/Background_2017ext3/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=0 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2017/Backgrounds/HHWWgg_bkg_v5.json -D -P -n 500 -d Background_2017ext4 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/Background_2017ext4/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2017_RR-31Mar2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=0 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1
```

**Resubmit Commands**

```bash
python Systematics/scripts/resubmit_jobs.py --dir Background_2017 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/Background_2017/
```

**Hadd commands**

```bash
. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/ --inFolder  Background_2017 --outFolder  Background_2017_Hadded -b -t
. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/ --inFolder  Background_2017ext1 --outFolder  Background_2017ext1_Hadded -b -t
. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2017/ --inFolder  Background_2017ext4 --outFolder  Background_2017ext4_Hadded -b -t
```

# 2018


```bash
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2018/Signal/FH_NLO_2018.json    -D -P -n 500 -d FHWW_NLO_2018 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Signal/FHWW_NLO_2018/  -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2018_RR-17Sep2018_v1-HHWWgg.json  doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=1 useParentDataset=1 recalculatePDFWeights=1
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2018/Signal/FHZZ_NLO_2018.json  -D -P -n 500 -d FHZZ_NLO_2018 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Signal/FHZZ_NLO_2018/  -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2018_RR-17Sep2018_v1-HHWWgg.json  doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=1 useParentDataset=1 recalculatePDFWeights=1

fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2018/Signal/FH_LO_2018.json     -D -P -n 500 -d FHWW_LO_2018  --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Signal/FHWW_LO_2018/   -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2018_RR-17Sep2018_v1-HHWWgg.json  doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=1 useParentDataset=1 recalculatePDFWeights=1
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2018/Signal/FHZZ_LO_2018.json   -D -P -n 500 -d FHZZ_LO_2018  --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Signal/FHZZ_LO_2018/   -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2018_RR-17Sep2018_v1-HHWWgg.json  doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=1 useParentDataset=1 recalculatePDFWeights=1
```

**Resubmit commands**

```bash
python Systematics/scripts/resubmit_jobs.py --dir FHWW_LO_2018 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Signal/FHWW_LO_2018/
python Systematics/scripts/resubmit_jobs.py --dir FHZZ_LO_2018 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Signal/FHZZ_LO_2018/

python Systematics/scripts/resubmit_jobs.py --dir FHWW_NLO_2018 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Signal/FHWW_NLO_2018/
python Systematics/scripts/resubmit_jobs.py --dir FHZZ_NLO_2018 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Signal/FHZZ_NLO_2018/
```


**Hadd commands**

```bash
. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Signal/ --inFolder  FHWW_LO_2018 --outFolder  FHWW_LO_2018_Hadded -s --signalType NORES -t
. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Signal/ --inFolder  FHZZ_LO_2018 --outFolder  FHZZ_LO_2018_Hadded -s --signalType NORES -t
```


## Single_H

```bash
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2018/Single_H/Higgs_bkg_2018_120125130.json -D -P -n 500 -d Single_H_2018 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Single_H_2018/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=1 useParentDataset=1 recalculatePDFWeights=1
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2018/Single_H/Higgs_bkg_2018_missingTemp.json -D -P -n 500 -d Single_H_2018_ext1 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Single_H_2018_ext1/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2016_RR-17Jul2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1 doGranularJEC=1 doPdfWeights=1 useParentDataset=1 recalculatePDFWeights=1
```

**resubmit_jobs**

```bash
python Systematics/scripts/resubmit_jobs.py --dir Single_H_2018 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Single_H_2018/
python Systematics/scripts/resubmit_jobs.py --dir Single_H_2018_ext1 -s /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Single_H_2018_ext1/
```

**Hadd Commands**

```bash
. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/ --inFolder Single_H_2018Fixed      --outFolder  Single_H_2018Fixed_Hadded -b -t
```


## Data

```bash
fggRunJobs.py --load Taggers/test/HHWWgg/January_2021_Production/2018/Data/HHWWgg_Data_All_2018.json   -D -P -n 500 -d Data_Trees_2018 --stage-to=/eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Data_Trees_2018/ -x cmsRun Systematics/test/workspaceStd.py maxEvents=-1 -q tomorrow --no-use-tarball --no-copy-proxy metaConditions=$CMSSW_BASE/src/flashgg/MetaData/data/MetaConditions/Era2018_RR-17Sep2018_v1-HHWWgg.json doHHWWggTag=1 HHWWggTagsOnly=1 doSystematics=1 dumpTrees=1 dumpWorkspace=0 doHHWWggTagCutFlow=1 saveHHWWggFinalStateVars=1 copyInputMicroAOD=1
```

**Resubmit commands**
```bash
python Systematics/scripts/resubmit_jobs.py --dir  Data_Trees_2018 -s  /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/Data_Trees_2018
```


**Hadd Commands**

```bash
. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/ --inFolder Data_Trees_2018 --outFolder Data_Trees_2018_Hadded -d -t

. HHWWgg_Process_Files.sh --nTupleDir /eos/user/r/rasharma/post_doc_ihep/double-higgs/ntuples/January_2021_Production/2018/ --inFolder Data_Trees_2018_Hadded --outFolder  Data_Trees_2018_Hadded_Combined -d -c -t
```


