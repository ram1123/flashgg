# General Introduction

# Code setup

```bash
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_5_0 
cd CMSSW_10_5_0/src
cmsenv
git cms-init
cd $CMSSW_BASE/src 
git clone -b HHWWgg_dev https://github.com/atishelmanch/flashgg 
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

after the voms command, you should receive an output similar to:

    Created proxy in /tmp/x509up_u95168

to set this proxy to your `X509_USER_PROXY` environment variable for the example above, simply use the command:

```bash
. proxy.sh x509up_u95168
```

where `x590up_u95168` would be replaced by whatever your proxy name is. 


# General Workflow/Information of Package

## 