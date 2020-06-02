grep -r Run2017E-31Mar2018-v1/190606_095510 ../Era2017_RR-31Mar2018_v2/datasets_8.json |awk '{print$2}' >e.log
sed -i 's#"/store#root://xrootd-cms.infn.it//store#g' e.log
sed -i 's#root",#root#g' e.log
cp e.log /eos/user/c/chuw/data/
