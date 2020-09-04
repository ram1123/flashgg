# Motivation of this script: 
#fggManageSamples.py -C rasharma -V  94X_mc2017-RunIIFall18_ggF_nodes_ZZ import

import os

#if os.path.exists("MetaData/data/rasharma-94X_mc2017"):

Campaign = 'rasharma'
#Version = '94X_mc2017-RunIIFall18_Radion_Bqrk_wPU'
Version = '94X_mc2017-RunIIFall18_ggF_nodes_WW'

while True:
    # os.system("fggManageSamples.py -C rasharma -V 94X_mc2017-RunIIFall18_Radion_Bqrk_wPU  import")
    os.system("fggManageSamples.py -C "+Campaign+" -V "+Version+"  import")
    # if os.path.exists("MetaData/data/rasharma"):
    if os.path.exists("MetaData/data/"+Campaign):
        print "exits"
        break;
    else:
        print "="*51,"\n\n\n\nDir does not exists.. Continue the job...\n\n\n"

