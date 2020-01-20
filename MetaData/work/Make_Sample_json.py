# Abe Tishelman-Charny
# 5 December 2019 
# The purpose of this module is to create an MC_Configs.json file to run with . main.sh

import os 
# import subprocess

outputName = 'sampleJSON.json' # output json file path 

masses = ['260', '270', '280', '300', '320', '350', '400', '500', '550', '600', '650', '700', '800', '850', '900', '1000']

# Begin writing file 
Sample_JSON = '{' 
Sample_JSON += '\n' 
Sample_JSON += '   "data" : [],\n'
Sample_JSON += '   "sig" : [],\n'
Sample_JSON += '   "bkg"  : ['

Sample_JSON += '"/ggF_X250_WWgg_qqlnugg/atishelm-100000events_wPU_MINIAOD-5f646ecd4e1c7a39ab0ed099ff55ceb9/USER",\n'
Sample_JSON += '	"/ggF_X750_WWgg_qqlnugg/atishelm-100000events_wPU_MINIAOD-5f646ecd4e1c7a39ab0ed099ff55ceb9/USER",\n'
Sample_JSON += '	"/ggF_SM_WWgg_qqlnugg/atishelm-100000events_wPU_MINIAOD-5f646ecd4e1c7a39ab0ed099ff55ceb9/USER",\n'
Sample_JSON += '	"/ggF_X1250_WWgg_qqlnugg/atishelm-100000events_wPU_MINIAOD-5f646ecd4e1c7a39ab0ed099ff55ceb9/USER",'

for im,mass in enumerate(masses):
	Sample_JSON += """
	"/ggF_X{mass}_WWgg_qqlnugg/atishelm-100000events_wPU_MINIAOD-5f646ecd4e1c7a39ab0ed099ff55ceb9/USER" """
	Sample_JSON = Sample_JSON.replace("{mass}",str(mass))
	if im is not len(masses)-1: 
		Sample_JSON += ',' # need comma separation 
		# Sample_JSON += '\n'

	else: continue # no comma at end of last object 

Sample_JSON += '	]\n'
Sample_JSON += '}'


# Sample_JSON += '\n]\n' # finish json 

with open(outputName, "w") as output:
		output.write(Sample_JSON) # write json file 

print 
print'[Make_Sample_json] - sampleJSON.json created'
print'[Make_Sample_json] - Make sure sampleJSON.json looks good before submitting !' 
print 
