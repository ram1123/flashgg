import os 

finalStateDict = {
  #"SL" : "2Qlnu",
  "FH" : "4Q",
  #"FL" : "2l2nu"
}

nodes = [node for node in range(1,13)]
nodes.append("SM")

for finalState in ["FH"]:
#for finalState in ["FL"]:
  finalStateParticles = finalStateDict[finalState]
  for node in nodes:
    command = "hadd %s_LO_2016_hadded/GluGluToHHTo2G%s_node_%s_2016.root %s_LO_2016/output_GluGluToHHTo2G%s_node_%s_*.root"%(finalState,finalStateParticles,node,finalState,finalStateParticles,node)
    print"==========COMMAND============"
    print command 
    os.system(command)
