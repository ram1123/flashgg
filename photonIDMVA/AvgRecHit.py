from ROOT import *

gROOT.SetBatch(1) # Do not plot upon draw statement 

def DrawSave(h_tmp, name):
    output_path = "/eos/user/a/atishelm/www/PhotonIDStudy/RecHits/"
    c1 = TCanvas()
    h_tmp.Draw("COLZ1")
    # h_tmp.GetZaxis().SetRangeUser(0,)
    c1.SaveAs(output_path + name)
    h_tmp.SaveAs(output_path + name[:-4] + ".root")

Avg_sigRecHitMapEB = TH2F("Avg_sigRecHitMapEB","Avg_sigRecHitMapEB",10,0,10,10,0,10)  
Avg_sigRecHitMapEB.GetXaxis().SetTitle("DOF1");  
Avg_sigRecHitMapEB.GetYaxis().SetTitle("DOF2")

# Avg_sigRecHitMapEEp = TH2F("Avg_sigRecHitMapEEp","Avg_sigRecHitMapEEp",10,0,10,10,0,10)  
# Avg_sigRecHitMapEEp.GetXaxis().SetTitle("DOF1");  
# Avg_sigRecHitMapEEp.GetYaxis().SetTitle("DOF2")

# Avg_sigRecHitMapEEm = TH2F("Avg_sigRecHitMapEEm","Avg_sigRecHitMapEEm",10,0,10,10,0,10)  
# Avg_sigRecHitMapEEm.GetXaxis().SetTitle("DOF1");  
# Avg_sigRecHitMapEEm.GetYaxis().SetTitle("DOF2")

Avg_FakeRecHitMapEB = TH2F("Avg_FakeRecHitMapEB","Avg_FakeRecHitMapEB",10,0,10,10,0,10)  
Avg_FakeRecHitMapEB.GetXaxis().SetTitle("DOF1");  
Avg_FakeRecHitMapEB.GetYaxis().SetTitle("DOF2")

# Avg_FakeRecHitMapEEp = TH2F("Avg_FakeRecHitMapEEp","Avg_FakeRecHitMapEEp",10,0,10,10,0,10)  
# Avg_FakeRecHitMapEEp.GetXaxis().SetTitle("DOF1");  
# Avg_FakeRecHitMapEEp.GetYaxis().SetTitle("DOF2")

# Avg_FakeRecHitMapEEm = TH2F("Avg_FakeRecHitMapEEm","Avg_FakeRecHitMapEEm",10,0,10,10,0,10)  
# Avg_FakeRecHitMapEEm.GetXaxis().SetTitle("DOF1");  
# Avg_FakeRecHitMapEEm.GetYaxis().SetTitle("DOF2")

input_file = "/afs/cern.ch/work/a/atishelm/21JuneFlashgg/CMSSW_10_5_0/src/flashgg/output_numEvent1000.root"

f = TFile(input_file)
Signal_t = f.Get("photonViewDumper/trees/promptPhotons") 
Fake_t = f.Get("photonViewDumper/trees/fakePhotons") 

# print't.Print() = ',t.Print()

num_rec_hits = 100

DOF1s = []
DOF2s = []
DOF3s = []
irechitEs = []
N_SigEB_rechits = 0 

width = 10 

# 10x10 hit window 
# https://stackoverflow.com/questions/9459337/assign-value-to-an-individual-cell-in-a-two-dimensional-python-array
# doing [[0]*10]*10 copies the same list object 10 times and then if you update one, it updates all, becasue they are all copies of the same instance 

SigEB_rechits = [ [ 0 for i in range(width) ] for j in range(width) ]
NSigEB_rechits = [ [ 0 for i in range(width) ] for j in range(width) ]
AvgSigEB_rechits = [ [ 0 for i in range(width) ] for j in range(width) ]

# print'SigEB_rechits = ',SigEB_rechits

#"""
for i in range(num_rec_hits):
    tmp_name = "DOF1s_"
    tmp_name += str(i)
    DOF1s.append(tmp_name)

    tmp_name = "DOF2s_"
    tmp_name += str(i)
    DOF2s.append(tmp_name)

    tmp_name = "DOF3s_"
    tmp_name += str(i)
    DOF3s.append(tmp_name)

    tmp_name = "recHit_"
    tmp_name += str(i)
    irechitEs.append(tmp_name)

h1 = TH1F("h1","h1",20,-10,10)
h2 = TH1F("h2","h2",20,-10,10)
rechits = TH1F("rechits","rechits",50,0,50)
rel_position = TH2F("rel_position","rel_position",10,0,10,10,0,10)

ie = 0
for entry in Signal_t:

    # Given entry has one photon object with some DOF1,2,3 
    pho_object_DOF1 = eval("Signal_t.DOF1")
    pho_object_DOF2 = eval("Signal_t.DOF2")
    pho_object_DOF3 = eval("Signal_t.DOF3")

    for rh in range(0,num_rec_hits):
        # Given entry has many rec hits to make map surrounding one photon object position 
    
        DOF1_eval = "Signal_t." + DOF1s[rh]
        DOF2_eval = "Signal_t." + DOF2s[rh]
        DOF3_eval = "Signal_t." + DOF3s[rh]
        recHitEeval = "Signal_t." + irechitEs[rh]

        DOF1 = eval(DOF1_eval)
        DOF2 = eval(DOF2_eval)
        DOF3 = eval(DOF3_eval)
        recHitE = eval(recHitEeval)

        if DOF1 is not -9999:
            # if DOF3 == -1:  
                # sigRecHitMapEEm.Fill(DOF1,DOF2,recHitE)
            if DOF3 == 0:  
                # plot relative to seed DOF1, DOF2 
                # 1,4 is getting filled if the second rel dof is 4 ... 
                rel_dof1 = int((DOF1 - pho_object_DOF1)) + 4  
                rel_dof2 = int((DOF2 - pho_object_DOF2)) + 4  
                # print'rel_dof1 = ',rel_dof1
                # print'rel_dof2 = ',rel_dof2
                # print'recHitE = ',recHitE
                h1.Fill(rel_dof1)
                h2.Fill(rel_dof2)
                # print'(rel_dof1,rel_dof2) = (',rel_dof1,',',rel_dof2,')'
                rel_position.Fill(rel_dof1,rel_dof2)
                # print'SigEB_rechits[rel_dof1][rel_dof2] = ',SigEB_rechits[rel_dof1][rel_dof2]
                # print'before adding this: SigEB_rechits[1][4] = ',SigEB_rechits[1][4]
                # SigEB_rechits[rel_dof1][rel_dof2] += recHitE  
                SigEB_rechits[rel_dof1][rel_dof2] += recHitE  
                # print'before adding this: SigEB_rechits[1][4] = ',SigEB_rechits[1][4]

                NSigEB_rechits[rel_dof1][rel_dof2] += 1 

                rechits.Fill(recHitE)

                # ie += 1 
                

                # SigEB_rel_dof1s.append(rel_dof1)
                # SigEB_rel_dof2s.append(rel_dof2)
                # SigEB_rechits.append(rechit)
                # Avg_sigRecHitMapEB.Fill(rel_dof1,rel_dof2,rechit)
            # elif DOF3 == 1:  
                # sigRecHitMapEEp.Fill(DOF1,DOF2,recHitE) 
        # if ie == 1: break 


DrawSave(rechits,"rechits.png")
DrawSave(rel_position,"rel_position.png")

print'SigEB_rechits = ',SigEB_rechits
print'NSigEB_rechits = ',NSigEB_rechits

for i in range(10):
    for j in range(10):

        cell_totalRecHit = float(SigEB_rechits[i][j])
        cell_NRecHit = float(NSigEB_rechits[i][j])

        if cell_NRecHit == 0:
            avgRecHit = 0
        else: 
            avgRecHit = cell_totalRecHit / cell_NRecHit

        AvgSigEB_rechits[i][j] = avgRecHit 

        # if j == 5:
        #     print'cell_totalRecHit = ',cell_totalRecHit
        #     print'cell_NRecHit = ',cell_NRecHit
        #     print'avgRecHit = ',avgRecHit


# print 'AvgSigEB_rechits = ',AvgSigEB_rechits

h3 = TH1F("h3","h3",100,0,10)

for i in range(10):
    for j in range(10):
        # print 'i = ',i
        # print 'j = ',j
        avg_rechit = AvgSigEB_rechits[i][j]
        # print 'avg_rechit = ',avg_rechit
        h3.Fill(avg_rechit)
        Avg_sigRecHitMapEB.Fill(i,j,avg_rechit)

DrawSave(h1, "reldof1.png")
DrawSave(h2, "reldof2.png")
DrawSave(h3, "avgrechitdistro.png")
DrawSave(Avg_sigRecHitMapEB, "Avg_sigRecHitMapEB.png")

# h1.SaveAs("rel_dof1.root")
# h2.SaveAs("rel_dof2.root")

# print'sigRecHitMapEEm.GetEntries() = ',sigRecHitMapEEm.GetEntries()

# DrawSave(sigRecHitMapEEm, "sigRecHitMapEEm.png")
# DrawSave(sigRecHitMapEB, "sigRecHitMapEB.png")
# DrawSave(sigRecHitMapEEp, "sigRecHitMapEEp.png")

FakeEB_rechits = [ [ 0 for i in range(width) ] for j in range(width) ]
NFakeEB_rechits = [ [ 0 for i in range(width) ] for j in range(width) ]
AvgFakeEB_rechits = [ [ 0 for i in range(width) ] for j in range(width) ]

for entry in Fake_t:
 
    # Given entry has one photon object with some DOF1,2,3 
    pho_object_DOF1 = eval("Fake_t.DOF1")
    pho_object_DOF2 = eval("Fake_t.DOF2")
    pho_object_DOF3 = eval("Fake_t.DOF3")

    for rh in range(0,num_rec_hits):
    
        DOF1_eval = "Fake_t." + DOF1s[rh]
        DOF2_eval = "Fake_t." + DOF2s[rh]
        DOF3_eval = "Fake_t." + DOF3s[rh]
        recHitEeval = "Fake_t." + irechitEs[rh]

        DOF1 = eval(DOF1_eval)
        DOF2 = eval(DOF2_eval)
        DOF3 = eval(DOF3_eval)
        recHitE = eval(recHitEeval)

        # if DOF1 is not -9999:
        #     if DOF3 == -1:  FakeRecHitMapEEm.Fill(DOF1,DOF2,recHitE)
        #     elif DOF3 == 0:  FakeRecHitMapEB.Fill(DOF1,DOF2,recHitE) 
        #     elif DOF3 == 1:  FakeRecHitMapEEp.Fill(DOF1,DOF2,recHitE) 

        if DOF1 is not -9999:
            # if DOF3 == -1:  
                # sigRecHitMapEEm.Fill(DOF1,DOF2,recHitE)
            if DOF3 == 0:  
                # plot relative to seed DOF1, DOF2 
                rel_dof1 = int((DOF1 - pho_object_DOF1)) + 4  
                rel_dof2 = int((DOF2 - pho_object_DOF2)) + 4  
                
                FakeEB_rechits[rel_dof1][rel_dof2] += recHitE  
                NFakeEB_rechits[rel_dof1][rel_dof2] += 1 

                # SigEB_rel_dof1s.append(rel_dof1)
                # SigEB_rel_dof2s.append(rel_dof2)
                # SigEB_rechits.append(rechit)
                # Avg_sigRecHitMapEB.Fill(rel_dof1,rel_dof2,rechit)
            # elif DOF3 == 1:  
                # sigRecHitMapEEp.Fill(DOF1,DOF2,recHitE) 
        # if ie == 1: break 


print'FakeEB_rechits = ',FakeEB_rechits
print'NFakeEB_rechits = ',NFakeEB_rechits

for i in range(10):
    for j in range(10):

        cell_totalRecHit = float(FakeEB_rechits[i][j])
        cell_NRecHit = float(NFakeEB_rechits[i][j])

        if cell_NRecHit == 0:
            avgRecHit = 0
        else: 
            avgRecHit = cell_totalRecHit / cell_NRecHit

        AvgFakeEB_rechits[i][j] = avgRecHit 

        # if j == 5:
        #     print'cell_totalRecHit = ',cell_totalRecHit
        #     print'cell_NRecHit = ',cell_NRecHit
        #     print'avgRecHit = ',avgRecHit


# print 'AvgSigEB_rechits = ',AvgSigEB_rechits

# h3 = TH1F("h3","h3",100,0,10)

for i in range(10):
    for j in range(10):
        # print 'i = ',i
        # print 'j = ',j
        avg_rechit = AvgFakeEB_rechits[i][j]
        # print 'avg_rechit = ',avg_rechit
        h3.Fill(avg_rechit)
        Avg_FakeRecHitMapEB.Fill(i,j,avg_rechit)

DrawSave(Avg_FakeRecHitMapEB, "Avg_FakeRecHitMapEB.png")

# DrawSave(FakeRecHitMapEEm, "FakeRecHitMapEEm.png")
# DrawSave(FakeRecHitMapEB, "FakeRecHitMapEB.png")
# DrawSave(FakeRecHitMapEEp, "FakeRecHitMapEEp.png")

# Want average 10x10 window for sigRecHitMapEB, FakeRecHitMapEB and side by side 
