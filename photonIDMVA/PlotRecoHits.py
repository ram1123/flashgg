from ROOT import *

gROOT.SetBatch(1) # Do not plot upon draw statement 

def DrawSave(h_tmp, name):
    output_path = "/eos/user/a/atishelm/www/PhotonIDStudy/RecHits/"
    c1 = TCanvas()
    h_tmp.Draw("COLZ1")
    c1.SaveAs(output_path + name)

sigRecHitMapEB = TH2F("sigRecHitMapEB","sigRecHitMapEB",170,-85,85,360,0,360)  
sigRecHitMapEB.GetXaxis().SetTitle("i#eta");  
sigRecHitMapEB.GetYaxis().SetTitle("i#phi")

sigRecHitMapEEp = TH2F("sigRecHitMapEEp","sigRecHitMapEEp",100,0,100,100,0,100)  
sigRecHitMapEEp.GetXaxis().SetTitle("ix");  
sigRecHitMapEEp.GetYaxis().SetTitle("iy")

sigRecHitMapEEm = TH2F("sigRecHitMapEEm","sigRecHitMapEEm",100,0,100,100,0,100)  
sigRecHitMapEEm.GetXaxis().SetTitle("ix");  
sigRecHitMapEEm.GetYaxis().SetTitle("iy")

FakeRecHitMapEB = TH2F("FakeRecHitMapEB","FakeRecHitMapEB",170,-85,85,360,0,360)  
FakeRecHitMapEB.GetXaxis().SetTitle("i#eta");  
FakeRecHitMapEB.GetYaxis().SetTitle("i#phi")

FakeRecHitMapEEp = TH2F("FakeRecHitMapEEp","FakeRecHitMapEEp",100,0,100,100,0,100)  
FakeRecHitMapEEp.GetXaxis().SetTitle("ix");  
FakeRecHitMapEEp.GetYaxis().SetTitle("iy")

FakeRecHitMapEEm = TH2F("FakeRecHitMapEEm","FakeRecHitMapEEm",100,0,100,100,0,100)  
FakeRecHitMapEEm.GetXaxis().SetTitle("ix");  
FakeRecHitMapEEm.GetYaxis().SetTitle("iy")



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

for entry in Signal_t:
    for rh in range(0,num_rec_hits):
    
        DOF1_eval = "Signal_t." + DOF1s[rh]
        DOF2_eval = "Signal_t." + DOF2s[rh]
        DOF3_eval = "Signal_t." + DOF3s[rh]
        recHitEeval = "Signal_t." + irechitEs[rh]

        DOF1 = eval(DOF1_eval)
        DOF2 = eval(DOF2_eval)
        DOF3 = eval(DOF3_eval)
        recHitE = eval(recHitEeval)

        if DOF1 is not -9999:
            if DOF3 == -1:  sigRecHitMapEEm.Fill(DOF1,DOF2,recHitE)
            elif DOF3 == 0:  sigRecHitMapEB.Fill(DOF1,DOF2,recHitE) 
            elif DOF3 == 1:  sigRecHitMapEEp.Fill(DOF1,DOF2,recHitE) 

# print'sigRecHitMapEEm.GetEntries() = ',sigRecHitMapEEm.GetEntries()

DrawSave(sigRecHitMapEEm, "sigRecHitMapEEm.png")
DrawSave(sigRecHitMapEB, "sigRecHitMapEB.png")
DrawSave(sigRecHitMapEEp, "sigRecHitMapEEp.png")

for entry in Fake_t:
    for rh in range(0,num_rec_hits):
    
        DOF1_eval = "Fake_t." + DOF1s[rh]
        DOF2_eval = "Fake_t." + DOF2s[rh]
        DOF3_eval = "Fake_t." + DOF3s[rh]
        recHitEeval = "Fake_t." + irechitEs[rh]

        DOF1 = eval(DOF1_eval)
        DOF2 = eval(DOF2_eval)
        DOF3 = eval(DOF3_eval)
        recHitE = eval(recHitEeval)

        if DOF1 is not -9999:
            if DOF3 == -1:  FakeRecHitMapEEm.Fill(DOF1,DOF2,recHitE)
            elif DOF3 == 0:  FakeRecHitMapEB.Fill(DOF1,DOF2,recHitE) 
            elif DOF3 == 1:  FakeRecHitMapEEp.Fill(DOF1,DOF2,recHitE) 

DrawSave(FakeRecHitMapEEm, "FakeRecHitMapEEm.png")
DrawSave(FakeRecHitMapEB, "FakeRecHitMapEB.png")
DrawSave(FakeRecHitMapEEp, "FakeRecHitMapEEp.png")

# Want average 10x10 window for sigRecHitMapEB, FakeRecHitMapEB and side by side 

