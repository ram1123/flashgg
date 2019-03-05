# Can call created functions here to get variables 
# For example: HHWWggCandidate.cc defines function:
# std::vector<float> HHWWggCandidate::CosThetaAngles()
#
# In this file, you can call it with: 
# CosThetaAngles()[1]
#
# This will return whatever variable you're interested in 

pho_variables = [
      "npho                    := phoVector.size()",
      "pho1_pt                 := phoP4Corrected[0].pt()",
      "pho2_pt                 := phoP4Corrected[1].pt()",
      "pho3_pt                 := ? phoP4Corrected.size() > 2 ? phoP4Corrected[2].pt() : -1",
      "pho4_pt                 := ? phoP4Corrected.size() > 3 ? phoP4Corrected[3].pt() : -1",
      "pho1_eta                := phoP4Corrected[0].eta()",
      "pho2_eta                := phoP4Corrected[1].eta()",
      "pho3_eta                := ? phoP4Corrected.size() > 2 ? phoP4Corrected[2].eta() : -999",
      "pho4_eta                := ? phoP4Corrected.size() > 3 ? phoP4Corrected[3].eta() : -999",
      "pho1_phi                := phoP4Corrected[0].phi()",
      "pho2_phi                := phoP4Corrected[1].phi()",
      "pho3_phi                := ? phoP4Corrected.size() > 2 ? phoP4Corrected[2].phi() : -999",
      "pho4_phi                := ? phoP4Corrected.size() > 3 ? phoP4Corrected[3].phi() : -999",
      "pho1_energy             := phoP4Corrected[0].energy()",
      "pho2_energy             := phoP4Corrected[1].energy()",
      "pho3_energy             := ? phoP4Corrected.size() > 2 ? phoP4Corrected[2].energy() : -1",
      "pho4_energy             := ? phoP4Corrected.size() > 3 ? phoP4Corrected[3].energy() : -1",
      "pho1_old_r9             := phoP4Corrected[0].old_r9()",
      "pho2_old_r9             := phoP4Corrected[1].old_r9()",
      "pho3_old_r9             := ? phoP4Corrected.size() > 2 ? phoP4Corrected[2].old_r9() : -999",
      "pho4_old_r9             := ? phoP4Corrected.size() > 3 ? phoP4Corrected[3].old_r9() : -999",
      "pho1_full5x5_r9         := phoP4Corrected[0].full5x5_r9()",
      "pho2_full5x5_r9         := phoP4Corrected[1].full5x5_r9()",
      "pho3_full5x5_r9         := ? phoP4Corrected.size() > 2 ? phoP4Corrected[2].full5x5_r9() : -999",
      "pho4_full5x5_r9         := ? phoP4Corrected.size() > 3 ? phoP4Corrected[3].full5x5_r9() : -999",
      "pho1_EGMVA              := phoP4Corrected[0].userFloat('PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values')",
      "pho2_EGMVA              := phoP4Corrected[1].userFloat('PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values')",
      "pho3_EGMVA              := ? phoP4Corrected.size() > 2 ? phoP4Corrected[2].userFloat('PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values') : -999",
      "pho4_EGMVA              := ? phoP4Corrected.size() > 3 ? phoP4Corrected[3].userFloat('PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values') : -999",
      "pho1_MVA                := pho1_MVA",
      "pho2_MVA                := pho2_MVA",
      "pho3_MVA                := pho3_MVA",
      "pho4_MVA                := pho4_MVA",
      "pho1_match              := phoP4Corrected[0].genMatchType()",
      "pho2_match              := phoP4Corrected[1].genMatchType()",
      "pho3_match              := ? phoP4Corrected.size() > 2 ? phoP4Corrected[2].genMatchType() : -999",
      "pho4_match              := ? phoP4Corrected.size() > 3 ? phoP4Corrected[3].genMatchType() : -999",
      "pho1_pixelseed          := phoP4Corrected[0].hasPixelSeed()",
      "pho2_pixelseed          := phoP4Corrected[1].hasPixelSeed()",
      "pho3_pixelseed          := ? phoP4Corrected.size() > 2 ? phoP4Corrected[2].hasPixelSeed() : -999",
      "pho4_pixelseed          := ? phoP4Corrected.size() > 3 ? phoP4Corrected[3].hasPixelSeed() : -999",
      "pho1_electronveto       := phoP4Corrected[0].passElectronVeto()",
      "pho2_electronveto       := phoP4Corrected[1].passElectronVeto()",
      "pho3_electronveto       := ? phoP4Corrected.size() > 2 ? phoP4Corrected[2].passElectronVeto() : -999",
      "pho4_electronveto       := ? phoP4Corrected.size() > 3 ? phoP4Corrected[3].passElectronVeto() : -999",
      "pho12_dR                := deltaR( phoP4Corrected[0].eta, phoP4Corrected[0].phi, phoP4Corrected[1].eta, phoP4Corrected[1].phi )",
      "pho13_dR                := ? phoP4Corrected.size() > 2 ? deltaR( phoP4Corrected[0].eta, phoP4Corrected[0].phi, phoP4Corrected[2].eta, phoP4Corrected[2].phi ) : -1",
      "pho14_dR                := ? phoP4Corrected.size() > 3 ? deltaR( phoP4Corrected[0].eta, phoP4Corrected[0].phi, phoP4Corrected[3].eta, phoP4Corrected[3].phi ) : -1",
      "pho23_dR                := ? phoP4Corrected.size() > 2 ? deltaR( phoP4Corrected[1].eta, phoP4Corrected[1].phi, phoP4Corrected[2].eta, phoP4Corrected[2].phi ) : -1",
      "pho24_dR                := ? phoP4Corrected.size() > 3 ? deltaR( phoP4Corrected[1].eta, phoP4Corrected[1].phi, phoP4Corrected[3].eta, phoP4Corrected[3].phi ) : -1",
      "pho34_dR                := ? phoP4Corrected.size() > 3 ? deltaR( phoP4Corrected[2].eta, phoP4Corrected[2].phi, phoP4Corrected[3].eta, phoP4Corrected[3].phi ) : -1",
      "pho12_m                 := HHWWggPho12.mass",
      "pho13_m                 := HHWWggPho13.mass",
      "pho14_m                 := HHWWggPho14.mass",
      "pho23_m                 := HHWWggPho23.mass",
      "pho24_m                 := HHWWggPho24.mass",
      "pho34_m                 := HHWWggPho34.mass",
      "dZ                      := genVertex.z() - vertex.z()"
    ]

dipho_variables = [
    "dp1_mass                := HHWWggDiPho1.mass() ",
    "dp2_mass                := HHWWggDiPho2.mass() ",
    "avg_dp_mass             := (HHWWggDiPho1.mass()+HHWWggDiPho2.mass())/2",
    "dp1_pt                  := HHWWggDiPho1.pt() ",
    "dp2_pt                  := HHWWggDiPho2.pt() ",
    "dp1_eta                 := HHWWggDiPho1.eta() ",
    "dp2_eta                 := HHWWggDiPho2.eta() ",
    "dp1_phi                 := HHWWggDiPho1.phi() ",
    "dp2_phi                 := HHWWggDiPho2.phi() ",
    "dp1_e                   := HHWWggDiPho1.energy() ",
    "dp2_e                   := HHWWggDiPho2.energy() ",
    "dp1_dr                  := deltaR(HHWWggDiPho1_Pho1.eta(), HHWWggDiPho1_Pho1.phi(), HHWWggDiPho1_Pho2.eta(), HHWWggDiPho1_Pho2.phi())",
    "dp2_dr                  := deltaR(HHWWggDiPho2_Pho1.eta(), HHWWggDiPho2_Pho1.phi(), HHWWggDiPho2_Pho2.eta(), HHWWggDiPho2_Pho2.phi())",
    "dp1_p1i                 := HHWWggDiPho1_iPho1",
    "dp1_p2i                 := HHWWggDiPho1_iPho2",
    "dp2_p1i                 := HHWWggDiPho2_iPho1",
    "dp2_p2i                 := HHWWggDiPho2_iPho2",
    "dp1_PtoverMass          := HHWWggDiPho1.mass()/HHWWggDiPho1.pt()",
    "dp2_PtoverMass          := HHWWggDiPho2.mass()/HHWWggDiPho2.pt()",
    #"absCosThetaStar_CS      := abs(getCosThetaStar_CS(6500))",
    #"absCosTheta_pho_a1      := abs(CosThetaAngles()[1])",
    #"absCosTheta_pho_a2      := abs(CosThetaAngles()[0])"
]

tp_variables = [
    "tp_mass                 := HHWWggFourVect.mass()",
    "tp_pt                   := HHWWggFourVect.pt()",
    "tp_eta                  := HHWWggFourVect.eta()",
    "tp_phi                  := HHWWggFourVect.phi()",
    "tp_PtoverMass           := HHWWggFourVect.mass()/HHWWggFourVect.pt()"
]

ws_variables = [
   "tp_mass                 := HHWWggFourVect.mass()",
   "dp1_mass                := HHWWggDiPho1.mass() ",
   "dp2_mass                := HHWWggDiPho2.mass() ",
   "avg_dp_mass             := (HHWWggDiPho1.mass()+HHWWggDiPho2.mass())/2",
   "dZ                      := 0"
]

abe_variables = [

    # One entry per event 
    # Plot from header file utils 
    "n_ps_dipho                   := diphoVector.size()",
    "ps_dipho_mass                := ? Abe_HHWWggDiPho.mass() != 0 ? Abe_HHWWggDiPho.mass() : -999 ",
    "elec1_pt                     := elec1.pt()", # Leading pT
    "elec2_pt                     := elec2.pt()",  # Subleading pT 
    "muon1_pt                     := muon1.pt()", 
    "MET                          := MET_fourvec.pt()",
    "gen_lepton_pt                := gen_lepton_pt",
    "gen_neutrino_pt              := gen_neutrino_pt"
    # Change to gen_electron_pt
    # gen_muon_pt 
    #"Transverse_W_Mass            := (MET + elec1).Mt()"
    #"W1_TM                        := W1_TM",
    #"W2_TM                        := (MET + elec1).Mt()"
    #"MET                          := theMET.mPt()"
    #"MET                          := theMET.mPt()"
    # "theMETcorpt                  := theMETcorpt"
    #"ps_dipho_mass                :=                      "
    #"dipho_mass	              := "
    #"npho_                    := phoVector.size()"

    # GEN variables 

]
