import AnaProd.HH_bbWW.baseline as AnaBaseline
import Common.BaselineSelection as CommonBaseline
from Corrections.Corrections import Corrections

loadTF = False
loadHHBtag = False
lepton_legs = [ "tau1", "tau2" ]

deepTauScores= ["rawDeepTau2017v2p1VSe","rawDeepTau2017v2p1VSmu",
            "rawDeepTau2017v2p1VSjet", "rawDeepTau2018v2p5VSe", "rawDeepTau2018v2p5VSmu",
            "rawDeepTau2018v2p5VSjet",
            "idDeepTau2017v2p1VSe", "idDeepTau2017v2p1VSjet", "idDeepTau2017v2p1VSmu",
            "idDeepTau2018v2p5VSe","idDeepTau2018v2p5VSjet","idDeepTau2018v2p5VSmu",
            "decayMode"]
Muon_observables = ["Muon_tkRelIso", "Muon_pfRelIso04_all"]
Electron_observables = ["Electron_mvaNoIso_WP90", "Electron_mvaIso_WP90", "Electron_pfRelIso03_all"]
JetObservables = ["particleNetAK4_B", "particleNetAK4_CvsB",
                "particleNetAK4_CvsL","particleNetAK4_QvsG","particleNetAK4_puIdDisc",
                "btagDeepFlavB","btagDeepFlavCvB","btagDeepFlavCvL", "bRegCorr", "bRegRes", "idbtagDeepFlavB",
                "btagPNetB", "btagPNetCvL", "btagPNetCvB", "btagPNetQvG", "btagPNetTauVJet", "PNetRegPtRawCorr", "PNetRegPtRawCorrNeutrino", "PNetRegPtRawRes"] # 2016]
JetObservablesMC = ["hadronFlavour","partonFlavour"]
FatJetObservables = ["area", "btagCSVV2", "btagDDBvLV2", "btagDeepB", "btagHbb", "deepTagMD_HbbvsQCD",
                     "deepTagMD_ZHbbvsQCD", "deepTagMD_ZbbvsQCD", "deepTagMD_bbvsLight", "deepTag_H",
                     "jetId", "msoftdrop", "nBHadrons", "nCHadrons", "nConstituents","rawFactor",
                      "particleNetMD_QCD", "particleNetMD_Xbb", "particleNet_HbbvsQCD", "particleNet_mass", # 2018
                     "particleNet_QCD","particleNet_XbbVsQCD", # 2016
                     "particleNetLegacy_QCD", "particleNetLegacy_Xbb", "particleNetLegacy_mass", # 2016
                     "particleNetWithMass_QCD", "particleNetWithMass_HbbvsQCD", "particleNet_massCorr" # 2016
                     ]

# in this PR https://github.com/cms-sw/cmssw/commit/17457a557bd75ab479dfb78013edf9e551ecd6b7, particleNet MD have been removed therefore we will switch to take


# # New ParticleNet trainings are not available in MiniAOD until Run3 13X
'''
particleNetWithMass_QCD YES
particleNetWithMass_TvsQCD NO
particleNetWithMass_WvsQCD NO
particleNetWithMass_ZvsQCD NO
particleNetWithMass_H4qvsQCD NO
particleNetWithMass_HbbvsQCD YES
particleNetWithMass_HccvsQCD NO
particleNet_QCD YES
particleNet_QCD2HF NO
particleNet_QCD1HF NO
particleNet_QCD0HF NO
particleNet_massCorr YES
particleNet_XbbVsQCD YES
particleNet_XccVsQCD NO
particleNet_XqqVsQCD NO
particleNet_XggVsQCD NO
particleNet_XttVsQCD NO
particleNet_XtmVsQCD NO
particleNet_XteVsQCD NO
'''
# Restore taggers that were decommisionned for Run-3
'''
particleNetLegacy_mass YES
particleNetLegacy_Xbb YES
particleNetLegacy_Xcc NO
particleNetLegacy_Xqq NO
particleNetLegacy_QCD YES
'''

# ParticleNet legacy jet tagger is already in 106Xv2 MINIAOD,
# add PartlceNet legacy mass regression and new combined tagger + mass regression

# for more info: ---> https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/jetsAK8_cff.py

FatJetObservablesMC = ["hadronFlavour","partonFlavour"]

SubJetObservables = ["btagDeepB", "eta", "mass", "phi", "pt", "rawFactor"]
SubJetObservablesMC = ["hadronFlavour","partonFlavour"]

defaultColToSave = ["entryIndex","luminosityBlock", "run","event", "sample_type", "sample_name", "period", "X_mass", "X_spin", "isData","PuppiMET_pt", "PuppiMET_phi", "nJet","DeepMETResolutionTune_pt", "DeepMETResolutionTune_phi","DeepMETResponseTune_pt", "DeepMETResponseTune_phi","PV_npvs"]

def getDefaultColumnsToSave(isData):
    colToSave = defaultColToSave.copy()
    if not isData:
        colToSave.extend(['Pileup_nTrueInt'])
    return colToSave


# rewrite this function
# add reco p4, create W (ntuple will contain only W)
def addAllVariables(dfw, syst_name, isData, trigger_class, lepton_legs, isSignal, global_params):
    print(f"Adding variables for {syst_name}")
    dfw.Apply(CommonBaseline.SelectRecoP4, syst_name, global_params["nano_version"])
    dfw.Apply(AnaBaseline.RecoHWWCandidateSelection)
    dfw.Apply(AnaBaseline.RecoHWWJetSelection)

    dfw.DefineAndAppend("lep1_pt", "lep1_p4.pt()")
    dfw.DefineAndAppend("lep1_eta", "lep1_p4.eta()")
    dfw.DefineAndAppend("lep1_phi", "lep1_p4.phi()")
    dfw.DefineAndAppend("lep1_mass", "lep1_p4.mass()")

    dfw.DefineAndAppend("lep2_pt", "lep2_p4.pt()")
    dfw.DefineAndAppend("lep2_eta", "lep2_p4.eta()")
    dfw.DefineAndAppend("lep2_phi", "lep2_p4.phi()")
    dfw.DefineAndAppend("lep2_mass", "lep2_p4.mass()")

    dfw.colToSave.extend(["lep1_type", "lep2_type"])

    kinematic_variables = ["pt", "eta", "phi", "mass"]
    for i in range(10):
        for var in kinematic_variables:
            name = f"Jet_{i}_{var}"
            dfw.DefineAndAppend(name, f"""  if (n_Jet_Sel > {i})
                                                return Jet_p4[Jet_sel][{i}].{var}();
                                            return 0.0; """)

    dfw.Define("H_to_VV", """GetGenHVVCandidate(event, genLeptons, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, true)""")
    dfw.DefineAndAppend(f"genHVV_E", f"H_to_VV.cand_p4.E()")
    dfw.DefineAndAppend(f"genHVV_px", f"H_to_VV.cand_p4.px()")
    dfw.DefineAndAppend(f"genHVV_py", f"H_to_VV.cand_p4.py()")
    dfw.DefineAndAppend(f"genHVV_pz", f"H_to_VV.cand_p4.pz()")

    dfw.DefineAndAppend(f"genV1_E", f"H_to_VV.legs[0].cand_p4.E()")
    dfw.DefineAndAppend(f"genV1_px", f"H_to_VV.legs[0].cand_p4.px()")
    dfw.DefineAndAppend(f"genV1_py", f"H_to_VV.legs[0].cand_p4.py()")
    dfw.DefineAndAppend(f"genV1_pz", f"H_to_VV.legs[0].cand_p4.pz()")

    dfw.DefineAndAppend(f"genV2_E", f"H_to_VV.legs[1].cand_p4.E()")
    dfw.DefineAndAppend(f"genV2_px", f"H_to_VV.legs[1].cand_p4.px()")
    dfw.DefineAndAppend(f"genV2_py", f"H_to_VV.legs[1].cand_p4.py()")
    dfw.DefineAndAppend(f"genV2_pz", f"H_to_VV.legs[1].cand_p4.pz()")

    dfw.DefineAndAppend(f"genV1prod1_E", f"H_to_VV.legs[0].leg_p4[0].E()")
    dfw.DefineAndAppend(f"genV1prod1_px", f"H_to_VV.legs[0].leg_p4[0].px()")
    dfw.DefineAndAppend(f"genV1prod1_py", f"H_to_VV.legs[0].leg_p4[0].py()")
    dfw.DefineAndAppend(f"genV1prod1_pz", f"H_to_VV.legs[0].leg_p4[0].pz()")

    dfw.DefineAndAppend(f"genV1prod2_E", f"H_to_VV.legs[0].leg_p4[1].E()")
    dfw.DefineAndAppend(f"genV1prod2_px", f"H_to_VV.legs[0].leg_p4[1].px()")
    dfw.DefineAndAppend(f"genV1prod2_py", f"H_to_VV.legs[0].leg_p4[1].py()")
    dfw.DefineAndAppend(f"genV1prod2_pz", f"H_to_VV.legs[0].leg_p4[1].pz()")

    dfw.DefineAndAppend(f"genV2prod1_E", f"H_to_VV.legs[1].leg_p4[0].E()")
    dfw.DefineAndAppend(f"genV2prod1_px", f"H_to_VV.legs[1].leg_p4[0].px()")
    dfw.DefineAndAppend(f"genV2prod1_py", f"H_to_VV.legs[1].leg_p4[0].py()")
    dfw.DefineAndAppend(f"genV2prod1_pz", f"H_to_VV.legs[1].leg_p4[0].pz()")

    dfw.DefineAndAppend(f"genV2prod2_E", f"H_to_VV.legs[1].leg_p4[1].E()")
    dfw.DefineAndAppend(f"genV2prod2_px", f"H_to_VV.legs[1].leg_p4[1].px()")
    dfw.DefineAndAppend(f"genV2prod2_py", f"H_to_VV.legs[1].leg_p4[1].py()")
    dfw.DefineAndAppend(f"genV2prod2_pz", f"H_to_VV.legs[1].leg_p4[1].pz()")

    dfw.Define("H_to_bb", """GetGenHBBCandidate(event, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, true)""")

    dfw.DefineAndAppend(f"genHbb_E", f"H_to_bb.cand_p4.E()")
    dfw.DefineAndAppend(f"genHbb_px", f"H_to_bb.cand_p4.px()")
    dfw.DefineAndAppend(f"genHbb_py", f"H_to_bb.cand_p4.py()")
    dfw.DefineAndAppend(f"genHbb_pz", f"H_to_bb.cand_p4.pz()")

    dfw.DefineAndAppend(f"genb1_E", f"H_to_bb.leg_p4[0].E()")
    dfw.DefineAndAppend(f"genb1_px", f"H_to_bb.leg_p4[0].px()")
    dfw.DefineAndAppend(f"genb1_py", f"H_to_bb.leg_p4[0].py()")
    dfw.DefineAndAppend(f"genb1_pz", f"H_to_bb.leg_p4[0].pz()")

    dfw.DefineAndAppend(f"genb2_E", f"H_to_bb.leg_p4[1].E()")
    dfw.DefineAndAppend(f"genb2_px", f"H_to_bb.leg_p4[1].px()")
    dfw.DefineAndAppend(f"genb2_py", f"H_to_bb.leg_p4[1].py()")
    dfw.DefineAndAppend(f"genb2_pz", f"H_to_bb.leg_p4[1].pz()")

