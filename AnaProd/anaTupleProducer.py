import copy
import datetime
import os
import sys
import ROOT
import shutil
import zlib

if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])

import Common.BaselineSelection as Baseline
import Common.Utilities as Utilities
import Common.ReportTools as ReportTools
import Common.triggerSel as Triggers
import Corrections.Corrections as Corrections
from Corrections.lumi import LumiFilter

#ROOT.EnableImplicitMT(1)
ROOT.EnableThreadSafety()

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

def SelectBTagShapeSF(df,weight_name):
    df = df.Define("weight_bTagShapeSF", weight_name)
    return df

def addAllVariables(dfw, syst_name, isData, trigger_class, mode, nLegs, isSignal):
    dfw.Apply(Baseline.SelectRecoP4, syst_name)
    dfw.Apply(Baseline.DefineMETCuts,80, ["MET", "DeepMETResolutionTune", "DeepMETResponseTune", "PuppiMET"])
    # qua va Select btagShapeWeight
    if mode == "HH":
        dfw.Apply(Baseline.RecoHttCandidateSelection, config["GLOBAL"])
        dfw.Apply(Baseline.RecoJetSelection)
        dfw.Apply(Baseline.ThirdLeptonVeto)
    elif mode == 'ttHH':
        dfw.Apply(Baseline.RecottHttCandidateSelection_ttHH)
        dfw.Apply(Baseline.RecoJetSelection_ttHH)

    '''
    if not isData:
        print(f"before PassGenAcceptance : {dfw.df.Count().GetValue()}")
        dfw.DefineAndAppend("n_GenJet", "GenJet_idx.size()")
        dfw.Apply(Baseline.PassGenAcceptance)
        print(f"after PassGenAcceptance : {dfw.df.Count().GetValue()}")
        dfw.Apply(Baseline.GenJetSelection)
        print(f"after GenJetSelection : {dfw.df.Count().GetValue()}")
        dfw.Apply(Baseline.GenJetHttOverlapRemoval)
        print(f"after GenJetHttOverlapRemoval : {dfw.df.Count().GetValue()}")
        dfw.Apply(Baseline.RequestOnlyResolvedGenJets)
        print(f"after RequestOnlyResolvedGenJets : {dfw.df.Count().GetValue()}")
    '''
    dfw.Apply(Baseline.DefineHbbCand)
    dfw.DefineAndAppend("Hbb_isValid" , "HbbCandidate.has_value()")
    dfw.Apply(Baseline.ExtraRecoJetSelection)
    dfw.Apply(Corrections.jet.getEnergyResolution)
    #dfw.Apply(Corrections.fatjet.getEnergyResolution)
    dfw.Apply(Corrections.btag.getWPid)
    jet_obs = []
    jet_obs.extend(JetObservables)
    dfw.Apply(Baseline.ApplyJetSelection)
    if not isData:
        dfw.Define(f"Jet_genJet_idx", f" FindMatching(Jet_p4,GenJet_p4,0.3)")
        jet_obs.extend(JetObservablesMC)


    dfw.DefineAndAppend(f"nBJets", f"Jet_p4[Jet_bCand].size()")
    if config["GLOBAL"]["storeExtraJets"]:
        dfw.DefineAndAppend(f"ExtraJet_pt", f"v_ops::pt(Jet_p4[ExtraJet_B1])")
        dfw.DefineAndAppend(f"ExtraJet_eta", f"v_ops::eta(Jet_p4[ExtraJet_B1])")
        dfw.DefineAndAppend(f"ExtraJet_phi", f"v_ops::phi(Jet_p4[ExtraJet_B1])")
        dfw.DefineAndAppend(f"ExtraJet_mass", f"v_ops::mass(Jet_p4[ExtraJet_B1])")
        dfw.DefineAndAppend(f"ExtraJet_ptRes", f"Jet_ptRes[ExtraJet_B1]")
        for jetVar in jet_obs:
            if(f"Jet_{jetVar}" not in dfw.df.GetColumnNames()): continue
            dfw.DefineAndAppend(f"ExtraJet_{jetVar}", f"Jet_{jetVar}[ExtraJet_B1]")
        dfw.DefineAndAppend(f"ExtraJet_HHbtag", f"Jet_HHBtagScore[ExtraJet_B1]")
    else:
        dfw.DefineAndAppend(f"nExtraJets", f"Jet_p4[ExtraJet_B1].size()")

    dfw.DefineAndAppend(f"met_pt_nano", f"static_cast<float>(MET_p4_nano.pt())")
    dfw.DefineAndAppend(f"met_phi_nano", f"static_cast<float>(MET_p4_nano.phi())")
    dfw.DefineAndAppend("met_pt", "static_cast<float>(MET_p4.pt())")
    dfw.DefineAndAppend("met_phi", "static_cast<float>(MET_p4.phi())")
    dfw.DefineAndAppend("metnomu_pt_nano", "static_cast<float>(GetMetNoMu(HttCandidate, MET_p4_nano).pt())")
    dfw.DefineAndAppend("metnomu_phi_nano", "static_cast<float>(GetMetNoMu(HttCandidate, MET_p4_nano).phi())")
    dfw.DefineAndAppend("metnomu_pt", "static_cast<float>(GetMetNoMu(HttCandidate, MET_p4).pt())")
    dfw.DefineAndAppend("metnomu_phi", "static_cast<float>(GetMetNoMu(HttCandidate, MET_p4).phi())")
    for var in ["covXX", "covXY", "covYY"]:
        dfw.DefineAndAppend(f"met_{var}", f"static_cast<float>(MET_{var})")

    if trigger_class is not None:
        hltBranches = dfw.Apply(trigger_class.ApplyTriggers, nLegs, isData, isSignal)
        dfw.colToSave.extend(hltBranches)
    dfw.Define(f"Tau_recoJetMatchIdx", f"FindMatching(Tau_p4, Jet_p4, 0.5)")
    dfw.Define(f"Muon_recoJetMatchIdx", f"FindMatching(Muon_p4, Jet_p4, 0.5)")
    dfw.Define( f"Electron_recoJetMatchIdx", f"FindMatching(Electron_p4, Jet_p4, 0.5)")
    dfw.DefineAndAppend("channelId","static_cast<int>(HttCandidate.channel())")
    if mode == "HH":
        channel_to_select = " || ".join(f"HttCandidate.channel()==Channel::{ch}" for ch in config["GLOBAL"]["channelSelection"])
        dfw.Filter(channel_to_select, "select channels")
        #print(f"after {channel_to_select}: {dfw.df.Count().GetValue()}")
    fatjet_obs = []
    fatjet_obs.extend(FatJetObservables)
    if not isData:
        dfw.Define(f"FatJet_genJet_idx", f" FindMatching(FatJet_p4[FatJet_bbCand],GenJetAK8_p4,0.3)")
        fatjet_obs.extend(JetObservablesMC)
        if isSignal:
            dfw.DefineAndAppend("genchannelId","static_cast<int>(genHttCandidate->channel())")
        #dfw.df.Display({"channelId","genchannelId"}).Print()
    dfw.DefineAndAppend(f"SelectedFatJet_pt", f"v_ops::pt(FatJet_p4[FatJet_bbCand])")
    dfw.DefineAndAppend(f"SelectedFatJet_eta", f"v_ops::eta(FatJet_p4[FatJet_bbCand])")
    dfw.DefineAndAppend(f"SelectedFatJet_phi", f"v_ops::phi(FatJet_p4[FatJet_bbCand])")
    dfw.DefineAndAppend(f"SelectedFatJet_mass", f"v_ops::mass(FatJet_p4[FatJet_bbCand])")

    for fatjetVar in fatjet_obs:
        if(f"FatJet_{fatjetVar}" not in dfw.df.GetColumnNames()): continue
        dfw.DefineAndAppend(f"SelectedFatJet_{fatjetVar}", f"FatJet_{fatjetVar}[FatJet_bbCand]")
    subjet_obs = []
    subjet_obs.extend(SubJetObservables)
    if not isData:
        dfw.Define(f"SubJet1_genJet_idx", f" FindMatching(SubJet_p4[FatJet_subJetIdx1],SubGenJetAK8_p4,0.3)")
        dfw.Define(f"SubJet2_genJet_idx", f" FindMatching(SubJet_p4[FatJet_subJetIdx2],SubGenJetAK8_p4,0.3)")
        fatjet_obs.extend(SubJetObservablesMC)
    for subJetIdx in [1,2]:
        dfw.Define(f"SelectedFatJet_subJetIdx{subJetIdx}", f"FatJet_subJetIdx{subJetIdx}[FatJet_bbCand]")
        dfw.Define(f"FatJet_SubJet{subJetIdx}_isValid", f" FatJet_subJetIdx{subJetIdx} >=0 && FatJet_subJetIdx{subJetIdx} < nSubJet")
        dfw.DefineAndAppend(f"SelectedFatJet_SubJet{subJetIdx}_isValid", f"FatJet_SubJet{subJetIdx}_isValid[FatJet_bbCand]")
        for subJetVar in subjet_obs:
            dfw.DefineAndAppend(f"SelectedFatJet_SubJet{subJetIdx}_{subJetVar}", f"""
                                RVecF subjet_var(SelectedFatJet_pt.size(), 0.f);
                                for(size_t fj_idx = 0; fj_idx<SelectedFatJet_pt.size(); fj_idx++) {{
                                    auto sj_idx = SelectedFatJet_subJetIdx{subJetIdx}.at(fj_idx);
                                    if(sj_idx >= 0 && sj_idx < SubJet_{subJetVar}.size()){{
                                        subjet_var[fj_idx] = SubJet_{subJetVar}.at(sj_idx);
                                    }}
                                }}
                                return subjet_var;
                                """)

    n_legs = 2 if mode == "HH" else 4
    for leg_idx in range(n_legs):
        def LegVar(var_name, var_expr, var_type=None, var_cond=None, check_leg_type=True, default=0):
            cond = var_cond
            if check_leg_type:
                type_cond = f"HttCandidate.leg_type[{leg_idx}] != Leg::none"
                cond = f"{type_cond} && ({cond})" if cond else type_cond
            define_expr = f'static_cast<{var_type}>({var_expr})' if var_type else var_expr
            if cond:
                define_expr = f'{cond} ? ({define_expr}) : {default}'
            dfw.DefineAndAppend( f"tau{leg_idx+1}_{var_name}", define_expr)

        LegVar('legType', f"HttCandidate.leg_type[{leg_idx}]", var_type='int', check_leg_type=False)
        for var in [ 'pt', 'eta', 'phi', 'mass' ]:
            LegVar(var, f'HttCandidate.leg_p4[{leg_idx}].{var}()', var_type='float', default='-1.f')
        LegVar('charge', f'HttCandidate.leg_charge[{leg_idx}]', var_type='int')

        dfw.Define(f"tau{leg_idx+1}_recoJetMatchIdx", f"""HttCandidate.leg_type[{leg_idx}] != Leg::none
                                                          ? FindMatching(HttCandidate.leg_p4[{leg_idx}], Jet_p4, 0.3)
                                                          : -1""")
        LegVar('iso', f"HttCandidate.leg_rawIso.at({leg_idx})")
        for deepTauScore in deepTauScores:
            LegVar(deepTauScore, f"Tau_{deepTauScore}.at(HttCandidate.leg_index[{leg_idx}])",
                   var_cond=f"HttCandidate.leg_type[{leg_idx}] == Leg::tau", default='-1.f')
        for muon_obs in Muon_observables:
            LegVar(muon_obs, f"{muon_obs}.at(HttCandidate.leg_index[{leg_idx}])",
                   var_cond=f"HttCandidate.leg_type[{leg_idx}] == Leg::mu", default='-1')
        for ele_obs in Electron_observables:
            LegVar(ele_obs, f"{ele_obs}.at(HttCandidate.leg_index[{leg_idx}])",
                   var_cond=f"HttCandidate.leg_type[{leg_idx}] == Leg::e", default='-1')
        if not isData:
            dfw.Define(f"tau{leg_idx+1}_genMatchIdx",
                       f"HttCandidate.leg_type[{leg_idx}] != Leg::none ? HttCandidate.leg_genMatchIdx[{leg_idx}] : -1")
            LegVar('gen_kind', f'genLeptons.at(tau{leg_idx+1}_genMatchIdx).kind()',
                   var_type='int', var_cond=f"tau{leg_idx+1}_genMatchIdx>=0",
                   default='static_cast<int>(GenLeptonMatch::NoMatch)')
            for var in [ 'pt', 'eta', 'phi', 'mass' ]:
                LegVar(f'gen_vis_{var}', f'genLeptons.at(tau{leg_idx+1}_genMatchIdx).visibleP4().{var}()',
                       var_type='float', var_cond=f"tau{leg_idx+1}_genMatchIdx>=0", default='-1.f')
            LegVar('gen_nChHad', f'genLeptons.at(tau{leg_idx+1}_genMatchIdx).nChargedHadrons()',
                   var_type='int', var_cond=f"tau{leg_idx+1}_genMatchIdx>=0", default='-1')
            LegVar('gen_nNeutHad', f'genLeptons.at(tau{leg_idx+1}_genMatchIdx).nNeutralHadrons()',
                   var_type='int', var_cond=f"tau{leg_idx+1}_genMatchIdx>=0", default='-1')
            LegVar('gen_charge', f'genLeptons.at(tau{leg_idx+1}_genMatchIdx).charge()',
                   var_type='int', var_cond=f"tau{leg_idx+1}_genMatchIdx>=0", default='-10')
            LegVar('seedingJet_partonFlavour', f'Jet_partonFlavour.at(tau{leg_idx+1}_recoJetMatchIdx)',
                   var_type='int', var_cond=f"tau{leg_idx+1}_recoJetMatchIdx>=0", default='-10')
            LegVar('seedingJet_hadronFlavour', f'Jet_hadronFlavour.at(tau{leg_idx+1}_recoJetMatchIdx)',
                   var_type='int', var_cond=f"tau{leg_idx+1}_recoJetMatchIdx>=0", default='-10')

        for var in [ 'pt', 'eta', 'phi', 'mass' ]:
            LegVar(f'seedingJet_{var}', f"Jet_p4.at(tau{leg_idx+1}_recoJetMatchIdx).{var}()",
                   var_type='float', var_cond=f"tau{leg_idx+1}_recoJetMatchIdx>=0", default='-1.f')
        if mode == "HH":
            dfw.Define(f"b{leg_idx+1}_idx", f"Hbb_isValid ? HbbCandidate->leg_index[{leg_idx}] : -100")
            dfw.DefineAndAppend(f"b{leg_idx+1}_ptRes",f"Hbb_isValid ? static_cast<float>(Jet_ptRes.at(HbbCandidate->leg_index[{leg_idx}])) : -100.f")
            dfw.DefineAndAppend(f"b{leg_idx+1}_pt", f"Hbb_isValid ? static_cast<float>(HbbCandidate->leg_p4[{leg_idx}].Pt()) : -100.f")
            dfw.DefineAndAppend(f"b{leg_idx+1}_pt_raw", f"Hbb_isValid ? static_cast<float>(Jet_pt.at(HbbCandidate->leg_index[{leg_idx}])) : - 100.f")
            dfw.DefineAndAppend(f"b{leg_idx+1}_eta", f"Hbb_isValid ? static_cast<float>(HbbCandidate->leg_p4[{leg_idx}].Eta()) : -100.f")
            dfw.DefineAndAppend(f"b{leg_idx+1}_phi", f"Hbb_isValid ? static_cast<float>(HbbCandidate->leg_p4[{leg_idx}].Phi()) : -100.f")
            dfw.DefineAndAppend(f"b{leg_idx+1}_mass", f"Hbb_isValid ? static_cast<float>(HbbCandidate->leg_p4[{leg_idx}].M()) : -100.f")
            if not isData:
                dfw.Define(f"b{leg_idx+1}_genJet_idx", f" Hbb_isValid ?  Jet_genJet_idx.at(HbbCandidate->leg_index[{leg_idx}]) : -100")
                for var in [ 'pt', 'eta', 'phi', 'mass' ]:
                    dfw.DefineAndAppend(f"b{leg_idx+1}_genJet_{var}", f"Hbb_isValid && b{leg_idx+1}_genJet_idx>=0 ? static_cast<float>(GenJet_p4.at(b{leg_idx+1}_genJet_idx).{var}()) : -100.f")
            for jetVar in jet_obs:
                if(f"Jet_{jetVar}" not in dfw.df.GetColumnNames()): continue
                dfw.DefineAndAppend(f"b{leg_idx+1}_{jetVar}", f"Hbb_isValid ? Jet_{jetVar}.at(HbbCandidate->leg_index[{leg_idx}]) : -100")
            dfw.DefineAndAppend(f"b{leg_idx+1}_HHbtag", f"Hbb_isValid ?  static_cast<float>(Jet_HHBtagScore.at(HbbCandidate->leg_index[{leg_idx}])) : -100.f")

def createAnatuple(inFile, treeName, outDir, config, sample_name, anaCache, snapshotOptions,range, evtIds,
                   store_noncentral, compute_unc_variations, uncertainties, print_cutflow, mode):
    start_time = datetime.datetime.now()
    compression_settings = snapshotOptions.fCompressionAlgorithm * 100 + snapshotOptions.fCompressionLevel
    period = config["GLOBAL"]["era"]
    mass = -1 if 'mass' not in config[sample_name] else config[sample_name]['mass']
    spin = -100 if 'spin' not in config[sample_name] else config[sample_name]['spin']
    isHH = True if mass > 0 else False
    isData = True if config[sample_name]['sampleType'] == 'data' else False
    loadTF = mode == "HH"
    loadHHBtag = mode == "HH"
    nLegs = 4 if mode == "ttHH" else 2
    Baseline.Initialize(loadTF, loadHHBtag)
    Corrections.Initialize(config=config['GLOBAL'],isData=isData)
    triggerFile = config['GLOBAL'].get('triggerFile')
    #print(triggerFile)
    if triggerFile is not None:
        triggerFile = os.path.join(os.environ['ANALYSIS_PATH'], triggerFile)
        trigger_class = Triggers.Triggers(triggerFile)
    else:
        trigger_class = None
    inFiles = Utilities.ListToVector(inFile.split(','))
    df = ROOT.RDataFrame(treeName, inFiles)
    if range is not None:
        df = df.Range(range)
    #print(f"at the beginning: {df.Count().GetValue()}")
    if len(evtIds) > 0:
        df = df.Filter(f"static const std::set<ULong64_t> evts = {{ {evtIds} }}; return evts.count(event) > 0;")


    if isData and 'lumiFile' in config['GLOBAL']:
        lumiFilter = LumiFilter(config['GLOBAL']['lumiFile'])
        df = lumiFilter.filter(df)
    df = Baseline.applyMETFlags(df, config["GLOBAL"]["MET_flags"])
    df = df.Define("sample_type", f"static_cast<int>(SampleType::{config[sample_name]['sampleType']})")
    df = df.Define("sample_name", f"{zlib.crc32(sample_name.encode())}")
    isSignal = config[sample_name]['sampleType'] in config['GLOBAL']['signal_types']
    print("isSignal? ", isSignal)
    df = df.Define("period", f"static_cast<int>(Period::{period})")
    df = df.Define("X_mass", f"static_cast<int>({mass})")
    df = df.Define("X_spin", f"static_cast<int>({spin})")
    df = df.Define("entryIndex", "static_cast<int>(rdfentry_)")
    is_data = 'true' if isData else 'false'
    df = df.Define("isData", is_data)

    df = Baseline.CreateRecoP4(df)
    df = Baseline.DefineGenObjects(df, isData=isData, isHH=isHH)
    if isData:
        syst_dict = { 'nano' : 'Central' }
    else:
        df, syst_dict = Corrections.applyScaleUncertainties(df)
    df_empty = df
    snaps = []
    reports = []
    outfilesNames = []
    k=0
    if not isData:
        defaultColToSave.extend(['Pileup_nTrueInt'])
    for syst_name, source_name in syst_dict.items():
        if source_name not in uncertainties and "all" not in uncertainties: continue
        is_central = syst_name in [ 'Central', 'nano' ]
        if not is_central and not compute_unc_variations: continue
        suffix = '' if is_central else f'_{syst_name}'
        if len(suffix) and not store_noncentral: continue
        dfw = Utilities.DataFrameWrapper(df_empty,defaultColToSave)
        #print(syst_name)

        addAllVariables(dfw, syst_name, isData, trigger_class, mode, nLegs, isSignal)

        dfw.DefineAndAppend("weight_L1PreFiring_Central","L1PreFiringWeight_Nom")
        dfw.DefineAndAppend("weight_L1PreFiring_ECAL_Central","L1PreFiringWeight_ECAL_Nom")
        dfw.DefineAndAppend("weight_L1PreFiring_Muon_Central","L1PreFiringWeight_Muon_Nom")
        if is_central and compute_unc_variations:
            dfw.DefineAndAppend("weight_L1PreFiringDown_rel","L1PreFiringWeight_Dn/L1PreFiringWeight_Nom")
            dfw.DefineAndAppend("weight_L1PreFiringUp_rel","L1PreFiringWeight_Up/L1PreFiringWeight_Nom")
            dfw.DefineAndAppend("weight_L1PreFiring_ECALDown_rel","L1PreFiringWeight_ECAL_Dn/L1PreFiringWeight_ECAL_Nom")
            dfw.DefineAndAppend("weight_L1PreFiring_ECALUp_rel","L1PreFiringWeight_ECAL_Up/L1PreFiringWeight_ECAL_Nom")
            dfw.DefineAndAppend("weight_L1PreFiring_Muon_StatDown_rel", "L1PreFiringWeight_Muon_StatDn/L1PreFiringWeight_Muon_Nom")
            dfw.DefineAndAppend("weight_L1PreFiring_Muon_StatUp_rel", "L1PreFiringWeight_Muon_StatUp/L1PreFiringWeight_Muon_Nom")
            dfw.DefineAndAppend("weight_L1PreFiring_Muon_SystDown_rel", "L1PreFiringWeight_Muon_SystDn/L1PreFiringWeight_Muon_Nom")
            dfw.DefineAndAppend("weight_L1PreFiring_Muon_SystUp_rel", "L1PreFiringWeight_Muon_SystUp/L1PreFiringWeight_Muon_Nom")
        if not isData:
            weight_branches = dfw.Apply(Corrections.getNormalisationCorrections, config, sample_name, nLegs,
                                        return_variations=is_central and compute_unc_variations, isCentral=is_central,
                                        ana_cache=anaCache)
            weight_branches.extend(dfw.Apply(Corrections.trg.getTrgSF, trigger_class.trigger_dict.keys(), nLegs,
                                             is_central and compute_unc_variations, is_central))
            SF_branches_core,SF_weight_jes=dfw.Apply(Corrections.jet.getBtagShapeSFs, syst_name, is_central, compute_unc_variations)
            syst_name_selected = 'Central' if SF_weight_jes=="" else syst_name
            SF_branches_core.remove(f'weight_bTagShapeSF_Central')
            weight_name = f'weight_bTagShapeSF_{syst_name_selected}'
            dfw.df = SelectBTagShapeSF(dfw.df, weight_name)
            if is_central: weight_branches.extend(SF_branches_core)
            weight_branches.extend([f'weight_bTagShapeSF'])
            puIDbranches = ["weight_Jet_PUJetID_Central_tmp", "weight_Jet_PUJetID_effUp_rel_tmp", "weight_Jet_PUJetID_effDown_rel_tmp"]
            for puIDbranch in puIDbranches:
                if puIDbranch in dfw.df.GetColumnNames():
                    new_branch_name= puIDbranch.strip("_tmp")
                    dfw.Define(f"""ExtraJet_{new_branch_name}""", f"{puIDbranch}[ExtraJet_B1]")
                    if config["GLOBAL"]["storeExtraJets"]:
                        dfw.colToSave.append(f"""ExtraJet_{new_branch_name}""")
                    for bjet_idx in [1,2]:
                        dfw.DefineAndAppend(f"{new_branch_name}_b{bjet_idx}", f"Hbb_isValid ? {puIDbranch}[b{bjet_idx}_idx] : -100.f")
                if puIDbranch in weight_branches: weight_branches.remove(puIDbranch)
            dfw.colToSave.extend(weight_branches)
        reports.append(dfw.df.Report())
        varToSave = Utilities.ListToVector(dfw.colToSave)
        outfile_prefix = inFiles[0].split('/')[-1]
        outfile_prefix = outfile_prefix.split('.')[0]
        outFileName = os.path.join(outDir, f"{outfile_prefix}{suffix}.root")
        outfilesNames.append(outFileName)
        if os.path.exists(outFileName):
            os.remove(outFileName)
        snaps.append(dfw.df.Snapshot(f"Events", outFileName, varToSave, snapshotOptions))
    if snapshotOptions.fLazy == True:
        ROOT.RDF.RunGraphs(snaps)
    hist_time = ROOT.TH1D(f"time", f"time", 1, 0, 1)
    end_time = datetime.datetime.now()
    hist_time.SetBinContent(1, (end_time - start_time).total_seconds())
    for index,fileName in enumerate(outfilesNames):
        outputRootFile= ROOT.TFile(fileName, "UPDATE", "", compression_settings)
        rep = ReportTools.SaveReport(reports[index].GetValue(), reoprtName=f"Report")
        outputRootFile.WriteTObject(rep, f"Report", "Overwrite")
        if index==0:
            outputRootFile.WriteTObject(hist_time, f"runtime", "Overwrite")
        outputRootFile.Close()
        # if print_cutflow:
        #     report.Print()
    #print(f"number of loops is {df_empty.GetNRuns()}")


if __name__ == "__main__":
    import argparse
    import os
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True, type=str)
    parser.add_argument('--inFile', required=True, type=str)
    parser.add_argument('--treeName', required=False, type=str, default="Events")
    parser.add_argument('--outDir', required=True, type=str)
    parser.add_argument('--sample', required=True, type=str)
    parser.add_argument('--anaCache', required=True, type=str)
    parser.add_argument('--compressionLevel', type=int, default=4)
    parser.add_argument('--compressionAlgo', type=str, default="ZLIB")
    parser.add_argument('--nEvents', type=int, default=None)
    parser.add_argument('--evtIds', type=str, default='')
    parser.add_argument('--store-noncentral', action="store_true", help="Store ES variations.")
    parser.add_argument('--compute_unc_variations', type=bool, default=False)
    parser.add_argument('--print-cutflow', type=bool, default=False)
    parser.add_argument('--customisations', type=str, default="")
    parser.add_argument('--uncertainties', type=str, default="all")
    parser.add_argument('--mode', type=str, default="HH")
    parser.add_argument('--particleFile', type=str,
                        default=f"{os.environ['ANALYSIS_PATH']}/config/pdg_name_type_charge.txt")
    args = parser.parse_args()

    ROOT.gROOT.ProcessLine(".include "+ os.environ['ANALYSIS_PATH'])
    ROOT.gROOT.ProcessLine('#include "include/GenTools.h"')
    ROOT.gInterpreter.ProcessLine(f"ParticleDB::Initialize(\"{args.particleFile}\");")
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    if len(args.customisations)>0:
        Utilities.ApplyConfigCustomisations(config['GLOBAL'], args.customisations)
    with open(args.anaCache, 'r') as f:
        anaCache = yaml.safe_load(f)

    if os.path.isdir(args.outDir):
        shutil.rmtree(args.outDir)
    os.makedirs(args.outDir, exist_ok=True)
    snapshotOptions = ROOT.RDF.RSnapshotOptions()
    snapshotOptions.fOverwriteIfExists=False
    snapshotOptions.fLazy = True
    snapshotOptions.fMode="RECREATE"
    snapshotOptions.fCompressionAlgorithm = getattr(ROOT.ROOT, 'k' + args.compressionAlgo)
    snapshotOptions.fCompressionLevel = args.compressionLevel
    createAnatuple(args.inFile,args.treeName, args.outDir, config, args.sample, anaCache, snapshotOptions, args.nEvents,
                   args.evtIds, args.store_noncentral, args.compute_unc_variations, args.uncertainties.split(","), args.print_cutflow, args.mode)
