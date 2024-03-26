import ROOT
import os
from .Utilities import *

initialized = False

#ana_reco_object_collections = [ "Electron", "Muon", "Tau", "Jet", "FatJet", "boostedTau", "MET", "PuppiMET", "DeepMETResponseTune", "DeepMETResolutionTune","SubJet"]
ana_reco_object_collections = [ "Electron", "Muon", "Tau", "Jet", "FatJet", "MET", "PuppiMET", "DeepMETResponseTune", "DeepMETResolutionTune","SubJet"]
deepTauVersions = {"2p1":"2017", "2p5":"2018"}

def Initialize(loadTF=False, loadHHBtag=False):
    global initialized
    if not initialized:
        headers_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gROOT.ProcessLine(".include "+ os.environ['ANALYSIS_PATH'])
        header_path_RootExt ="include/RootExt.h"
        header_path_GenLepton ="include/GenLepton.h"
        header_path_Gen ="include/BaselineGenSelection.h"
        header_path_Reco ="include/BaselineRecoSelection.h"
        header_path_HHbTag ="include/HHbTagScores.h"
        ROOT.gInterpreter.Declare(f'#include "{header_path_RootExt}"')
        ROOT.gInterpreter.Declare(f'#include "{header_path_GenLepton}"')
        ROOT.gInterpreter.Declare(f'#include "{header_path_Gen}"')
        ROOT.gInterpreter.Declare(f'#include "{header_path_Reco}"')
        for wpcl in [WorkingPointsTauVSe,WorkingPointsTauVSmu,WorkingPointsTauVSjet,WorkingPointsbTag, WorkingPointsMuonID]:
            ROOT.gInterpreter.Declare(f'{generate_enum_class(wpcl)}')
        if(loadTF):
            import RunKit.includeCMSSWlibs as IncludeLibs
            IncludeLibs.includeLibTool("tensorflow")
        if(loadHHBtag):
            ROOT.gInterpreter.Declare(f'#include "{header_path_HHbTag}"')
            ROOT.gROOT.ProcessLine(f'HHBtagWrapper::Initialize("{os.environ["CMSSW_BASE"]}/src/HHTools/HHbtag/models/", 2)')
        initialized = True

#leg_names = [ "Electron", "Muon", "Tau", "boostedTau" ]
leg_names = [ "Electron", "Muon", "Tau" ]

channels = [ 'muMu', 'eMu', 'eE', 'muTau', 'eTau', 'tauTau' ] # in order of importance during the channel selection

def getChannelLegs(channel):
    ch_str = channel.lower()
    legs = []
    while len(ch_str) > 0:
        name_idx = None
        obj_name = None
        for idx, obj in enumerate(['e', 'mu', 'tau']):
            if ch_str.startswith(obj):
                name_idx = idx
                obj_name = obj
                break
        if name_idx is None:
            raise RuntimeError(f"Invalid channel name {channel}")
        legs.append(leg_names[name_idx])
        ch_str = ch_str[len(obj_name):]
    return legs


def applyMETFlags(df, MET_flags):
    MET_flags_string = ' && '.join(MET_flags)
    return df.Filter(MET_flags_string, "MET filters")

def DefineGenObjects(df, isData=False, isHH=False, Hbb_AK4mass_mpv=125., p4_suffix='nano'):
    if isData:
        df = df.Define("genLeptons", "std::vector<reco_tau::gen_truth::GenLepton>()")
    else:
        df = df.Define("GenPart_daughters", "GetDaughters(GenPart_genPartIdxMother)")
        df = df.Define("genLeptons", """reco_tau::gen_truth::GenLepton::fromNanoAOD(GenPart_pt, GenPart_eta,
                                        GenPart_phi, GenPart_mass, GenPart_genPartIdxMother, GenPart_pdgId,
                                        GenPart_statusFlags, event)""")

    for lep in ["Electron", "Muon", "Tau"]:
        df = df.Define(f"{lep}_genMatchIdx",  f"MatchGenLepton({lep}_p4_{p4_suffix}, genLeptons, 0.2)")
        df = df.Define(f"{lep}_genMatch",  f"GetGenLeptonMatch({lep}_genMatchIdx, genLeptons)")
    if isData:
        return df

    if isHH:
        df = df.Define("genHttCandidate", """GetGenHTTCandidate(event, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, false)""")
        df = df.Define("genHbbIdx", """GetGenHBBIndex(event, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags)""")
        df = df.Define("genHbb_isBoosted", "GenPart_pt[genHbbIdx]>550")
    for var in ["GenJet", "GenJetAK8", "SubGenJetAK8"]:
        df = df.Define(f"{var}_idx", f"CreateIndexes({var}_pt.size())")
        df = df.Define(f"{var}_p4", f"GetP4({var}_pt,{var}_eta,{var}_phi,{var}_mass, {var}_idx)")

    df = df.Define("GenJet_b_PF", "abs(GenJet_partonFlavour)==5")
    df = df.Define("GenJetAK8_b_PF", "abs(GenJetAK8_partonFlavour)==5")
    df = df.Define("GenJet_Hbb",f"FindTwoJetsClosestToMPV({Hbb_AK4mass_mpv}, GenJet_p4, GenJet_b_PF)")
    df = df.Define("GenJetAK8_Hbb", "FindGenJetAK8(GenJetAK8_mass, GenJetAK8_b_PF)")

    return df

def PassGenAcceptance(df):
    df = df.Filter("genHttCandidate.get() != nullptr", "genHttCandidate present")
    return df.Filter("PassGenAcceptance(*genHttCandidate)", "genHttCandidate Acceptance")

def GenJetSelection(df):
    df = df.Define("GenJet_B1","GenJet_pt > 20 && abs(GenJet_eta) < 2.5 && GenJet_Hbb")
    df = df.Define("GenJetAK8_B1","GenJetAK8_pt > 170 && abs(GenJetAK8_eta) < 2.5 && GenJetAK8_Hbb")
    return df.Filter("GenJet_idx[GenJet_B1].size()==2 || (GenJetAK8_idx[GenJetAK8_B1].size()==1 && genHbb_isBoosted)", "(One)Two b-parton (Fat)jets at least")

def GenJetHttOverlapRemoval(df):
    for var in ["GenJet", "GenJetAK8"]:
        df = df.Define(f"{var}_B2", f"RemoveOverlaps({var}_p4, {var}_B1,{{{{genHttCandidate->leg_p4[0], genHttCandidate->leg_p4[1]}},}}, 2, 0.5)" )
    return df.Filter("GenJet_idx[GenJet_B2].size()==2 || (GenJetAK8_idx[GenJetAK8_B2].size()==1 && genHbb_isBoosted)", "No overlap between genJets and genHttCandidates")

def RequestOnlyResolvedGenJets(df):
    return df.Filter("GenJet_idx[GenJet_B2].size()==2", "Resolved topology")


def SelectRecoP4(df, syst_name='nano'):
    for obj in ana_reco_object_collections:
        df = df.Define(f"{obj}_p4", f"{obj}_p4_{syst_name}")
    return df

def CreateRecoP4(df, suffix='nano'):
    if len(suffix) > 0:
        suffix = "_" + suffix
    if("TrigObj_pt" in df.GetColumnNames()):
        df = df.Define(f"TrigObj_idx", f"CreateIndexes(TrigObj_pt.size())")
        if "TrigObj_mass" not in df.GetColumnNames():
            df = df.Define("TrigObj_mass", "RVecF(TrigObj_pt.size(), 0.f)")
        df = df.Define(f"TrigObj_p4", f"GetP4(TrigObj_pt,TrigObj_eta,TrigObj_phi, TrigObj_mass, TrigObj_idx)")
    for obj in ana_reco_object_collections:
        if "MET" in obj:
            df = df.Define(f"{obj}_p4{suffix}", f"LorentzVectorM({obj}_pt, 0., {obj}_phi, 0.)")
        else:
            df = df.Define(f"{obj}_idx", f"CreateIndexes({obj}_pt.size())")
            df = df.Define(f"{obj}_p4{suffix}", f"GetP4({obj}_pt, {obj}_eta, {obj}_phi, {obj}_mass, {obj}_idx)")
    return df

def DefineMETCuts(df, met_thr, met_collections):
    cut = ' || '.join([f'{v}_pt > {met_thr}' for v in met_collections ])
    df = df.Filter(cut, "MET cuts")
    return df


def RecoHttCandidateSelection(df, config):

    df = df.Define("Electron_B0", f"""
        v_ops::pt(Electron_p4) > 10 && abs(v_ops::eta(Electron_p4)) < 2.5 && abs(Electron_dz) < 0.2 && abs(Electron_dxy) < 0.045  """) # && (Electron_mvaIso_WP90 || (Electron_mvaNoIso_WP90 && Electron_pfRelIso03_all < 0.5))

    df = df.Define("Muon_B0", f"""
        v_ops::pt(Muon_p4) > 15 && abs(v_ops::eta(Muon_p4)) < 2.4 && abs(Muon_dz) < 0.2 && abs(Muon_dxy) < 0.045
    """) # && ( ((Muon_tightId || Muon_mediumId) && Muon_pfRelIso04_all < 0.5) || (Muon_highPtId && Muon_tkRelIso < 0.5) )

    eta_cut = 2.3 if config["deepTauVersion"] == '2p1' else 2.5
    df = df.Define("Tau_B0", f"""
        v_ops::pt(Tau_p4) > 20 && abs(v_ops::eta(Tau_p4)) < {eta_cut} && abs(Tau_dz) < 0.2 && Tau_decayMode != 5 && Tau_decayMode != 6 && ( Tau_idDeepTau{deepTauVersions[config["deepTauVersion"]]}v{config["deepTauVersion"]}VSjet >= {WorkingPointsTauVSjet.VVVLoose.value} )
    """)
    df = df.Define("Electron_iso", "Electron_pfRelIso03_all") \
           .Define("Muon_iso", "Muon_pfRelIso04_all") \
           .Define("Tau_iso", f"""-Tau_rawDeepTau{deepTauVersions[config["deepTauVersion"]]}v{config["deepTauVersion"]}VSjet""")


    df = df.Define("Electron_B2_eTau_1", f"Electron_B0 && Electron_mvaIso_WP80 ")
    #  Electron_mvaNoIso_WP80 && Electron_pfRelIso03_all < 0.3


    df = df.Define("Muon_B2_muTau_1", f"""
        Muon_B0 &&  ( (Muon_tightId && Muon_pfRelIso04_all < 0.15) || (Muon_highPtId && Muon_tkRelIso < 0.15) )
    """)

    for ch in [ "eTau", "muTau", "tauTau" ]:
        cut_str = f'''Tau_B0
            && (Tau_idDeepTau{deepTauVersions[config["deepTauVersion"]]}v{config["deepTauVersion"]}VSe >= {getattr(WorkingPointsTauVSe, config["deepTauWPs"][ch]["VSe"]).value})
            && (Tau_idDeepTau{deepTauVersions[config["deepTauVersion"]]}v{config["deepTauVersion"]}VSmu >= {getattr(WorkingPointsTauVSmu, config["deepTauWPs"][ch]["VSmu"]).value})'''
        df = df.Define(f'Tau_B2_{ch}_2', cut_str)
        if ch == 'tauTau':
            cut_str += f' && (Tau_idDeepTau{deepTauVersions[config["deepTauVersion"]]}v{config["deepTauVersion"]}VSjet >= {getattr(WorkingPointsTauVSjet, config["deepTauWPs"]["tauTau"]["VSjet"]).value})'
            df = df.Define(f'Tau_B2_{ch}_1', cut_str)


    df = df.Define("Muon_B2_muMu_1", f"""
        Muon_B0 && ( (Muon_tightId && Muon_pfRelIso04_all < 0.15) || (Muon_highPtId && Muon_tkRelIso < 0.15) )
    """)
    df = df.Define("Muon_B2_muMu_2", f"""
        Muon_B0 && ( (Muon_tightId && Muon_pfRelIso04_all < 0.3) || (Muon_highPtId && Muon_tkRelIso < 0.3) )
    """)

    df = df.Define("Electron_B2_eMu_1",f"Electron_B0 && Electron_mvaIso_WP80 ")
    #  Electron_mvaNoIso_WP80 && Electron_pfRelIso03_all < 0.3
    df = df.Define("Muon_B2_eMu_2", f"""
        Muon_B0 && ( (Muon_tightId && Muon_pfRelIso04_all < 0.3) || (Muon_highPtId && Muon_tkRelIso < 0.3) )
    """)

    df = df.Define("Electron_B2_eE_1",f"Electron_B0 && Electron_mvaIso_WP80 ")
    #  Electron_mvaNoIso_WP80 && Electron_pfRelIso03_all < 0.3
    df = df.Define("Electron_B2_eE_2", f""" Electron_B0 && Electron_mvaNoIso_WP80 && Electron_pfRelIso03_all < 0.3 """)

    cand_columns = []
    for ch in channels:
        leg1, leg2 = getChannelLegs(ch)
        cand_column = f"HttCandidates_{ch}"
        df = df.Define(cand_column, f"""
            GetHTTCandidates<2>(Channel::{ch}, 0.5, {leg1}_B2_{ch}_1, {leg1}_p4, {leg1}_iso, {leg1}_charge, {leg1}_genMatchIdx,
                                                 {leg2}_B2_{ch}_2, {leg2}_p4, {leg2}_iso, {leg2}_charge, {leg2}_genMatchIdx)
        """)
        cand_columns.append(cand_column)
    cand_filters = [ f'{c}.size() > 0' for c in cand_columns ]
    stringfilter = " || ".join(cand_filters)
    df = df.Filter(" || ".join(cand_filters), "Reco Baseline 2")
    cand_list_str = ', '.join([ '&' + c for c in cand_columns])
    return df.Define('HttCandidate', f'GetBestHTTCandidate<2>({{ {cand_list_str} }}, event)')

def ThirdLeptonVeto(df):
    df = df.Define("Electron_vetoSel",
                   f"""v_ops::pt(Electron_p4) > 10 && abs(v_ops::eta(Electron_p4)) < 2.5 && abs(Electron_dz) < 0.2 && abs(Electron_dxy) < 0.045 && ( Electron_mvaIso_WP90 == true || ( Electron_mvaNoIso_WP90 && Electron_pfRelIso03_all<0.3) )
                     && (HttCandidate.isLeg(Electron_idx, Leg::e)== false)""")
    df = df.Filter("Electron_idx[Electron_vetoSel].size() == 0", "No extra electrons")
    df = df.Define("Muon_vetoSel",
                   f"""v_ops::pt(Muon_p4) > 10 && abs(v_ops::eta(Muon_p4)) < 2.4 && abs(Muon_dz) < 0.2 && abs(Muon_dxy) < 0.045
                      && ( Muon_mediumId || Muon_tightId ) && Muon_pfRelIso04_all<0.3
                      && (HttCandidate.isLeg(Muon_idx, Leg::mu) == false)""")
    df = df.Filter("Muon_idx[Muon_vetoSel].size() == 0", "No extra muons")
    return df

def RecoJetSelection(df):
    df = df.Define("Jet_bIncl", f"v_ops::pt(Jet_p4)>20 && abs(v_ops::eta(Jet_p4)) < 2.5 && ( Jet_jetId & 2 ) && (Jet_puId>0 || v_ops::pt(Jet_p4)>50)")
    df = df.Define("FatJet_bbIncl", "FatJet_msoftdrop > 30 && abs(v_ops::eta(FatJet_p4)) < 2.5")
    df = df.Define("Jet_bCand", "RemoveOverlaps(Jet_p4, Jet_bIncl,{{HttCandidate.leg_p4[0], HttCandidate.leg_p4[1]},}, 2, 0.5)")
    df = df.Define("FatJet_bbCand", "RemoveOverlaps(FatJet_p4, FatJet_bbIncl, {{HttCandidate.leg_p4[0], HttCandidate.leg_p4[1]},}, 2, 0.5)")
    return df

def ExtraRecoJetSelection(df):
    df = df.Define("ExtraJet_B0", f"v_ops::pt(Jet_p4)>20 && abs(v_ops::eta(Jet_p4)) < 5 && ( Jet_jetId & 2 ) && (Jet_puId>0 || v_ops::pt(Jet_p4)>50)")
    df = df.Define(f"ObjectsToRemoveOverlap", "if(Hbb_isValid){return std::vector<RVecLV>({{HttCandidate.leg_p4[0], HttCandidate.leg_p4[1],HbbCandidate->leg_p4[0],HbbCandidate->leg_p4[1]}}); } return std::vector<RVecLV>({{HttCandidate.leg_p4[0], HttCandidate.leg_p4[1]}})")
    df = df.Define(f"ExtraJet_B1", """ RemoveOverlaps(Jet_p4, ExtraJet_B0,ObjectsToRemoveOverlap, 2, 0.5)""")
    return df


def ApplyJetSelection(df):
    return df.Filter("Jet_idx[Jet_bCand].size()>=2 || FatJet_idx[FatJet_bbCand].size()>=1", "Reco bjet candidates")

def GenRecoJetMatching(df):
    df = df.Define("Jet_genJetIdx_matched", "GenRecoJetMatching(event,Jet_idx, GenJet_idx, Jet_bCand, GenJet_B2, GenJet_p4, Jet_p4 , 0.3)")
    df = df.Define("Jet_genMatched", "Jet_genJetIdx_matched>=0")
    return df.Filter("Jet_genJetIdx_matched[Jet_genMatched].size()>=2", "Two different gen-reco jet matches at least")

def DefineHbbCand(df):
    df = df.Define("Jet_HHBtagScore", "GetHHBtagScore(Jet_bCand, Jet_idx, Jet_p4,Jet_btagDeepFlavB, MET_pt,  MET_phi, HttCandidate, period, event)")
    df = df.Define("HbbCandidate", "GetHbbCandidate(Jet_HHBtagScore, Jet_bCand, Jet_p4, Jet_idx)")
    return df

def RecottHttCandidateSelection_ttHH(df):
    df = df.Define("Electron_iso", "Electron_pfRelIso03_all") \
           .Define("Muon_iso", "Muon_pfRelIso04_all") \
           .Define("Tau_iso", "-Tau_rawDeepTau2018v2p5VSjet")

    df = df.Define("Electron_ttHH", f"""
        v_ops::pt(Electron_p4) > 20 && abs(v_ops::eta(Electron_p4)) < 2.5
        && abs(Electron_dz) < 0.2 && abs(Electron_dxy) < 0.045
        && Electron_mvaNoIso_WP90 && Electron_pfRelIso03_all < 0.5
    """)
    df = df.Define("Electron_ttHH_tight", "Electron_ttHH && Electron_pfRelIso03_all < 0.15")

    df = df.Define("Muon_ttHH", f"""
        v_ops::pt(Muon_p4) > 20 && abs(v_ops::eta(Muon_p4)) < 2.4
        && abs(Muon_dz) < 0.2 && abs(Muon_dxy) < 0.045
        && Muon_mediumId && Muon_pfRelIso04_all < 0.5
    """)
    df = df.Define("Muon_ttHH_tight", "Muon_ttHH && Muon_pfRelIso04_all < 0.15")

    df = df.Define("Tau_ttHH", f"""
        v_ops::pt(Tau_p4) > 20 && abs(v_ops::eta(Tau_p4)) < 2.5 && abs(Tau_dz) < 0.2
        && Tau_decayMode != 5 && Tau_decayMode != 6
        && Tau_idDeepTau2018v2p5VSe >= {WorkingPointsTauVSe.VVLoose.value}
        && Tau_idDeepTau2018v2p5VSmu >= {WorkingPointsTauVSmu.VLoose.value}
        && Tau_idDeepTau2018v2p5VSjet >= {WorkingPointsTauVSjet.VVVLoose.value}
    """)
    df = df.Define("Tau_ttHH_tight", f"Tau_ttHH && Tau_idDeepTau2018v2p5VSjet >= {WorkingPointsTauVSjet.VLoose.value}")

    cand_columns = []
    ttHH_exclueded_channels = [
        'tauTauTauTau', 'muTauTauTau', 'eTauTauTau', 'tauTauTau',
    ]
    ttHH_channels = [
        'muMuMuMu', 'eMuMuMu', 'muMuMuTau', 'eEMuMu', 'eMuMuTau', 'muMuTauTau',
        'eEEMu', 'eEMuTau', 'eMuTauTau', 'muTauTauTau',
        'eEEE', 'eEETau', 'eETauTau', 'eTauTauTau', 'tauTauTauTau',
        'muMuMu', 'eMuMu', 'muMuTau', 'eEMu', 'eMuTau', 'muTauTau', 'eEE', 'eETau', 'eTauTau', 'tauTauTau',
        'muMu', 'eMu', 'muTau', 'eE', 'eTau', 'tauTau'
    ] # in order of importance during the channel selection
    for ch in ttHH_channels:
        if ch in ttHH_exclueded_channels: continue
        legs = getChannelLegs(ch)
        cand_column = f"HttCandidates_{ch}"
        leg_inputs = []
        for leg_idx, leg in enumerate(legs):
            sel_suffix = '_tight' if leg_idx < len(legs) - 1 else ''
            leg_inputs.extend([
                f"{leg}_ttHH{sel_suffix}", f"{leg}_p4", f"{leg}_iso", f"{leg}_charge", f"{leg}_genMatchIdx"
            ])
        leg_inputs_str = ', '.join(leg_inputs)
        df = df.Define(cand_column, f"GetHTTCandidates<4>(Channel::{ch}, 0.5, {leg_inputs_str})")
        cand_columns.append(cand_column)
    cand_filters = [ f'{c}.size() > 0' for c in cand_columns ]
    df = df.Filter(" || ".join(cand_filters), "At lease one HTT candidate")
    cand_list_str = ', '.join([ '&' + c for c in cand_columns])
    return df.Define('HttCandidate', f'GetBestHTTCandidate<4>({{ {cand_list_str} }}, event)')

def RecoJetSelection_ttHH(df):
    df = df.Define("Jet_ttHH_sel", f"v_ops::pt(Jet_p4)>20 && abs(v_ops::eta(Jet_p4)) < 5 && ( Jet_jetId & 2 ) && (Jet_puId>0 || v_ops::pt(Jet_p4)>50)")
    df = df.Define("FatJet_ttHH_sel", "FatJet_msoftdrop > 30 && abs(v_ops::eta(FatJet_p4)) < 5")
    df = df.Define("Jet_ttHH", "RemoveOverlaps(Jet_p4, Jet_ttHH_sel, HttCandidate.getLegP4s(), 0.5)")
    df = df.Define("Jet_bCand", "Jet_ttHH && abs(v_ops::eta(Jet_p4)) < 2.5")
    df = df.Define("FatJet_ttHH", "RemoveOverlaps(FatJet_p4, FatJet_ttHH_sel, HttCandidate.getLegP4s(), 0.5)")
    return df.Filter("Jet_idx[Jet_ttHH].size() + FatJet_idx[FatJet_ttHH].size() * 2 >= 4", "Reco jet candidates")
