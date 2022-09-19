import ROOT
import os
from scipy import stats
import numpy as np
import enum

initialized = False

def Initialize():
    global initialized
    if not initialized:
        import os
        headers_dir = os.path.dirname(os.path.abspath(__file__))
        header_path_Gen = os.path.join(headers_dir, "BaselineGenSelection.h")
        header_path_Reco = os.path.join(headers_dir, "BaselineRecoSelection.h")
        ROOT.gInterpreter.Declare(f'#include "{header_path_Gen}"')
        ROOT.gInterpreter.Declare(f'#include "{header_path_Reco}"')
        initialized = True

leg_names = [ "Electron", "Muon", "Tau", "boostedTau" ]

channels = [ 'muMu', 'eMu', 'eE', 'muTau', 'eTau', 'tauTau' ] # in order of importance during the channel selection

channelLegs = {
    "eTau": [ "Electron", "Tau" ],
    "muTau": [ "Muon", "Tau" ],
    "tauTau": [ "Tau", "Tau" ],
    "muMu": [ "Muon", "Muon" ],
    "eMu": [ "Electron", "Muon" ],
    "eE": [ "Electron", "Electron" ],
}

class WorkingPointsTauVSmu:
    VLoose = 1
    Loose = 2
    Medium = 3
    Tight = 4

class WorkingPointsTauVSjet:
   VVVLoose = 1
   VVLoose = 2
   VLoose = 3
   Loose = 4
   Medium = 5
   Tight = 6
   VTight = 7
   VVTight = 8

class WorkingPointsTauVSe:
    VVVLoose = 1
    VVLoose = 2
    VLoose = 3
    Loose = 4
    Medium = 5
    Tight = 6
    VTight = 7
    VVTight = 8

class WorkingPointsBoostedTauVSjet:
   VVLoose = 1
   VLoose = 2
   Loose = 3
   Medium = 4
   Tight = 5
   VTight = 6
   VVTight = 7


def DefineGenObjects(df, Hbb_AK4mass_mpv):
    df = df.Define("GenPart_daughters", "GetDaughters(GenPart_genPartIdxMother)")
    df = df.Define("genHttCand", """GetGenHTTCandidate(event, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags,
                                                       GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)""")
    df = df.Define("genHbbIdx", """GetGenHBBIndex(event, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags)""")
    for var in ["GenJet", "GenJetAK8"]:
        df = df.Define(f"{var}_idx", f"CreateIndexes({var}_pt.size())")
        df = df.Define(f"{var}_p4", f"GetP4({var}_pt,{var}_eta,{var}_phi,{var}_mass, {var}_idx)")

    df = df.Define("GenJet_b_PF", "abs(GenJet_partonFlavour)==5")
    df = df.Define("GenJetAK8_b_PF", "abs(GenJetAK8_partonFlavour)==5")
    df = df.Define("GenJet_Hbb",f"FindTwoJetsClosestToMPV({Hbb_AK4mass_mpv}, GenJet_p4, GenJet_b_PF)")
    df = df.Define("GenJetAK8_Hbb", "FindGenJetAK8(GenJetAK8_mass, GenJetAK8_b_PF)")
    df = df.Define("genHbb_isBoosted", "GenPart_pt[genHbbIdx]>550")
    return df

def ApplyGenBaseline0(df):
    return df.Filter("PassGenAcceptance(genHttCand)", "GenHttCand Acceptance")

def ApplyGenBaseline1(df):
    df = df.Define("GenJet_B1","GenJet_pt > 20 && abs(GenJet_eta) < 2.5 && GenJet_Hbb")
    df = df.Define("GenJetAK8_B1","GenJetAK8_pt > 170 && abs(GenJetAK8_eta) < 2.5 && GenJetAK8_Hbb")
    return df.Filter("GenJet_idx[GenJet_B1].size()==2 || (GenJetAK8_idx[GenJetAK8_B1].size()==1 && genHbb_isBoosted)", "(One)Two b-parton (Fat)jets at least")

def ApplyGenBaseline2(df):
    for var in ["GenJet", "GenJetAK8"]:
        df = df.Define(f"{var}_B2", f"RemoveOverlaps({var}_p4, {var}_B1,{{{{genHttCand.leg_p4[0], genHttCand.leg_p4[1]}},}}, 2, 0.5)" )
    return df.Filter("GenJet_idx[GenJet_B2].size()==2 || (GenJetAK8_idx[GenJetAK8_B2].size()==1 && genHbb_isBoosted)", "No overlap between genJets and genHttCands")

def ApplyGenBaseline3(df):
    return df.Filter("GenJet_idx[GenJet_B2].size()==2", "Resolved topology")


def ApplyRecoBaseline0(df, apply_filter=True):
    for obj in [ "Electron", "Muon", "Tau", "Jet", "FatJet", "boostedTau" ]:
        df = df.Define(f"{obj}_idx", f"CreateIndexes({obj}_pt.size())") \
               .Define(f"{obj}_p4", f"GetP4({obj}_pt, {obj}_eta, {obj}_phi, {obj}_mass, {obj}_idx)")

    df = df.Define("Electron_B0", """
        Electron_pt > 18 && abs(Electron_eta) < 2.3 && abs(Electron_dz) < 0.2 && abs(Electron_dxy) < 0.045
        && (Electron_mvaIso_WP90 || (Electron_mvaNoIso_WP90 && Electron_pfRelIso03_all < 0.5))
    """)

    df = df.Define("Muon_B0", """
        Muon_pt > 18 && abs(Muon_eta) < 2.3 && abs(Muon_dz) < 0.2 && abs(Muon_dxy) < 0.045
        && ( ((Muon_tightId || Muon_mediumId) && Muon_pfRelIso04_all < 0.5) || (Muon_highPtId && Muon_tkRelIso < 0.5) )
    """)

    df = df.Define("Tau_B0", f"""
        Tau_pt > 15 && abs(Tau_eta) < 2.5 && abs(Tau_dz) < 0.2 && Tau_decayMode != 5 && Tau_decayMode != 6
        && (    (    Tau_idDeepTau2017v2p1VSe >= {WorkingPointsTauVSe.VVLoose}
                  && Tau_idDeepTau2017v2p1VSmu >= {WorkingPointsTauVSmu.VLoose}
                  && Tau_idDeepTau2017v2p1VSjet >= {WorkingPointsTauVSjet.VVVLoose} )
             || (    Tau_idDeepTau2018v2p5VSe >= {WorkingPointsTauVSe.VVLoose}
                  && Tau_idDeepTau2018v2p5VSmu >= {WorkingPointsTauVSmu.VLoose}
                  && Tau_idDeepTau2018v2p5VSjet >= {WorkingPointsTauVSjet.VVVLoose} )
           )
    """)

    df = df.Define("boostedTau_B0", f"""
        boostedTau_pt > 40 && abs(boostedTau_eta) < 2.3 && abs(boostedTau_dz) < 0.2 && boostedTau_decayMode != 5
        && boostedTau_decayMode != 6 && boostedTau_idMVAnewDM2017v2 >= {WorkingPointsBoostedTauVSjet.VVLoose}
    """)

    df = df.Define("Electron_B0T", """
        Electron_B0 && (Electron_mvaIso_WP80
                        || (Electron_mvaNoIso_WP80 && Electron_pfRelIso03_all < 0.15))
    """)

    df = df.Define("Muon_B0T", """
        Muon_B0 && ( ((Muon_tightId || Muon_mediumId) && Muon_pfRelIso04_all < 0.15)
                    || (Muon_highPtId && Muon_tkRelIso < 0.15) )
    """)

    df = df.Define("Tau_B0T", f"""
        Tau_B0 && (
                      Tau_idDeepTau2017v2p1VSjet >= {WorkingPointsTauVSjet.Medium}
                   || Tau_idDeepTau2018v2p5VSjet >= {WorkingPointsTauVSjet.Medium} )
    """)

    df = df.Define("boostedTau_B0T", f"""
        boostedTau_B0 && boostedTau_idMVAnewDM2017v2 >= {WorkingPointsBoostedTauVSjet.Medium}
    """)

    ch_filters = []
    for leg1_idx in range(len(leg_names)):
        for leg2_idx in range(max(1, leg1_idx), len(leg_names)):
            leg1, leg2 = leg_names[leg1_idx], leg_names[leg2_idx]
            ch_filter = f"{leg1}{leg2}_B0"
            ch_filters.append(ch_filter)
            if leg1 == leg2:
                ch_filter_def = f"{leg1}_idx[{leg1}_B0].size() > 1 && {leg1}_idx[{leg1}_B0T].size() > 0"
            else:
                ch_filter_def = f"""
                    ({leg1}_idx[{leg1}_B0].size() > 0 && {leg2}_idx[{leg2}_B0T].size() > 0)
                    || ({leg1}_idx[{leg1}_B0T].size() > 0 && {leg2}_idx[{leg2}_B0].size() > 0)
                """
            df = df.Define(ch_filter, ch_filter_def)
    filter_expr = " || ".join(ch_filters)
    if apply_filter:
        return df.Filter(filter_expr, "Reco leptons requirements")
    else:
        return df, filter_expr


def ApplyRecoBaseline1(df, apply_filter=True):
    df = df.Define("Jet_B1", f"Jet_pt>15 && abs(Jet_eta) < 2.5 && ( Jet_jetId & 2 )")
    df = df.Define("FatJet_B1", "FatJet_msoftdrop > 30 && abs(FatJet_eta) < 2.5")

    df = df.Define("Lepton_p4_B0", "std::vector<RVecLV>{Electron_p4[Electron_B0], Muon_p4[Muon_B0], Tau_p4[Tau_B0]}")
    df = df.Define("Jet_B1T", "RemoveOverlaps(Jet_p4, Jet_B1, Lepton_p4_B0, 2, 0.5)")
    df = df.Define("FatJet_B1T", "RemoveOverlaps(FatJet_p4, FatJet_B1, Lepton_p4_B0, 2, 0.5)")

    filter_expr = "Jet_idx[Jet_B1T].size() >= 1 || FatJet_idx[FatJet_B1T].size() >= 1"
    if apply_filter:
        return df.Filter(filter_expr, "Reco Jet Acceptance")
    else:
        return df, filter_expr


def ApplyRecoBaseline2(df):
    df = df.Define("Electron_iso", "Electron_pfRelIso03_all") \
           .Define("Muon_iso", "Muon_pfRelIso04_all") \
           .Define("Tau_iso", "-Tau_rawDeepTau2017v2p1VSjet")

    df = df.Define("Electron_B2_eTau_1", "Electron_B0 && Electron_pt > 20 && Electron_mvaIso_WP80")
    df = df.Define("Tau_B2_eTau_2", f"""
        Tau_B0 && Tau_pt > 20
        && (Tau_idDeepTau2017v2p1VSe & {WorkingPointsTauVSe.VLoose})
        && (Tau_idDeepTau2017v2p1VSmu & {WorkingPointsTauVSmu.Tight})
    """)

    df = df.Define("Muon_B2_muTau_1", """
        Muon_B0 && Muon_pt > 20 && (   (Muon_tightId && Muon_pfRelIso04_all < 0.15)
                                    || (Muon_highPtId && Muon_tkRelIso < 0.15) )
    """)
    df = df.Define("Tau_B2_muTau_2", f"""
        Tau_B0 && Tau_pt > 20
        && (Tau_idDeepTau2017v2p1VSe & {WorkingPointsTauVSe.VLoose})
        && (Tau_idDeepTau2017v2p1VSmu & {WorkingPointsTauVSmu.Tight})
    """)

    df = df.Define("Tau_B2_tauTau_1", f"""
        Tau_B0 && Tau_pt > 20
        && (Tau_idDeepTau2017v2p1VSe & {WorkingPointsTauVSe.VVLoose})
        && (Tau_idDeepTau2017v2p1VSmu & {WorkingPointsTauVSmu.VLoose})
        && (Tau_idDeepTau2017v2p1VSjet & {WorkingPointsTauVSjet.Medium})
    """)

    df = df.Define("Tau_B2_tauTau_2", f"""
        Tau_B0 && Tau_pt > 20
        && (Tau_idDeepTau2017v2p1VSe & {WorkingPointsTauVSe.VVLoose})
        && (Tau_idDeepTau2017v2p1VSmu & {WorkingPointsTauVSmu.VLoose})
    """)

    df = df.Define("Muon_B2_muMu_1", """
        Muon_B0 && Muon_pt > 20 && (   (Muon_tightId && Muon_pfRelIso04_all < 0.15)
                                    || (Muon_highPtId && Muon_tkRelIso < 0.15) )
    """)
    df = df.Define("Muon_B2_muMu_2", """
        Muon_B0 && Muon_pt > 20 && (   (Muon_tightId && Muon_pfRelIso04_all < 0.3)
                                    || (Muon_highPtId && Muon_tkRelIso < 0.3) )
    """)

    df = df.Define("Electron_B2_eMu_1", """
        Electron_B0 && Electron_pt > 20 && Electron_mvaNoIso_WP80 && Electron_pfRelIso03_all < 0.3
    """)
    df = df.Define("Muon_B2_eMu_2", """
        Muon_B0 && Muon_pt > 20 && (   (Muon_tightId && Muon_pfRelIso04_all < 0.15)
                                    || (Muon_highPtId && Muon_tkRelIso < 0.15) )
    """)

    df = df.Define("Electron_B2_eE_1", """
        Electron_B0 && Electron_pt > 20
        && (Electron_mvaIso_WP80 || Electron_mvaNoIso_WP80 && Electron_pfRelIso03_all < 0.15)
    """)
    df = df.Define("Electron_B2_eE_2", """
        Electron_B0 && Electron_pt > 20 && Electron_mvaNoIso_WP80 && Electron_pfRelIso03_all < 0.3
    """)

    cand_columns = []
    for ch in channels:
        leg1, leg2 = channelLegs[ch]
        cand_column = f"httCands_{ch}"
        df = df.Define(cand_column, f"""
            GetHTTCandidates(Channel::{ch}, 0.4, {leg1}_B2_{ch}_1, {leg1}_p4, {leg1}_iso, {leg1}_charge,
                                                 {leg2}_B2_{ch}_2, {leg2}_p4, {leg2}_iso, {leg2}_charge)
        """)
        cand_columns.append(cand_column)
    cand_filters = [ f'{c}.size() > 0' for c in cand_columns ]
    df = df.Filter(" || ".join(cand_filters), "Reco Baseline 2")
    cand_list_str = ', '.join([ '&' + c for c in cand_columns])
    return df.Define('httCand', f'GetBestHTTCandidate({{ {cand_list_str} }})')



def ApplyRecoBaseline3(df):
    df = df.Define("Jet_B3T", "RemoveOverlaps(Jet_p4, Jet_B1T,{{httCand.leg_p4[0], httCand.leg_p4[1]},}, 2, 0.5)")
    df = df.Define("FatJet_B3T", "RemoveOverlaps(FatJet_p4, FatJet_B1T,{{httCand.leg_p4[0], httCand.leg_p4[1]},}, 2, 0.5)")
    return df.Filter("Jet_idx[Jet_B3T].size()>=2 || FatJet_idx[FatJet_B3T].size()>=1", "Reco Baseline 3")


def ApplyRecoBaseline4(df):
    return df.Filter("Jet_idx[Jet_B3T].size()>=2", "Reco Baseline 4")


def ApplyGenRecoJetMatching(df):
    df = df.Define("Jet_genJetIdx_matched", "GenRecoJetMatching(event,Jet_idx, GenJet_idx, Jet_B3T, GenJet_B2, GenJet_p4, Jet_p4 , 0.3)")
    df = df.Define("Jet_genMatched", "Jet_genJetIdx_matched>=0")
    return df.Filter("Jet_genJetIdx_matched[Jet_genMatched].size()>=2", "Two different gen-reco jet matches at least")

