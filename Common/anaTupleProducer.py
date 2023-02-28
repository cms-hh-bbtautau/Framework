import copy
import os
import sys
import ROOT

if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])

import Common.BaselineSelection as Baseline
import Common.Utilities as Utilities
import Common.ReportTools as ReportTools
import Common.triggerSel as Triggers
import Corrections.Corrections as Corrections
from Corrections.lumi import LumiFilter


#ROOT.EnableImplicitMT(1)

deepTauScores= ["rawDeepTau2017v2p1VSe","rawDeepTau2017v2p1VSmu",
            "rawDeepTau2017v2p1VSjet", "rawDeepTau2018v2p5VSe", "rawDeepTau2018v2p5VSmu",
            "rawDeepTau2018v2p5VSjet",
            "idDeepTau2017v2p1VSe", "idDeepTau2017v2p1VSjet", "idDeepTau2017v2p1VSmu",
            "idDeepTau2018v2p5VSe","idDeepTau2018v2p5VSjet","idDeepTau2018v2p5VSmu",
            "decayMode"]
JetObservables = ["particleNetAK4_B", "particleNetAK4_CvsB",
                "particleNetAK4_CvsL","particleNetAK4_QvsG","particleNetAK4_puIdDisc",
                "btagDeepFlavB","btagDeepFlavCvB","btagDeepFlavCvL"]
JetObservablesMC = ["hadronFlavour","partonFlavour"]

defaultColToSave = ["event","luminosityBlock","run", "sample_type", "period", "X_mass", "isData",
                "MET_pt", "MET_phi","PuppiMET_pt", "PuppiMET_phi",
                "DeepMETResolutionTune_pt", "DeepMETResolutionTune_phi","DeepMETResponseTune_pt", "DeepMETResponseTune_phi",
                "MET_covXX", "MET_covXY", "MET_covYY", "PV_npvs"]

class DataFrameWrapper:
    def __init__(self, df):
        self.df = df
        self.colToSave = copy.deepcopy(defaultColToSave)

    def Define(self, varToDefine, varToCall):
        self.df = self.df.Define(f"{varToDefine}", f"{varToCall}")

    def DefineAndAppend(self, varToDefine, varToCall):
        self.Define(varToDefine, varToCall)
        self.colToSave.append(varToDefine)

    def Apply(self, func, *args, **kwargs):
        result = func(self.df, *args, **kwargs)
        if isinstance(result, tuple):
            self.df = result[0]
            if len(result) == 2:
                return result[1]
            return result[1:]
        else:
            self.df = result

def addAllVariables(dfw, syst_name, isData, trigger_class):
    dfw.Apply(Baseline.SelectRecoP4, syst_name)
    dfw.Apply(Baseline.RecoLeptonsSelection)
    dfw.Apply(Baseline.RecoJetAcceptance)
    dfw.Apply(Baseline.RecoHttCandidateSelection, config["GLOBAL"])
    dfw.Apply(Baseline.RecoJetSelection)
    dfw.Apply(Baseline.RequestOnlyResolvedRecoJets)
    dfw.Apply(Baseline.ThirdLeptonVeto)
    dfw.Apply(Baseline.DefineHbbCand)
    if trigger_class is not None:
        hltBranches = dfw.Apply(trigger_class.ApplyTriggers, isData)
        dfw.colToSave.extend(hltBranches)
    dfw.Define(f"Tau_recoJetMatchIdx", f"FindMatching(Tau_p4, Jet_p4, 0.5)")
    dfw.Define(f"Muon_recoJetMatchIdx", f"FindMatching(Muon_p4, Jet_p4, 0.5)")
    dfw.Define( f"Electron_recoJetMatchIdx", f"FindMatching(Electron_p4, Jet_p4, 0.5)")
    dfw.DefineAndAppend("channelId","static_cast<int>(httCand.channel())")
    jet_obs = JetObservables
    if not isData:
        jet_obs.extend(JetObservablesMC)
        dfw.colToSave.append("LHE_HT")
    for leg_idx in [0,1]:
        dfw.DefineAndAppend( f"tau{leg_idx+1}_pt", f"static_cast<float>(httCand.leg_p4[{leg_idx}].Pt())")
        dfw.DefineAndAppend( f"tau{leg_idx+1}_eta", f"static_cast<float>(httCand.leg_p4[{leg_idx}].Eta())")
        dfw.DefineAndAppend(f"tau{leg_idx+1}_phi", f"static_cast<float>(httCand.leg_p4[{leg_idx}].Phi())")
        dfw.DefineAndAppend(f"tau{leg_idx+1}_mass", f"static_cast<float>(httCand.leg_p4[{leg_idx}].M())")
        dfw.DefineAndAppend(f"tau{leg_idx+1}_charge", f"httCand.leg_charge[{leg_idx}]")
        dfw.Define(f"tau{leg_idx+1}_idx", f"httCand.leg_index[{leg_idx}]")
        dfw.Define(f"tau{leg_idx+1}_genMatchIdx", f"httCand.leg_genMatchIdx[{leg_idx}]")
        dfw.Define(f"tau{leg_idx+1}_recoJetMatchIdx", f"FindMatching(httCand.leg_p4[{leg_idx}], Jet_p4, 0.3)")
        dfw.DefineAndAppend( f"tau{leg_idx+1}_iso", f"httCand.leg_rawIso.at({leg_idx})")
        for deepTauScore in deepTauScores:
            dfw.DefineAndAppend( f"tau{leg_idx+1}_{deepTauScore}",
                                     f"httCand.leg_type[{leg_idx}] == Leg::tau ? Tau_{deepTauScore}.at(httCand.leg_index[{leg_idx}]) : -1;")

        if not isData:
            dfw.DefineAndAppend(f"tau{leg_idx+1}_gen_kind", f"""tau{leg_idx+1}_genMatchIdx>=0 ? static_cast<int>(genLeptons.at(tau{leg_idx+1}_genMatchIdx).kind()) :
                                              static_cast<int>(GenLeptonMatch::NoMatch);""")
            dfw.DefineAndAppend(f"tau{leg_idx+1}_gen_vis_pt", f"""tau{leg_idx+1}_genMatchIdx>=0? static_cast<float>(genLeptons.at(tau{leg_idx+1}_genMatchIdx).visibleP4().Pt()) :
                                                    -1.f;""")
            dfw.DefineAndAppend(f"tau{leg_idx+1}_gen_vis_eta", f"""tau{leg_idx+1}_genMatchIdx>=0? static_cast<float>(genLeptons.at(tau{leg_idx+1}_genMatchIdx).visibleP4().Eta()) :
                                                        -1.f;""")
            dfw.DefineAndAppend(f"tau{leg_idx+1}_gen_vis_phi", f"""tau{leg_idx+1}_genMatchIdx>=0? static_cast<float>(genLeptons.at(tau{leg_idx+1}_genMatchIdx).visibleP4().Phi()) :
                                                        -1.f;""")
            dfw.DefineAndAppend(f"tau{leg_idx+1}_gen_vis_mass", f"""tau{leg_idx+1}_genMatchIdx>=0? static_cast<float>(genLeptons.at(tau{leg_idx+1}_genMatchIdx).visibleP4().M()) :
                                                    -1.f;""")
            dfw.DefineAndAppend(f"tau{leg_idx+1}_gen_nChHad", f"""tau{leg_idx+1}_genMatchIdx>=0? genLeptons.at(tau{leg_idx+1}_genMatchIdx).nChargedHadrons() : 0;""")
            dfw.DefineAndAppend(f"tau{leg_idx+1}_gen_nNeutHad", f"""tau{leg_idx+1}_genMatchIdx>=0? genLeptons.at(tau{leg_idx+1}_genMatchIdx).nNeutralHadrons() : 0;""")
            dfw.DefineAndAppend(f"tau{leg_idx+1}_gen_charge", f"""tau{leg_idx+1}_genMatchIdx>=0? genLeptons.at(tau{leg_idx+1}_genMatchIdx).charge() : -10;""")
            dfw.DefineAndAppend( f"tau{leg_idx+1}_seedingJet_partonFlavour",
                                        f"tau{leg_idx+1}_recoJetMatchIdx>=0 ? Jet_partonFlavour.at(tau{leg_idx+1}_recoJetMatchIdx) : -1;")
            dfw.DefineAndAppend( f"tau{leg_idx+1}_seedingJet_hadronFlavour",
                                    f"tau{leg_idx+1}_recoJetMatchIdx>=0 ? Jet_hadronFlavour.at(tau{leg_idx+1}_recoJetMatchIdx) : -1;")

        dfw.DefineAndAppend( f"tau{leg_idx+1}_seedingJet_pt",
                                    f"tau{leg_idx+1}_recoJetMatchIdx>=0 ? static_cast<float>(Jet_p4.at(tau{leg_idx+1}_recoJetMatchIdx).Pt()) : -1.f;")
        dfw.DefineAndAppend( f"tau{leg_idx+1}_seedingJet_eta",
                                    f"tau{leg_idx+1}_recoJetMatchIdx>=0 ? static_cast<float>(Jet_p4.at(tau{leg_idx+1}_recoJetMatchIdx).Eta()) : -1.f;")
        dfw.DefineAndAppend( f"tau{leg_idx+1}_seedingJet_phi",
                                    f"tau{leg_idx+1}_recoJetMatchIdx>=0 ? static_cast<float>(Jet_p4.at(tau{leg_idx+1}_recoJetMatchIdx).Phi()) : -1.f;")
        dfw.DefineAndAppend( f"tau{leg_idx+1}_seedingJet_mass",
                                    f"tau{leg_idx+1}_recoJetMatchIdx>=0 ? static_cast<float>(Jet_p4.at(tau{leg_idx+1}_recoJetMatchIdx).M()) : -1.f;")


        dfw.DefineAndAppend(f"b{leg_idx+1}_pt", f"static_cast<float>(HbbCandidate.leg_p4[{leg_idx}].Pt())")
        dfw.DefineAndAppend(f"b{leg_idx+1}_eta", f"static_cast<float>(HbbCandidate.leg_p4[{leg_idx}].Eta())")
        dfw.DefineAndAppend(f"b{leg_idx+1}_phi", f"static_cast<float>(HbbCandidate.leg_p4[{leg_idx}].Phi())")
        dfw.DefineAndAppend(f"b{leg_idx+1}_mass", f"static_cast<float>(HbbCandidate.leg_p4[{leg_idx}].M())")

        for jetVar in jet_obs:
            if(jetVar not in dfw.df.GetColumnNames()): continue
            dfw.DefineAndAppend(f"b{leg_idx+1}_{jetVar}", f"Jet_{jetVar}.at(HbbCandidate.leg_index[{leg_idx}])")
        dfw.DefineAndAppend(f"b{leg_idx+1}_HHbtag", f"static_cast<float>(Jet_HHBtagScore.at(HbbCandidate.leg_index[{leg_idx}]))")


def createAnatuple(inFile, outFile, config, sample_name, snapshotOptions,range, evtIds,
                   store_noncentral, compute_unc_variations):

    period = config["GLOBAL"]["era"]
    mass = -1 if 'mass' not in config[sample_name] else config[sample_name]['mass']
    isHH = True if mass > 0 else False
    isData = True if config[sample_name]['sampleType'] == 'data' else False
    Baseline.Initialize(True, True)
    if not isData:
        Corrections.Initialize(config=config['GLOBAL'])
    triggerFile = config['GLOBAL']['triggerFile']
    trigger_class = Triggers.Triggers(triggerFile) if triggerFile is not None else None
    df = ROOT.RDataFrame("Events", inFile)
    if range is not None:
        df = df.Range(range)
    if len(evtIds) > 0:
        df = df.Filter(f"static const std::set<ULong64_t> evts = {{ {evtIds} }}; return evts.count(event) > 0;")

    if isData and 'lumiFile' in config['GLOBAL']:
        lumiFilter = LumiFilter(config['GLOBAL']['lumiFile'])
        df = lumiFilter.filter(df)

    df = Baseline.applyMETFlags(df, config["GLOBAL"]["MET_flags"])
    df = df.Define("sample_type", f"static_cast<int>(SampleType::{config[sample_name]['sampleType']})")
    df = df.Define("period", f"static_cast<int>(Period::{period})")
    df = df.Define("X_mass", f"static_cast<int>({mass})")
    is_data = 'true' if isData else 'false'
    df = df.Define("isData", is_data)

    df = Baseline.CreateRecoP4(df)
    df = Baseline.DefineGenObjects(df, isData=isData, isHH=isHH)
    if isData:
        syst_dict = { 'nano' : 'Central' }
    else:
        df, syst_dict = Corrections.applyScaleUncertainties(df)
    for syst_name, source_name in syst_dict.items():
        is_central = syst_name in [ 'Central', 'nano' ]
        if not is_central and not compute_unc_variations: continue
        suffix = '' if is_central else f'_{syst_name}'
        if len(suffix) and not store_noncentral: continue
        dfw = DataFrameWrapper(df)
        addAllVariables(dfw, syst_name, isData, trigger_class)
        if not isData:
            weight_branches = dfw.Apply(Corrections.getNormalisationCorrections, config, sample_name, is_central and compute_unc_variations)
            weight_branches.extend(dfw.Apply(Corrections.trg.getTrgSF, trigger_class.trigger_dict.keys(), is_central and compute_unc_variations))
            dfw.colToSave.extend(weight_branches)
        report = dfw.df.Report()
        histReport = ReportTools.SaveReport(report.GetValue(), reoprtName=f"Report{suffix}")
        varToSave = Utilities.ListToVector(dfw.colToSave)
        dfw.df.Snapshot(f"Events{suffix}", outFile, varToSave, snapshotOptions)
        outputRootFile= ROOT.TFile(outFile, "UPDATE")
        outputRootFile.WriteTObject(histReport, f"Report{suffix}", "Overwrite")
        outputRootFile.Close()

if __name__ == "__main__":
    import argparse
    import os
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('--configFile', type=str)
    parser.add_argument('--inFile', type=str)
    parser.add_argument('--outFile', type=str)
    parser.add_argument('--sample', type=str)
    parser.add_argument('--compressionLevel', type=int, default=9)
    parser.add_argument('--compressionAlgo', type=str, default="LZMA")
    parser.add_argument('--nEvents', type=int, default=None)
    parser.add_argument('--evtIds', type=str, default='')
    parser.add_argument('--store-noncentral', action="store_true", help="Store ES variations.")
    parser.add_argument('--compute_unc_variations', type=bool, default=True)

    args = parser.parse_args()

    ROOT.gROOT.ProcessLine(".include "+ os.environ['ANALYSIS_PATH'])
    ROOT.gROOT.ProcessLine('#include "Common/GenTools.h"')
    isHH=False
    isData = False
    with open(args.configFile, 'r') as f:
        config = yaml.safe_load(f)

    if os.path.exists(args.outFile):
        os.remove(args.outFile)

    snapshotOptions = ROOT.RDF.RSnapshotOptions()
    snapshotOptions.fOverwriteIfExists=True
    snapshotOptions.fMode="UPDATE"
    snapshotOptions.fCompressionAlgorithm = getattr(ROOT.ROOT, 'k' + args.compressionAlgo)
    snapshotOptions.fCompressionLevel = args.compressionLevel
    createAnatuple(args.inFile, args.outFile, config, args.sample, snapshotOptions, args.nEvents,
                   args.evtIds, args.store_noncentral, args.compute_unc_variations)
