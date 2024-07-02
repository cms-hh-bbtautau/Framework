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
from Common.Setup import Setup
from Corrections.Corrections import Corrections
from Corrections.lumi import LumiFilter


#ROOT.EnableImplicitMT(1)
ROOT.EnableThreadSafety()

def SelectBTagShapeSF(df,weight_name):
    df = df.Define("weight_bTagShapeSF", weight_name)
    return df

def createAnatuple(inFile, treeName, outDir, setup, sample_name, anaCache, snapshotOptions,range, evtIds,
                   store_noncentral, compute_unc_variations, uncertainties, anaTupleDef):
    start_time = datetime.datetime.now()
    compression_settings = snapshotOptions.fCompressionAlgorithm * 100 + snapshotOptions.fCompressionLevel
    period = setup.global_params["era"]
    sample_config = setup.samples[sample_name]
    mass = -1 if 'mass' not in sample_config else sample_config['mass']
    spin = -100 if 'spin' not in sample_config else sample_config['spin']
    isHH = True if mass > 0 else False
    isData = True if sample_config['sampleType'] == 'data' else False
    loadTF = anaTupleDef.loadTF
    loadHHBtag = anaTupleDef.loadHHBtag
    lepton_legs = anaTupleDef.lepton_legs
    Baseline.Initialize(loadTF, loadHHBtag)
    Corrections.initializeGlobal(setup.global_params, isData=isData, load_corr_lib=True)
    corrections = Corrections.getGlobal()
    triggerFile = setup.global_params.get('triggerFile')
    if triggerFile is not None:
        triggerFile = os.path.join(os.environ['ANALYSIS_PATH'], triggerFile)
        trigger_class = Triggers.Triggers(triggerFile)
    else:
        trigger_class = None
    df = ROOT.RDataFrame(treeName, inFile)
    ROOT.RDF.Experimental.AddProgressBar(df)
    if range is not None:
        df = df.Range(range)
    if len(evtIds) > 0:
        df = df.Filter(f"static const std::set<ULong64_t> evts = {{ {evtIds} }}; return evts.count(event) > 0;")

    if isData and 'lumiFile' in setup.global_params:
        lumiFilter = LumiFilter(setup.global_params['lumiFile'])
        df = lumiFilter.filter(df)

    if "MET_flags" in setup.global_params:
        df = Baseline.applyMETFlags(df, setup.global_params["MET_flags"])
    df = df.Define("sample_type", f"static_cast<int>(SampleType::{sample_config['sampleType']})")
    df = df.Define("sample_name", f"static_cast<int>({zlib.crc32(sample_name.encode())})")
    isSignal = sample_config['sampleType'] in setup.global_params['signal_types']
    print("isSignal? ", isSignal)
    df = df.Define("period", f"static_cast<int>(Period::{period})")
    df = df.Define("X_mass", f"static_cast<int>({mass})")
    df = df.Define("X_spin", f"static_cast<int>({spin})")
    df = df.Define("entryIndex", "static_cast<int>(rdfentry_)")
    is_data = 'true' if isData else 'false'
    df = df.Define("isData", is_data)

    df = Baseline.CreateRecoP4(df, nano_version=setup.global_params['nano_version'])
    df = Baseline.DefineGenObjects(df, isData=isData, isHH=isHH)

    if isData:
        syst_dict = { 'nano' : 'Central' }
    else:
        ana_reco_objects = Baseline.ana_reco_object_collections[setup.global_params['nano_version']]
        df, syst_dict = corrections.applyScaleUncertainties(df, ana_reco_objects)
    df_empty = df
    snaps = []
    reports = []
    outfilesNames = []
    k=0
    for syst_name, source_name in syst_dict.items():
        if source_name not in uncertainties and "all" not in uncertainties: continue
        is_central = syst_name in [ 'Central', 'nano' ]
        if not is_central and not compute_unc_variations: continue
        suffix = '' if is_central else f'_{syst_name}'
        if len(suffix) and not store_noncentral: continue
        dfw = Utilities.DataFrameWrapper(df_empty, anaTupleDef.getDefaultColumnsToSave(isData))

        anaTupleDef.addAllVariables(dfw, syst_name, isData, trigger_class, lepton_legs, isSignal, setup.global_params)

        if setup.global_params['nano_version'] == 'v12':
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
            weight_branches = dfw.Apply(corrections.getNormalisationCorrections, setup.global_params,
                                        setup.samples, sample_name, lepton_legs,
                                        return_variations=is_central and compute_unc_variations, isCentral=is_central,
                                        ana_cache=anaCache)
            if 'trg' in corrections.to_apply:
                weight_branches.extend(dfw.Apply(corrections.trg.getTrgSF, trigger_class.trigger_dict.keys(), lepton_legs,
                                             is_central and compute_unc_variations, is_central))
            if 'btag' in corrections.to_apply:
                SF_branches_core,SF_weight_jes=dfw.Apply(corrections.jet.getBtagShapeSFs, syst_name, is_central,
                                                     compute_unc_variations)
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
                    if setup.global_params["storeExtraJets"]:
                        dfw.colToSave.append(f"""ExtraJet_{new_branch_name}""")
                    for bjet_idx in [1,2]:
                        dfw.DefineAndAppend(f"{new_branch_name}_b{bjet_idx}", f"Hbb_isValid ? {puIDbranch}[b{bjet_idx}_idx] : -100.f")
                if puIDbranch in weight_branches: weight_branches.remove(puIDbranch)
            dfw.colToSave.extend(weight_branches)
        reports.append(dfw.df.Report())
        varToSave = Utilities.ListToVector(dfw.colToSave)
        outfile_prefix = inFile.split('/')[-1]
        outfile_prefix = outfile_prefix.split('.')[0]
        outFileName = os.path.join(outDir, f"{outfile_prefix}{suffix}.root")
        outfilesNames.append(outFileName)
        snaps.append(dfw.df.Snapshot(f"Events", outFileName, varToSave, snapshotOptions))
    if snapshotOptions.fLazy == True:
        ROOT.RDF.RunGraphs(snaps)
    hist_time = ROOT.TH1D(f"time", f"time", 1, 0, 1)
    end_time = datetime.datetime.now()
    hist_time.SetBinContent(1, (end_time - start_time).total_seconds())
    # for index,fileName in enumerate(outfilesNames):
    #     outputRootFile= ROOT.TFile(fileName, "UPDATE", "", compression_settings)
    #     rep = ReportTools.SaveReport(reports[index].GetValue(), reoprtName=f"Report")
    #     outputRootFile.WriteTObject(rep, f"Report", "Overwrite")
    #     if index==0:
    #         outputRootFile.WriteTObject(hist_time, f"runtime", "Overwrite")
    #     outputRootFile.Close()
    #     # if print_cutflow:
    #     #     report.Print()

if __name__ == "__main__":
    import argparse
    import os
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('--period', required=True, type=str)
    parser.add_argument('--inFile', required=True, type=str)
    parser.add_argument('--outDir', required=True, type=str)
    parser.add_argument('--sample', required=True, type=str)
    parser.add_argument('--anaCache', required=True, type=str)
    parser.add_argument('--anaTupleDef', required=True, type=str)
    parser.add_argument('--store-noncentral', action="store_true", help="Store ES variations.")
    parser.add_argument('--compute-unc-variations', action="store_true")
    parser.add_argument('--uncertainties', type=str, default="all")
    parser.add_argument('--customisations', type=str, default=None)
    parser.add_argument('--treeName', required=False, type=str, default="Events")
    parser.add_argument('--particleFile', type=str,
                        default=f"{os.environ['ANALYSIS_PATH']}/config/pdg_name_type_charge.txt")
    parser.add_argument('--compressionLevel', type=int, default=4)
    parser.add_argument('--compressionAlgo', type=str, default="ZLIB")
    parser.add_argument('--nEvents', type=int, default=None)
    parser.add_argument('--evtIds', type=str, default='')

    args = parser.parse_args()

    ROOT.gROOT.ProcessLine(".include "+ os.environ['ANALYSIS_PATH'])
    ROOT.gROOT.ProcessLine('#include "include/GenTools.h"')
    ROOT.gInterpreter.ProcessLine(f"ParticleDB::Initialize(\"{args.particleFile}\");")
    setup = Setup.getGlobal(os.environ['ANALYSIS_PATH'], args.period, args.customisations)
    with open(args.anaCache, 'r') as f:
        anaCache = yaml.safe_load(f)

    anaTupleDef = Utilities.load_module(args.anaTupleDef)
    if os.path.isdir(args.outDir):
        shutil.rmtree(args.outDir)
    os.makedirs(args.outDir, exist_ok=True)
    snapshotOptions = ROOT.RDF.RSnapshotOptions()
    snapshotOptions.fOverwriteIfExists=False
    snapshotOptions.fLazy = True
    snapshotOptions.fMode="RECREATE"
    snapshotOptions.fCompressionAlgorithm = getattr(ROOT.ROOT, 'k' + args.compressionAlgo)
    snapshotOptions.fCompressionLevel = args.compressionLevel
    createAnatuple(args.inFile, args.treeName, args.outDir, setup, args.sample, anaCache, snapshotOptions,
                   args.nEvents, args.evtIds, args.store_noncentral, args.compute_unc_variations,
                   args.uncertainties.split(","), anaTupleDef)
