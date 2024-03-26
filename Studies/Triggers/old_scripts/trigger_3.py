
import copy
import datetime
import os
import sys
import ROOT
import shutil
import zlib
import time
import matplotlib.pyplot as plt

if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])


#ROOT.EnableImplicitMT(1)
ROOT.EnableThreadSafety()

import Common.Utilities as Utilities

ggR_samples = [
               "GluGluToRadionToHHTo2B2Tau_M-250", "GluGluToRadionToHHTo2B2Tau_M-260", "GluGluToRadionToHHTo2B2Tau_M-270", "GluGluToRadionToHHTo2B2Tau_M-280", "GluGluToRadionToHHTo2B2Tau_M-300", "GluGluToRadionToHHTo2B2Tau_M-320", "GluGluToRadionToHHTo2B2Tau_M-350", "GluGluToRadionToHHTo2B2Tau_M-450", "GluGluToRadionToHHTo2B2Tau_M-500", "GluGluToRadionToHHTo2B2Tau_M-550", "GluGluToRadionToHHTo2B2Tau_M-600", "GluGluToRadionToHHTo2B2Tau_M-650", "GluGluToRadionToHHTo2B2Tau_M-700", "GluGluToRadionToHHTo2B2Tau_M-750", "GluGluToRadionToHHTo2B2Tau_M-800", "GluGluToRadionToHHTo2B2Tau_M-850", "GluGluToRadionToHHTo2B2Tau_M-900", "GluGluToRadionToHHTo2B2Tau_M-1000", "GluGluToRadionToHHTo2B2Tau_M-1250", "GluGluToRadionToHHTo2B2Tau_M-1500", "GluGluToRadionToHHTo2B2Tau_M-1750", "GluGluToRadionToHHTo2B2Tau_M-2000", "GluGluToRadionToHHTo2B2Tau_M-2500", "GluGluToRadionToHHTo2B2Tau_M-3000"]

ggBG_samples = [
    'GluGluToBulkGravitonToHHTo2B2Tau_M-1000', 'GluGluToBulkGravitonToHHTo2B2Tau_M-1250', 'GluGluToBulkGravitonToHHTo2B2Tau_M-1500', 'GluGluToBulkGravitonToHHTo2B2Tau_M-1750', 'GluGluToBulkGravitonToHHTo2B2Tau_M-2000', 'GluGluToBulkGravitonToHHTo2B2Tau_M-250', 'GluGluToBulkGravitonToHHTo2B2Tau_M-2500', 'GluGluToBulkGravitonToHHTo2B2Tau_M-260', 'GluGluToBulkGravitonToHHTo2B2Tau_M-270', 'GluGluToBulkGravitonToHHTo2B2Tau_M-280', 'GluGluToBulkGravitonToHHTo2B2Tau_M-300', 'GluGluToBulkGravitonToHHTo2B2Tau_M-3000', 'GluGluToBulkGravitonToHHTo2B2Tau_M-320', 'GluGluToBulkGravitonToHHTo2B2Tau_M-350', 'GluGluToBulkGravitonToHHTo2B2Tau_M-400', 'GluGluToBulkGravitonToHHTo2B2Tau_M-450', 'GluGluToBulkGravitonToHHTo2B2Tau_M-500', 'GluGluToBulkGravitonToHHTo2B2Tau_M-550', 'GluGluToBulkGravitonToHHTo2B2Tau_M-600', 'GluGluToBulkGravitonToHHTo2B2Tau_M-650', 'GluGluToBulkGravitonToHHTo2B2Tau_M-700', 'GluGluToBulkGravitonToHHTo2B2Tau_M-750', 'GluGluToBulkGravitonToHHTo2B2Tau_M-800', 'GluGluToBulkGravitonToHHTo2B2Tau_M-850', 'GluGluToBulkGravitonToHHTo2B2Tau_M-900' ]

def addHistToFile(df, var, hist_name, out_file):
    x_bins = [0, 20, 40, 60, 80, 100, 125, 150,175, 200, 250, 300]
    x_bins_vec = Utilities.ListToVector(x_bins, "double")
    model = ROOT.RDF.TH1DModel(hist_name, hist_name, x_bins_vec.size()-1, x_bins_vec.data())
    histogram = df.Histo1D(model, var)
    histogram.Write()

regions = {
    'singleTau_region':
    {
        'eTau':'tau1_pt>20',
        'muTau':'tau1_pt>20',
        'tauTau':'(tau1_pt > 190 && abs(tau1_eta) < 2.1) || (tau2_pt > 190 && abs(tau2_eta) < 2.1)'
    },
    'other_trg_region':
    {
        'eTau':'!singleTau_region_eTau',
        'muTau':'!singleTau_region_muTau',
        #'tauTau':'tau1_pt > 40 && abs(tau1_eta) < 2.1 && tau2_pt > 40 && abs(tau2_eta) < 2.1 && !singleTau_region_tauTau'
        'tauTau':'tau1_pt > 40 && abs(tau1_eta) < 2.1 && tau2_pt > 40 && abs(tau2_eta) < 2.1'
    },
    'MET_region':
    {
        'eTau':'!singleTau_region_eTau',
        'muTau':'!singleTau_region_muTau',
        'tauTau':'met_pt > 180 && !(other_trg_region_tauTau || singleTau_region_tauTau)'
    }
}



if __name__ == "__main__":
    import argparse
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=False, type=str, default='GluGluToBulkGraviton')
    parser.add_argument('--deepTauVersion', required=False, type=str, default='v2p1')
    parser.add_argument('--mass', required=False, type=str, default='750')
    args = parser.parse_args()

    startTime = time.time()
    headers_dir = os.path.dirname(os.path.abspath(__file__))
    ROOT.gROOT.ProcessLine(f".include {os.environ['ANALYSIS_PATH']}")
    ROOT.gInterpreter.Declare(f'#include "include/Utilities.h"')

    #### useful stuff ####
    deepTauYear = '2017' if args.deepTauVersion=='v2p1' else '2018'
    tau2_iso_var = f"tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet"
    tau1_iso_var = f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet"
    filters = {
        "channel": "channelId==33",
        "gen_channel": "genchannelId ==33",
        "gen_kind": "tau1_gen_kind==5 && tau2_gen_kind==5",
        "OS":"tau1_charge*tau2_charge < 0",
        "Iso_2_loose": f"{tau2_iso_var} >=  {Utilities.WorkingPointsTauVSjet.Loose.value}" ,
        "Iso": f"{tau2_iso_var} >=  {Utilities.WorkingPointsTauVSjet.Medium.value}" ,
        #"Iso_1_tight": f"{tau1_iso_var} >=  {Utilities.WorkingPointsTauVSjet.Tight.value}" ,
    }

    inDir = "/afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/output/Run2_2018/{}/anaTuples/"
    outFileDir = "/afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/output/Run2_2018/{}/histograms/"
    masses = []

    sample_list = ggR_samples if args.sample == 'GluGluToRadion' else ggBG_samples

    for sample_name in sample_list:
        mass_string = sample_name.split('-')[-1]
        if mass_string != args.mass:continue
        mass_int = int(mass_string)
        masses.append(mass_int)
        inFile = os.path.join(inDir.format(sample_name), "nanoHTT_0.root")
        if not os.path.exists(inFile) :
            print(f"{inFile} does not exist")
            continue

        print(inFile)
        df_initial = ROOT.RDataFrame('Events', inFile)
        if not os.path.exists(outFileDir.format(sample_name)):
            os.mkdir(outFileDir.format(sample_name))
        outFile = os.path.join(outFileDir.format(sample_name), 'pt_spectrum.root')
        output_file = ROOT.TFile(outFile, "RECREATE")

        df = ROOT.RDataFrame('Events', inFile)

        addHistToFile(df, 'tau1_pt', f'tau1_pt_initial', output_file )
        addHistToFile(df, 'tau2_pt', f'tau2_pt_initial', output_file )

        print(f"nEventds initial : {df.Count().GetValue()}" )

        df_initial = df
        for filter_str in filters.keys():
            df = df.Filter(f""" {filters[filter_str]} """)
            print(f""" after {filter_str} = {df.Count().GetValue()}""")
            addHistToFile(df, 'tau1_pt', f'tau1_pt_{filter_str}', output_file )
            addHistToFile(df, 'tau2_pt', f'tau2_pt_{filter_str}', output_file )

        for region,reg_dict in regions.items():
            df = df.Define(f"{region}_tauTau", f"""{reg_dict["tauTau"]}""")




        pass_met_tauTau = "(HLT_MET && MET_region_tauTau) "
        pass_singleTau_tauTau = "(HLT_singleTau && singleTau_region_tauTau) "
        pass_diTau_tauTau = "(HLT_ditau && other_trg_region_tauTau ) "

        print(f""" after VALIDITY tauTau region {df.Filter("other_trg_region_tauTau").Count().GetValue()}""")
        addHistToFile(df.Filter("other_trg_region_tauTau"), 'tau1_pt', f'tau1_pt_validity', output_file )
        addHistToFile(df.Filter("other_trg_region_tauTau"), 'tau2_pt', f'tau2_pt_validity', output_file )

        print(f""" after TRG tauTau region {df.Filter("other_trg_region_tauTau && !singleTau_region_tauTau").Count().GetValue()}""")
        addHistToFile(df.Filter("other_trg_region_tauTau && !singleTau_region_tauTau"), 'tau1_pt', f'tau1_pt_application', output_file )
        addHistToFile(df.Filter("other_trg_region_tauTau && !singleTau_region_tauTau"), 'tau2_pt', f'tau2_pt_application', output_file )

        print(f""" after only HLT_ditau {df.Filter("HLT_ditau").Count().GetValue()}""")
        addHistToFile(df.Filter("HLT_ditau"), 'tau1_pt', f'tau1_pt_HLT_diTau', output_file )
        addHistToFile(df.Filter("HLT_ditau"), 'tau2_pt', f'tau2_pt_HLT_diTau', output_file )

        print(f""" after tauTau region && HLT_ditau {df.Filter(pass_diTau_tauTau).Count().GetValue()}""")
        addHistToFile(df.Filter(pass_diTau_tauTau), 'tau1_pt', f'tau1_pt_pass_diTau', output_file )
        addHistToFile(df.Filter(pass_diTau_tauTau), 'tau2_pt', f'tau2_pt_pass_diTau', output_file )

        output_file.Close()

    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))


