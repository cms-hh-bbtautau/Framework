
import copy
import datetime
import os
import sys
import ROOT
import shutil
import zlib
import time
import numpy as np

from array import array
import matplotlib.pyplot as plt

if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])


#ROOT.EnableImplicitMT(1)
ROOT.EnableThreadSafety()

import Common.Utilities as Utilities
from Analysis.HistHelper import *
from Analysis.hh_bbtautau import *
from Studies.Triggers.Trig_utilities import *
ROOT.gStyle.SetOptStat(0)

if __name__ == "__main__":
    import argparse
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=False, type=str, default='GluGluToBulkGraviton')
    parser.add_argument('--deepTauVersion', required=False, type=str, default='v2p1')
    parser.add_argument('--wantBigTau', required=False, type=bool, default=False)
    parser.add_argument('--wantLegend', required=False, type=bool, default=False)
    parser.add_argument('--mass', required=False, type=str, default="700")
    parser.add_argument('--deepTauWP', required=False, type=str, default='medium')


    args = parser.parse_args()

    startTime = time.time()
    headers_dir = os.path.dirname(os.path.abspath(__file__))
    ROOT.gROOT.ProcessLine(f".include {os.environ['ANALYSIS_PATH']}")
    ROOT.gInterpreter.Declare(f'#include "include/KinFitNamespace.h"')
    ROOT.gInterpreter.Declare(f'#include "include/HistHelper.h"')
    ROOT.gInterpreter.Declare(f'#include "include/Utilities.h"')
    ROOT.gROOT.ProcessLine('#include "include/AnalysisTools.h"')

    #### useful stuff ####

    #inDir = "/eos/home-k/kandroso/cms-hh-bbtautau/anaTuples/Run2_2018/v1_deepTau2p1_HTT/"
    inDir = "/eos/home-k/kandroso/cms-hh-bbtautau/anaTuples/Run2_2018/v2_deepTau2p1_HTT/"
    #inDir = "/afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/output/Run2_2018/{}/anaTuples/"
    masses = []
    x_bins = [40,50, 60,70, 80,90, 100,120, 140,160,180,200, 250, 300,400,800,1000]
    x_bins_vec = Utilities.ListToVector(x_bins, "double")

    #### OS_ISO requirements ####

    deepTauYear = '2017' if args.deepTauVersion=='v2p1' else '2018'

    os_iso_filtering = {
        'VLoose': f"""(tau1_charge * tau2_charge < 0) && tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VLoose.value}""",
        'Loose': f"""(tau1_charge * tau2_charge < 0) && tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Loose.value}""",
        'Medium':f"""(tau1_charge * tau2_charge < 0) && tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Medium.value}"""
        }


    ########## dataFrame creation ##########
    sample_list = ggR_samples if args.sample == 'GluGluToRadion' else ggBG_samples
    colors = ['red', 'blue']
    labels = ['HLT_IsoMu24','HLT_Mu50']
    inFiles = []
    for sample_name in sample_list:
        mass_string = sample_name.split('-')[-1]
        mass_int = int(mass_string)
        masses.append(mass_int)
        inFile = os.path.join(inDir,f"{sample_name}", "nanoHTT_0.root")
        if not os.path.exists(inFile) :
            print(f"{inFile} does not exist")
            continue
        inFiles.append(inFile)
    inFiles_vec= Utilities.ListToVector(inFiles)
    dfWrapped_central  = DataFrameBuilder(ROOT.RDataFrame('Events', inFiles_vec), args.deepTauVersion)
    PrepareDfWrapped(dfWrapped_central)
    df_muTau = dfWrapped_central.df.Filter('muTau').Filter("tau1_gen_kind==4 && tau2_gen_kind==5")
    nInitial_muTau = df_muTau.Count().GetValue()
    print(nInitial_muTau)

    df_muTau = df_muTau.Filter(os_iso_filtering['Medium'])

    regions_expr= get_regions_dict('muTau', [190,190], 180, False)
    for reg_key,reg_exp in regions_expr.items():
        df_muTau = df_muTau.Define(reg_key, reg_exp)

    pass_MET_application = "(HLT_MET && MET_region && !(singleTau_region) && ! (other_trg_region))"
    pass_singleTau_application = "(HLT_singleTau && singleTau_region && !(MET_region) && ! (other_trg_region))"
    pass_other_trg_application= "({} && other_trg_region  && !(MET_region) && ! (singleTau_region) )  "

    pass_singleTau_validity = f"""HLT_singleTau && ({GetTrgValidityRegion("HLT_singleTau")})"""
    pass_MET_validity = f"""HLT_MET && ({GetTrgValidityRegion("HLT_MET")})"""


    #### muTau efficiencies ####
    pass_mutau_application = pass_other_trg_application.format("HLT_mutau")
    pass_singleMu_application = pass_other_trg_application.format("HLT_singleMu")
    pass_singleMu50_application = pass_other_trg_application.format("HLT_singleMu50")
    pass_mutau_validity = f"""HLT_mutau && ({GetTrgValidityRegion("HLT_mutau")})"""
    pass_singleMu_validity = f"""HLT_singleMu && ({GetTrgValidityRegion("HLT_singleMu")})"""
    pass_singleMu50_validity = f"""HLT_singleMu50 && ({GetTrgValidityRegion("HLT_singleMu50")})"""
    k=0
    for trig in [pass_singleMu_validity, pass_singleMu50_validity]:
        color_plt = colors[k]
        label_plt = labels[k]
        model_hist_denumerator = ROOT.RDF.TH1DModel("tau1_pt_denum", "tau1_pt_denum", x_bins_vec.size()-1, x_bins_vec.data())
        histo_denumerator = df_muTau.Histo1D(model_hist_denumerator,"tau1_pt")
        model_hist_numerator = ROOT.RDF.TH1DModel("tau1_pt_num", "tau1_pt_num", x_bins_vec.size()-1, x_bins_vec.data())
        histo_numerator = df_muTau.Filter(trig).Histo1D(model_hist_numerator,"tau1_pt")
        #histo_efficiency = ROOT.TH1D("efficiency", "efficiency",x_bins_vec.size()-1, x_bins_vec.data())
        n = len(x_bins)
        outFile = ROOT.TFile(f"Studies/Triggers/Mu50_efficiencies/muPt50_{mass_string}.root", "RECREATE")
        histo_numerator.Write()
        histo_denumerator.Write()
        outFile.Close()
        efficiencies = []
        x_points = []
        errors_x = []
        errors_y_u = []
        errors_y_l = []

        for i in range(0, histo_numerator.GetNbinsX()+1):
            numerator = histo_numerator.GetBinContent(i)
            denominator = histo_denumerator.GetBinContent(i)
            x_pt = histo_denumerator.GetBinCenter(i)
            efficiency = numerator / denominator if denominator != 0 else 0
            if efficiency == 0 : continue
            conf_interval = proportion_confint(numerator, denominator, alpha=0.05, method='normal')
            error_low = efficiency - conf_interval[0]
            error_up = conf_interval[1] - efficiency
            error_low = np.minimum(error_low, 1)
            error_up = np.minimum(error_up, 1)
            efficiencies.append(efficiency)
            #errors_x.append(((x_bins[i] - x_bins[i-1]) / 2) if i > 0 else 0)
            #errors_y.append((error_low + error_up) / 2)
            errors_y_u.append(error_up)
            errors_y_l.append(error_low)
            x_points.append(x_pt)
        plt.errorbar(x_points, efficiencies, yerr=[errors_y_l, errors_y_u],color=color_plt, marker='o', linestyle='None', label=label_plt)
        plt.xlabel("#mu_pT")
        #plt.xlabel('tau2_pt')
        plt.ylabel('eff')
        plt.title('efficiency')
        k+=1
    plt.legend(loc='lower right')

    plt.savefig(f"Studies/Triggers/Mu50_efficiencies/muPt50_allmasses.pdf")
    plt.close()



    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))


