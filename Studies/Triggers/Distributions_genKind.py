
import copy
import datetime
import os
import sys
import ROOT
import shutil
import zlib
import time
import json
import matplotlib.pyplot as plt

if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])


#ROOT.EnableImplicitMT(1)
ROOT.EnableThreadSafety()

from Analysis.HistHelper import *
from Analysis.hh_bbtautau import *
from Studies.Triggers.Trig_utilities import *

import Common.BaselineSelection as Baseline
import Corrections.Corrections as Corrections
import Common.Utilities as Utilities
from Corrections.CorrectionsCore import *
from Corrections.tau import *
from Studies.Triggers.GetSOverSqrtB_nonreordered import *


def getHistogram(df, trg_reg_list,weight_expression, model, var):
    filter_expr = f' ('
    filter_expr+= ' || '.join(trg for trg in trg_reg_list)
    filter_expr+= ')'
    histo = df.Define("weight_for_histo",weight_expression).Filter(filter_expr).Histo1D(model, var, "weight_for_histo")
    return histo


if __name__ == "__main__":
    import argparse
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=False, type=str, default='GluGluToBulkGraviton')
    parser.add_argument('--deepTauVersion', required=False, type=str, default='v2p1')
    parser.add_argument('--want2b', required=False, type=bool, default=False)
    parser.add_argument('--bckg_sel', required=False, type=str, default="0")
    parser.add_argument('--bckg_samples', required=False, type=str, default="")



    args = parser.parse_args()

    startTime = time.time()
    headers_dir = os.path.dirname(os.path.abspath(__file__))
    ROOT.gROOT.ProcessLine(f".include {os.environ['ANALYSIS_PATH']}")
    ROOT.gInterpreter.Declare(f'#include "include/KinFitNamespace.h"')
    ROOT.gInterpreter.Declare(f'#include "include/HistHelper.h"')
    ROOT.gInterpreter.Declare(f'#include "include/Utilities.h"')
    ROOT.gROOT.ProcessLine('#include "include/AnalysisTools.h"')
    header_path_RootExt = os.path.join(os.environ['ANALYSIS_PATH'], "PlotKit/include/RootExt.h")
    ROOT.gInterpreter.Declare(f'#include "{header_path_RootExt}"')

    #### useful stuff ####
    config_file = "/afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/config/samples_Run2_2018.yaml"
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)


    inDir = "/eos/home-k/kandroso/cms-hh-bbtautau/anaTuples/Run2_2018/v2_deepTau2p1_HTT/"

    Corrections.Initialize(config=config['GLOBAL'], isData=False, load_pu=False, load_tau=True, load_trg=False, load_btag=False, loadBTagEff=False, load_met=False, load_mu = False, load_ele=False, load_puJetID=False, load_jet=False)

    MET_SF_fileName = os.path.join("/afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/Studies/Triggers/files", "eff_Data_Mu_MC_DY_TT_WJets_mumu_metnomu_et_TRG_METNoMu120_CUTS_mhtnomu_et_L_100_default_mutau.root")
    #MET_SF_file = ROOT.TFile.Open(MET_SF_fileName, "READ")
    ROOT.LoadSigmoidFunction(MET_SF_fileName, "SigmoidFuncSF")


    x_bins = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5 , 6.5, 7.5]
    x_bins_vec = Utilities.ListToVector(x_bins, "double")

    #### OS_ISO requirements ####

    deepTauYear = '2017' if args.deepTauVersion=='v2p1' else '2018'

    os_iso_filtering = {
        'VVLoose': f"""(tau1_charge * tau2_charge < 0) && tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VVLoose.value}""",
        'VLoose': f"""(tau1_charge * tau2_charge < 0) && tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VLoose.value}""",
        'Loose': f"""(tau1_charge * tau2_charge < 0) && tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Loose.value}""",
        'Medium':f"""(tau1_charge * tau2_charge < 0) && tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Medium.value}""",
        'Tight':f"""(tau1_charge * tau2_charge < 0) && tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}""",
        'VTight':f"""(tau1_charge * tau2_charge < 0) && tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VTight.value}"""
        }
    #### begin of loop over bckg samples ####
    bckg_samples = args.bckg_samples.split(",") if args.bckg_samples != "" else []
    #print(bckg_samples)
    #bckg_samples = ['TTToHadronic','TTTo2L2Nu']
    #bckg_samples = ['TTToSemiLeptonic']
    #bckg_samples = ['TTToSemiLeptonic','TTToHadronic','TTTo2L2Nu']

    signal_samples = ggR_samples if args.sample == 'GluGluToRadion' else ggBG_samples
    if bckg_samples:
        signal_samples = []
    #signal_samples = []
    b_values = {}
    s_values = {}

    #### begin of loop over signal samples ####

    background_string = {
        "0":"tau1_gen_kind!=5 || tau2_gen_kind!z=5",
        "1":"tau1_gen_kind == 5 && tau2_gen_kind == 5",
        "2":"tau1_gen_kind != 5 && tau2_gen_kind != 5"
        }
    for sample_name in signal_samples + bckg_samples :
        print(sample_name)
        all_inFiles = getInFiles(os.path.join(inDir,f"{sample_name}"))
        print(all_inFiles)
        df_initial = ROOT.RDataFrame('Events', Utilities.ListToVector(all_inFiles))
        dfWrapped_central  = DataFrameBuilder(df_initial, args.deepTauVersion)
        PrepareDfWrapped(dfWrapped_central)

        df_tauTau = dfWrapped_central.df.Filter('tauTau')

        #df_tautau_VVLoose = df_tauTau.Filter(os_iso_filtering['VVLoose'])
        #df_tautau_VLoose = df_tauTau.Filter(os_iso_filtering['VLoose'])
        df_tautau_Loose = df_tauTau.Filter(os_iso_filtering['Loose'])
        df_tautau_Medium = df_tauTau.Filter(os_iso_filtering['Medium'])
        df_tautau_Tight = df_tauTau.Filter(os_iso_filtering['Tight'])
        df_tautau_VTight = df_tauTau.Filter(os_iso_filtering['VTight'])
        dataframes_channel = {
            'tautau_Medium_Medium': df_tautau_Medium,
            'tautau_Medium_Loose': df_tautau_Loose,
            #'tautau_Medium_VLoose': df_tautau_VLoose,
            #'tautau_Medium_VVLoose': df_tautau_VVLoose,
            'tautau_Medium_Tight': df_tautau_Tight,
            'tautau_Medium_VTight': df_tautau_VTight,

            'tautau_Tight_Medium': df_tautau_Medium.Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}"),
            'tautau_Tight_Tight': df_tautau_Tight.Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}"),
            'tautau_Tight_Medium': df_tautau_Medium.Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}"),
            'tautau_Tight_Loose': df_tautau_Loose.Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}"),
            #'tautau_Tight_VLoose': df_tautau_VLoose.Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}"),
            #'tautau_Tight_VVLoose': df_tautau_VVLoose.Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}"),

            'tautau_VTight_VTight': df_tautau_VTight.Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VTight.value}"),
            'tautau_VTight_Medium': df_tautau_Medium.Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VTight.value}"),
            'tautau_VTight_Loose': df_tautau_Loose.Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VTight.value}"),
            #'tautau_VTight_VLoose': df_tautau_VLoose.Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VTight.value}"),
            #'tautau_VTight_VVLoose': df_tautau_VVLoose.Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VTight.value}"),

        }
        ##### application triggers #####

        pass_MET_application = "(HLT_MET && MET_region && !(singleTau_region) && ! (other_trg_region))"
        pass_singleTau_application = "(HLT_singleTau && singleTau_region && !(MET_region) && ! (other_trg_region))"
        pass_other_trg_application= "({} && other_trg_region  && !(MET_region) && ! (singleTau_region) )  "
        pass_ditau_application = pass_other_trg_application.format("HLT_ditau")
        trg_tauTau_list_application_fullOR = [pass_ditau_application,pass_singleTau_application,pass_MET_application]
        trg_tauTau_list_application_ditau_singleTau = [pass_ditau_application,pass_singleTau_application]
        trg_tauTau_list_application_ditau_MET = [pass_ditau_application,pass_MET_application]
        trg_tauTau_list_application_ditau = [pass_ditau_application] # ratio denumerator

        for channel in dataframes_channel.keys():


            wp_tau1 = channel.split('_')[1]
            wp_tau2 = channel.split('_')[2]
            #dataframes_channel[channel] = dataframes_channel[channel].Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSe >= {Utilities.WorkingPointsTauVSe.Loose.value}")
            #dataframes_channel[channel] = dataframes_channel[channel].Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSe >= {Utilities.WorkingPointsTauVSe.Medium.value}")
            #dataframes_channel[channel] = dataframes_channel[channel].Filter(f"tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSe >= {Utilities.WorkingPointsTauVSe.Loose.value}")
            #dataframes_channel[channel] = dataframes_channel[channel].Filter(f"tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSe >= {Utilities.WorkingPointsTauVSe.Medium.value}")
            dataframes_channel[channel] = getTauWPWeight(dataframes_channel[channel], wp_tau1, wp_tau2)
            dataframes_channel[channel] = dataframes_channel[channel].Define("weights_0","weight_total*weight_TauID_new*weight_tau1_EleidSF_Central*weight_tau1_MuidSF_Central*weight_tau2_EleidSF_Central*weight_tau2_MuidSF_Central*weight_L1PreFiring_Central*weight_L1PreFiring_ECAL_Central*weight_L1PreFiring_Muon_Central").Define("sel_id", "int(tau1_gen_kind == 5) + int(tau2_gen_kind == 5)")

            r_factor = getRFactor(dataframes_channel[channel]) if args.want2b else 1
            #print(f"rfactor for channel {channel} sample {sample_name} is {r_factor}")
            if args.want2b:
                dataframes_channel[channel] = dataframes_channel[channel].Filter("nSelBtag>1")
            #if sample_name in signal_samples:
            #    dataframes_channel[channel] = dataframes_channel[channel].Filter("sel_id == 2")
            #else:
            dataframes_channel[channel] = dataframes_channel[channel].Filter(f"sel_id == {args.bckg_sel}")
            btag_string = f"{r_factor}*weight_bTagShapeSF" if args.want2b else "1"

            total_weight_string = f'{btag_string}*weights_0'
            dataframes_channel[channel] = dataframes_channel[channel].Define('final_weight',total_weight_string)


            regions_expr= get_regions_dict('tauTau', [190,190], 180, False)
            for reg_key,reg_exp in regions_expr.items():
                dataframes_channel[channel] = dataframes_channel[channel].Define(reg_key, reg_exp)
            dataframes_channel[channel] = dataframes_channel[channel].Define(f"weight_ditau_application",f" {pass_ditau_application}?  weight_tau1_TrgSF_ditau_Central*weight_tau2_TrgSF_ditau_Central : 1.f ")
            dataframes_channel[channel] = dataframes_channel[channel].Define(f"weight_singleTau_application",f"""if({pass_singleTau_application}){{
                if (tau1_pt > 190 && abs(tau1_eta) < 2.1 && tau2_pt > 190 && abs(tau2_eta) < 2.1 ) return tau1_pt>tau2_pt ? weight_tau1_TrgSF_singleTau_Central : weight_tau2_TrgSF_singleTau_Central;
                if (tau1_pt > 190 && abs(tau1_eta) < 2.1) return weight_tau1_TrgSF_singleTau_Central;
                if (tau2_pt > 190 && abs(tau2_eta) < 2.1) return weight_tau2_TrgSF_singleTau_Central;
                throw std::runtime_error("problem with region"); }}
                return 1.f;""")
            dataframes_channel[channel] = dataframes_channel[channel].Define("weight_MET_application", f"{pass_MET_application}? getMETsf(metnomu_pt) : 1.f")


            model_1 = ROOT.RDF.TH1DModel("tau1_gen_kind", "tau1_gen_kind", x_bins_vec.size()-1, x_bins_vec.data())
            model_2 = ROOT.RDF.TH1DModel("tau2_gen_kind", "tau2_gen_kind", x_bins_vec.size()-1, x_bins_vec.data())

            histo_genKind_1 = getHistogram(dataframes_channel[channel], trg_tauTau_list_application_fullOR,f'final_weight*weight_ditau_application*weight_singleTau_application*weight_MET_application', model_1, "tau1_gen_kind")
            histo_genKind_2 = getHistogram(dataframes_channel[channel], trg_tauTau_list_application_fullOR,f'final_weight*weight_ditau_application*weight_singleTau_application*weight_MET_application', model_2, "tau2_gen_kind")
            print(channel)
            print(f"nEvents are {histo_genKind_1.GetEntries()}")
            outFile = ROOT.TFile.Open(f"Studies/Triggers/files/genkind_distrib_{args.bckg_sel}/gen_kinds_{sample_name}_{channel}_VSeVVLoose.root", "RECREATE")
            histo_genKind_1.Write()
            histo_genKind_2.Write()
            outFile.Close()

    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))


