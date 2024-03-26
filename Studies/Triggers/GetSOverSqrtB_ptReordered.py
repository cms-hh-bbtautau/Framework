
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

ROOT.gInterpreter.Declare("""
    std::unique_ptr<TF1> sigmoid_function;

    void LoadSigmoidFunction(const char* fileName, const char* functionName) {
        TFile file(fileName, "READ");
        sigmoid_function.reset((TF1*)file.Get(functionName)->Clone());
    }

    float getMETsf( const float metnomu_pt)
    {
        if(metnomu_pt>180 && metnomu_pt<300){
            if (sigmoid_function){
                return sigmoid_function->Eval(metnomu_pt);
            }
            else {
                std::cerr << "Sigmoid function is not loaded." << std::endl;
                return 0.0;
            }
        }
    return 1.f;
    }
    """)


def getTauWPWeight(df, working_point_tau1, working_point_tau2):
    WPs = ["VVLoose","VLoose", "Loose","Medium","Tight","VTight"]
    df = df.Define(f"weight_tau1new_TauID_SF_{working_point_tau1}",
                            f'''tau1new_gen_kind == 5?
                            ::correction::TauCorrProvider::getGlobal().getSF_WPStrings(tau1new_pt,  tau1new_eta, tau1new_decayMode,tau1new_gen_kind, "{working_point_tau1}") : 1.;
                            ''')
    df = df.Define(f"weight_tau2new_TauID_SF_{working_point_tau2}",
                            f'''tau2new_gen_kind == 5?
                            ::correction::TauCorrProvider::getGlobal().getSF_WPStrings(tau2new_pt,  tau2new_eta, tau2new_decayMode,tau2new_gen_kind,"{working_point_tau2}") : 1.;
                            ''')
    df = df.Define(f"weight_TauID_new",f"""weight_tau1new_TauID_SF_{working_point_tau1}*weight_tau2new_TauID_SF_{working_point_tau2}""")

    return df


def getSumWeights(df, trg_reg_list,weight_expression):
    filter_expr = f' ('
    filter_expr+= ' || '.join(trg for trg in trg_reg_list)
    filter_expr+= ')'
    sum_channel = df.Define("weight_for_sum",weight_expression).Filter(filter_expr).Sum("weight_for_sum").GetValue()
    sum_error = math.sqrt(df.Define("weight_for_sum_error",f'{weight_expression}*{weight_expression}').Filter(filter_expr).Sum("weight_for_sum_error").GetValue())
    return sum_channel,sum_error

def getInFiles(inDir):
    infiles= []
    for rootfile in os.listdir(inDir):
        if rootfile.endswith(".root"):
            infiles.append(os.path.join(inDir, rootfile))
    return infiles

def getRFactor(df):
    num = df.Sum("weights_0").GetValue()
    den = df.Define("weight_0_btag","weights_0*weight_bTagShapeSF").Sum("weight_0_btag").GetValue()
    return num/den


# ["Electron_mvaIso_WP90","Electron_mvaNoIso_WP90","Electron_pfRelIso03_all","HasMatching_MET","HasMatching_ditau","HasMatching_etau","HasMatching_mutau","HasMatching_singleEle","HasMatching_singleMu","HasMatching_singleMu50","HasMatching_singleTau","Muon_pfRelIso04_all","Muon_tkRelIso","charge","decayMode","eta","gen_charge","gen_kind","gen_nChHad","gen_nNeutHad","gen_vis_eta","gen_vis_mass","gen_vis_phi","gen_vis_pt","idDeepTau2017v2p1VSe","idDeepTau2017v2p1VSjet","idDeepTau2017v2p1VSmu","idDeepTau2018v2p5VSe","idDeepTau2018v2p5VSjet","idDeepTau2018v2p5VSmu","iso","legType","mass","phi","pt","rawDeepTau2017v2p1VSe","rawDeepTau2017v2p1VSjet","rawDeepTau2017v2p1VSmu","rawDeepTau2018v2p5VSe","rawDeepTau2018v2p5VSjet","rawDeepTau2018v2p5VSmu","seedingJet_eta","seedingJet_hadronFlavour","seedingJet_mass","seedingJet_partonFlavour","seedingJet_phi","seedingJet_pt"]
def reorderInPt(df, observables_list = ["charge","decayMode","eta","gen_charge","gen_kind","gen_nChHad","gen_nNeutHad","gen_vis_eta","gen_vis_mass","gen_vis_phi","gen_vis_pt","idDeepTau2017v2p1VSe","idDeepTau2017v2p1VSjet","idDeepTau2017v2p1VSmu","idDeepTau2018v2p5VSe","idDeepTau2018v2p5VSjet","idDeepTau2018v2p5VSmu","iso","legType","mass","phi","pt"]):
    for var in observables_list:
        df = df.Define(f'tau1new_{var}', f'tau1_pt > tau2_pt ? tau1_{var} : tau2_{var}')
        df = df.Define(f'tau2new_{var}', f'tau1_pt > tau2_pt ? tau2_{var} : tau1_{var}')
    return df

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


    inDir = "/eos/home-k/kandroso/cms-hh-bbtautau/anaTuples/Run2_2018/v3_deepTau2p1_HTT/"

    Corrections.Initialize(config=config['GLOBAL'], isData=False, load_pu=False, load_tau=True, load_trg=False, load_btag=False, loadBTagEff=False, load_met=False, load_mu = False, load_ele=False, load_puJetID=False, load_jet=False)

    MET_SF_fileName = os.path.join("/afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/Studies/Triggers/files", "eff_Data_Mu_MC_DY_TT_WJets_mumu_metnomu_et_TRG_METNoMu120_CUTS_mhtnomu_et_L_100_default_mutau.root")
    #MET_SF_file = ROOT.TFile.Open(MET_SF_fileName, "READ")
    ROOT.LoadSigmoidFunction(MET_SF_fileName, "SigmoidFuncSF")


    #### OS_ISO requirements ####

    deepTauYear = '2017' if args.deepTauVersion=='v2p1' else '2018'

    os_iso_filtering = {
        'VVLoose': f"""(tau1new_charge * tau2new_charge < 0) && tau2new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VVLoose.value}""",
        'VLoose': f"""(tau1new_charge * tau2new_charge < 0) && tau2new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VLoose.value}""",
        'Loose': f"""(tau1new_charge * tau2new_charge < 0) && tau2new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Loose.value}""",
        'Medium':f"""(tau1new_charge * tau2new_charge < 0) && tau2new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Medium.value}""",
        'Tight':f"""(tau1new_charge * tau2new_charge < 0) && tau2new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}""",
        'VTight':f"""(tau1new_charge * tau2new_charge < 0) && tau2new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VTight.value}"""
        }
    #### begin of loop over bckg samples ####
    bckg_samples = args.bckg_samples.split(",") if args.bckg_samples != "" else []
    ##print(bckg_samples)
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
        "0":"tau1new_gen_kind!=5 || tau2new_gen_kind!=5",
        "1":"tau1new_gen_kind == 5 && tau2new_gen_kind == 5",
        "2":"tau1new_gen_kind != 5 && tau2new_gen_kind != 5"
        }
    for sample_name in signal_samples + bckg_samples :
        #print(sample_name)
        all_inFiles = getInFiles(os.path.join(inDir,f"{sample_name}"))
        ##print(all_inFiles)
        df_initial = ROOT.RDataFrame('Events', Utilities.ListToVector(all_inFiles))
        df_initial = reorderInPt(df_initial)
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
            #'tautau_Medium_VLoose': df_tautau_VLoose,
            #'tautau_Medium_VVLoose': df_tautau_VVLoose,
            'tautau_Medium_VTight': df_tautau_VTight.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Medium.value}"),
            'tautau_Medium_Tight': df_tautau_Tight.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Medium.value}"),
            'tautau_Medium_Medium': df_tautau_Medium.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Medium.value}"),
            'tautau_Medium_Loose': df_tautau_Loose.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Medium.value}"),


            'tautau_Tight_VTight': df_tautau_VTight.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}"),
            'tautau_Tight_Tight': df_tautau_Tight.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}"),
            'tautau_Tight_Medium': df_tautau_Medium.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}"),
            'tautau_Tight_Loose': df_tautau_Loose.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}"),
            #'tautau_Tight_VLoose': df_tautau_VLoose.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}"),
            #'tautau_Tight_VVLoose': df_tautau_VVLoose.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}"),

            'tautau_VTight_VTight': df_tautau_VTight.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VTight.value}"),
            'tautau_VTight_Tight': df_tautau_Tight.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VTight.value}"),
            'tautau_VTight_Medium': df_tautau_Medium.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VTight.value}"),
            'tautau_VTight_Loose': df_tautau_Loose.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VTight.value}"),
            #'tautau_VTight_VLoose': df_tautau_VLoose.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VTight.value}"),
            #'tautau_VTight_VVLoose': df_tautau_VVLoose.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VTight.value}"),

            'tautau_Loose_VTight': df_tautau_VTight.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Loose.value}"),
            'tautau_Loose_Tight': df_tautau_Tight.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Loose.value}"),
            'tautau_Loose_Medium': df_tautau_Medium.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Loose.value}"),
            'tautau_Loose_Loose': df_tautau_Loose.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Loose.value}"),
            #'tautau_Loose_VLoose': df_tautau_VLoose.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Loose.value}"),
            #'tautau_Loose_VVLoose': df_tautau_VVLoose.Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Loose.value}"),


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
            #print(channel)
            #print(dataframes_channel[channel].Count().GetValue())
            dataframes_channel[channel] = dataframes_channel[channel].Filter(f"tau1new_idDeepTau{deepTauYear}{args.deepTauVersion}VSe >= {Utilities.WorkingPointsTauVSe.Loose.value}")
            #print(dataframes_channel[channel].Count().GetValue())
            #dataframes_channel[channel] = dataframes_channel[channel].Filter(f"tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSe >= {Utilities.WorkingPointsTauVSe.Medium.value}")
            dataframes_channel[channel] = dataframes_channel[channel].Filter(f"tau2new_idDeepTau{deepTauYear}{args.deepTauVersion}VSe >= {Utilities.WorkingPointsTauVSe.Loose.value}")
            #print(dataframes_channel[channel].Count().GetValue())
            #dataframes_channel[channel] = dataframes_channel[channel].Filter(f"tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSe >= {Utilities.WorkingPointsTauVSe.Medium.value}")
            dataframes_channel[channel] = getTauWPWeight(dataframes_channel[channel], wp_tau1, wp_tau2)
            dataframes_channel[channel] = dataframes_channel[channel].Define("weights_0","weight_total*weight_TauID_new*weight_tau1_EleidSF_Central*weight_tau1_MuidSF_Central*weight_tau2_EleidSF_Central*weight_tau2_MuidSF_Central*weight_L1PreFiring_Central*weight_L1PreFiring_ECAL_Central*weight_L1PreFiring_Muon_Central").Define("sel_id", "int(tau1new_gen_kind == 5) + int(tau2new_gen_kind == 5)")
            #print(dataframes_channel[channel].Sum('weights_0').GetValue())

            r_factor = getRFactor(dataframes_channel[channel]) if args.want2b else 1
            #print(f"rfactor for channel {channel} sample {sample_name} is {r_factor}")
            if args.want2b:
                dataframes_channel[channel] = dataframes_channel[channel].Filter("nSelBtag>1")
            if sample_name in signal_samples:
                dataframes_channel[channel] = dataframes_channel[channel].Filter("sel_id == 2")
            else:
                dataframes_channel[channel] = dataframes_channel[channel].Filter(f"sel_id == {args.bckg_sel}")
            btag_string = f"{r_factor}*weight_bTagShapeSF" if args.want2b else "1"

            total_weight_string = f'{btag_string}*weights_0'
            dataframes_channel[channel] = dataframes_channel[channel].Define('final_weight',total_weight_string)
            #print(dataframes_channel[channel].Sum('weights_0').GetValue())
            #print(dataframes_channel[channel].Sum('final_weight').GetValue())
            regions_expr= get_regions_dict_ptreordered('tauTau', [190,190], 180, False)
            for reg_key,reg_exp in regions_expr.items():
                dataframes_channel[channel] = dataframes_channel[channel].Define(reg_key, reg_exp)
            dataframes_channel[channel] = dataframes_channel[channel].Define(f"weight_ditau_application",f" {pass_ditau_application}?  weight_tau1_TrgSF_ditau_Central*weight_tau2_TrgSF_ditau_Central : 1.f ")
            dataframes_channel[channel] = dataframes_channel[channel].Define(f"weight_singleTau_application",f"""if({pass_singleTau_application}){{
                if (tau1_pt > 190 && abs(tau1_eta) < 2.1 && tau2_pt > 190 && abs(tau2_eta) < 2.1 ) return tau1_pt>tau2_pt ? weight_tau1_TrgSF_singleTau_Central : weight_tau2_TrgSF_singleTau_Central;
                else if (tau1new_pt > 190 && abs(tau1new_eta) < 2.1) return tau1new_pt == tau1_pt ? weight_tau1_TrgSF_singleTau_Central : weight_tau2_TrgSF_singleTau_Central ;
                else if (tau2new_pt > 190 && abs(tau2new_eta) < 2.1) return tau2new_pt == tau2_pt ? weight_tau2_TrgSF_singleTau_Central : weight_tau1_TrgSF_singleTau_Central ;
                }}
                return 1.f;""")
            dataframes_channel[channel] = dataframes_channel[channel].Define("weight_MET_application", f"{pass_MET_application}? getMETsf(metnomu_pt) : 1.f")

            sum_weights, sum_weights_error = getSumWeights(dataframes_channel[channel], trg_tauTau_list_application_fullOR,f'weights_0*weight_ditau_application*weight_singleTau_application*weight_MET_application')
            #print(dataframes_channel[channel].Define("provaweight", 'weights_0*weight_ditau_application*weight_singleTau_application*weight_MET_application').Sum("provaweight").GetValue())
            #print()

            if channel not in b_values.keys():
                b_values[channel] = {}
                b_values[channel]["all"] = 0

            if channel not in s_values.keys():
                s_values[channel] = {}
                s_values[channel]["all"] = 0
            if sample_name in signal_samples and sample_name not in s_values[channel].keys():
                s_values[channel][sample_name] = {}

            if sample_name in bckg_samples:
                if sample_name not in b_values[channel].keys():
                    b_values[channel][sample_name] = {}
                b_values[channel][sample_name]['sum'] = sum_weights
                b_values[channel][sample_name]['error'] = sum_weights_error
                b_values[channel]['all'] += sum_weights
            else:
                if sample_name not in s_values[channel].keys():
                    s_values[channel][sample_name] = {}
                s_values[channel][sample_name]['sum'] = sum_weights
                s_values[channel][sample_name]['error'] = sum_weights_error
                s_values[channel]['all'] += sum_weights

    ##print(f"b_values = {b_values}")
    sampleName = 'bulk_graviton'
    if args.sample == 'GluGluToRadion':
        sampleName='radion'
    suffix = '_2btag' if args.want2b == True else ''
    bckg_names = ""
    if bckg_samples:
        bckg_names = "_".join(b for b in bckg_samples)
        if 'DY' in bckg_samples[0]:
            bckg_names = 'DY'
        if 'TT' in bckg_samples[0]:
            bckg_names = 'TT'

    with open(f'Studies/Triggers/files/pt_reordered/signals_{sampleName}{suffix}_ptReordered.json', 'w') as f:
        json.dump(s_values, f, indent=4)
    if bckg_samples:
        with open(f'Studies/Triggers/files/pt_reordered/backgrounds_{bckg_names}_{args.bckg_sel}{suffix}_ptReordered.json', 'w') as f:
            json.dump(b_values, f, indent=4)

    ##print(f"s_values={s_values}")
    executionTime = (time.time() - startTime)
    #print('Execution time in seconds: ' + str(executionTime))


