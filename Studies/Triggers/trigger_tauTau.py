
import copy
import datetime
import os
import sys
import ROOT
import shutil
import zlib
import time

if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])


#ROOT.EnableImplicitMT(1)
ROOT.EnableThreadSafety()

import Common.Utilities as Utilities
from Analysis.HistHelper import *
from Analysis.hh_bbtautau import *
from Studies.Triggers.Trig_utilities import *

def GetTriggerWeight(trig_list):
    weights_str_list = []
    for trg in trig_list:
        if trg == 'MET' : continue
        weights_str_list.append(f'weight_tau1_TrgSF_{trg}_Central*weight_tau2_TrgSF_{trg}_Central')
    weight = " * ".join(weight_str for weight_str in weights_str_list)
    return weight


if __name__ == "__main__":
    import argparse
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=False, type=str, default='GluGluToBulkGraviton')
    parser.add_argument('--deepTauVersion', required=False, type=str, default='v2p1')
    parser.add_argument('--wantBigTau', required=False, type=bool, default=False)
    parser.add_argument('--wantTightFirstTau', required=False, type=bool, default=False)
    parser.add_argument('--wantLegend', required=False, type=bool, default=False)
    parser.add_argument('--wantNEvents', required=False, type=bool, default=False)
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
    if args.deepTauVersion=='v2p5':
        inDir = "/eos/home-k/kandroso/cms-hh-bbtautau/anaTuples/Run2_2018/v2_deepTau2p5_HTT/"
    #inDir = "/afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/output/Run2_2018/{}/anaTuples/"
    masses = []

    eff_tautau_validity = {}
    eff_tautau_application = {}
    eff_tautau_application_bigTau={}
    eff_tautau_application_bigTauOnly={}


    eff_tautau_validity_errors = {}
    eff_tautau_application_errors = {}
    eff_tautau_application_bigTau_errors={}
    eff_tautau_application_bigTauOnly_errors={}

    #### OS_ISO requirements ####

    deepTauYear = '2017' if args.deepTauVersion=='v2p1' else '2018'

    os_iso_filtering = {
        'VLoose': f"""(tau1_charge * tau2_charge < 0) && tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.VLoose.value}""",
        'Loose': f"""(tau1_charge * tau2_charge < 0) && tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Loose.value}""",
        'Medium':f"""(tau1_charge * tau2_charge < 0) && tau2_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Medium.value}"""
        }

    #### begin of loop over samples ####
    sample_list = ggR_samples if args.sample == 'GluGluToRadion' else ggBG_samples

    for sample_name in sample_list:
        tautau_labels_application = []
        tautau_labels_validity = []
        tautau_labels_application_bigTau = []
        tautau_labels_application_bigTauOnly = []

        tautau_linestyles_application = []
        tautau_linestyles_validity = []
        tautau_linestyles_application_bigTau = []
        tautau_linestyles_application_bigTauOnly = []

        tautau_colors_application = []
        tautau_colors_validity = []
        tautau_colors_application_bigTau = []
        tautau_colors_application_bigTauOnly = []

        mass_string = sample_name.split('-')[-1]
        mass_int = int(mass_string)
        masses.append(mass_int)
        inFile = os.path.join(inDir,f"{sample_name}", "nanoHTT_0.root")
        if not os.path.exists(inFile) :
            print(f"{inFile} does not exist")
            continue
        #print(sample_name)
        print(inFile)
        df_initial = ROOT.RDataFrame('Events', inFile)
        #print(f"nEventds initial : {df_initial.Count().GetValue()}" )
        dfWrapped_central  = DataFrameBuilder(ROOT.RDataFrame('Events', inFile), args.deepTauVersion)
        PrepareDfWrapped(dfWrapped_central)
        df_tauTau = dfWrapped_central.df.Filter('tauTau').Filter("tau1_gen_kind==5 && tau2_gen_kind==5")
        #nInitial_tauTau= df_tauTau.Count().GetValue()
        total_weight_string_denum = 'weight_total*weight_TauID_Central*weight_tau1_EleidSF_Central*weight_tau1_MuidSF_Central*weight_tau2_EleidSF_Central*weight_tau2_MuidSF_Central*weight_L1PreFiring_Central*weight_L1PreFiring_ECAL_Central*weight_L1PreFiring_Muon_Central'
        nInitial_tauTau = df_tauTau.Define('total_weight_string_denum', total_weight_string_denum).Sum('total_weight_string_denum').GetValue()

        df_tautau_vloose = df_tauTau.Filter(os_iso_filtering['VLoose'])
        df_tautau_loose = df_tauTau.Filter(os_iso_filtering['Loose'])
        df_tautau_medium = df_tauTau.Filter(os_iso_filtering['Medium'])
        if args.wantTightFirstTau:
            df_tautau_vloose.Filter("tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}")
            df_tautau_medium.Filter("tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}")
            df_tautau_loose.Filter("tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}")

        dataframes_channel = {
            'tautau_medium': df_tautau_medium,
            'tautau_loose': df_tautau_loose,
            'tautau_vloose': df_tautau_vloose,
        }
        dataframes_channel_bigTau = {
            'tautau_medium': df_tautau_medium,
            'tautau_loose': df_tautau_loose,
            'tautau_vloose': df_tautau_vloose,
        }

        for channel in ['tautau_medium','tautau_loose', 'tautau_vloose']:
            regions_expr= get_regions_dict('tauTau', [190,190], 150, False)
            regions_expr_bigTau= get_regions_dict('tauTau', [190,190], 150, True)
            for reg_key,reg_exp in regions_expr.items():
                dataframes_channel[channel] = dataframes_channel[channel].Define(reg_key, reg_exp)
            for reg_key,reg_exp in regions_expr_bigTau.items():
                dataframes_channel_bigTau[channel] = dataframes_channel_bigTau[channel].Define(reg_key, reg_exp)

        pass_MET_application = "(HLT_MET && MET_region && !(singleTau_region) && ! (other_trg_region))"
        pass_singleTau_application = "(HLT_singleTau && singleTau_region && !(MET_region) && ! (other_trg_region))"
        pass_other_trg_application= "({} && other_trg_region  && !(MET_region) && ! (singleTau_region) )  "
        pass_ditau_application = pass_other_trg_application.format("HLT_ditau")

        pass_singleTau_validity = f"""HLT_singleTau && ({GetTrgValidityRegion("tauTau","HLT_singleTau")})"""
        pass_MET_validity = f"""HLT_MET && ({GetTrgValidityRegion("tauTau","HLT_MET")})"""
        pass_ditau_validity = f"""HLT_ditau && ({GetTrgValidityRegion("tauTau","HLT_ditau")})"""


        ##### VALIDITY TRIGGERS, WEIGHTS #####
        trg_tauTau_list_validity_fullOR = [pass_ditau_validity,pass_singleTau_validity,pass_MET_validity]
        trg_tauTau_list_validity_diTau_singleTau = [pass_ditau_validity,pass_singleTau_validity]
        trg_tauTau_list_validity_diTau_MET = [pass_ditau_validity,pass_MET_validity]
        trg_tauTau_list_validity_diTau = [pass_ditau_validity] # ratio denumerator

        channel=f'tautau_{args.deepTauWP}'

        dataframes_channel[channel] = dataframes_channel[channel].Define(f"weight_ditau_validity",f" {pass_ditau_validity}?  weight_tau1_TrgSF_ditau_Central*weight_tau2_TrgSF_ditau_Central : 1.f ")
        dataframes_channel[channel] = dataframes_channel[channel].Define(f"weight_singleTau_validity",f"""if({pass_singleTau_validity}){{
            if (tau1_pt > 190 && abs(tau1_eta) < 2.1 && tau2_pt > 190 && abs(tau2_eta) < 2.1 ) return tau1_pt>tau2_pt ? weight_tau1_TrgSF_singleTau_Central : weight_tau2_TrgSF_singleTau_Central;
            if (tau1_pt > 190 && abs(tau1_eta) < 2.1) return weight_tau1_TrgSF_singleTau_Central;
            if (tau2_pt > 190 && abs(tau2_eta) < 2.1) return weight_tau2_TrgSF_singleTau_Central;

        }}
        return 1.f;""")

        dataframes_channel[channel] =  dataframes_channel[channel].Define("weight_singleTau_ditau_validity", "weight_ditau_validity*weight_singleTau_validity")

        ##### APPLICATION TRIGGERS AND WEIGHTS #####
        trg_tauTau_list_application_fullOR = [pass_ditau_application,pass_singleTau_application,pass_MET_application]
        trg_tauTau_list_application_diTau_singleTau = [pass_ditau_application,pass_singleTau_application]
        trg_tauTau_list_application_diTau_MET = [pass_ditau_application,pass_MET_application]
        trg_tauTau_list_application_diTau = [pass_ditau_application] # ratio denumerator

        dataframes_channel[channel] = dataframes_channel[channel].Define(f"weight_ditau_application",f" {pass_ditau_application}?  weight_tau1_TrgSF_ditau_Central*weight_tau2_TrgSF_ditau_Central : 1.f ")
        dataframes_channel[channel] = dataframes_channel[channel].Define(f"weight_singleTau_application",f"""if({pass_singleTau_application}){{
                if (tau1_pt > 190 && abs(tau1_eta) < 2.1 && tau2_pt > 190 && abs(tau2_eta) < 2.1 ) return tau1_pt>tau2_pt ? weight_tau1_TrgSF_singleTau_Central : weight_tau2_TrgSF_singleTau_Central;
                if (tau1_pt > 190 && abs(tau1_eta) < 2.1) return weight_tau1_TrgSF_singleTau_Central;
                if (tau2_pt > 190 && abs(tau2_eta) < 2.1) return weight_tau2_TrgSF_singleTau_Central;
                throw std::runtime_error("problem with region"); }}
                 return 1.f;""")
        dataframes_channel[channel] =  dataframes_channel[channel].Define("weight_singleTau_ditau_application", "weight_ditau_application*weight_singleTau_application")

        dataframes_channel_bigTau[channel] = dataframes_channel_bigTau[channel].Define(f"weight_ditau_application",f" {pass_ditau_application}?  weight_tau1_TrgSF_ditau_Central*weight_tau2_TrgSF_ditau_Central : 1.f ")
        dataframes_channel_bigTau[channel] = dataframes_channel_bigTau[channel].Define(f"weight_singleTau_application",f"""if({pass_singleTau_application}){{
                if (tau1_pt > 190 && abs(tau1_eta) < 2.1 && tau2_pt > 190 && abs(tau2_eta) < 2.1 ) return tau1_pt>tau2_pt ? weight_tau1_TrgSF_singleTau_Central : weight_tau2_TrgSF_singleTau_Central;
                if (tau1_pt > 190 && abs(tau1_eta) < 2.1) return weight_tau1_TrgSF_singleTau_Central;
                if (tau2_pt > 190 && abs(tau2_eta) < 2.1) return weight_tau2_TrgSF_singleTau_Central;
                throw std::runtime_error("problem with region"); }}
                 return 1.f;""")
        dataframes_channel_bigTau[channel] =  dataframes_channel_bigTau[channel].Define("weight_singleTau_ditau_application", "weight_ditau_application*weight_singleTau_application")

        ##### preparing dict for validity region plot #####

        AddEfficiencyToDict(dataframes_channel[channel].Define(f'total_weight_string',f'{total_weight_string_denum}*weight_ditau_validity'), trg_tauTau_list_validity_diTau, 'tauTau', nInitial_tauTau,'eff_validity_diTau', eff_tautau_validity, eff_tautau_validity_errors, args.wantNEvents)
        tautau_linestyles_validity.append('solid')
        tautau_labels_validity.append('Legacy')
        tautau_colors_validity.append("red")

        AddEfficiencyToDict(dataframes_channel[channel].Define(f'total_weight_string',f'{total_weight_string_denum}*weight_singleTau_ditau_validity'), trg_tauTau_list_validity_diTau_singleTau, 'tauTau', nInitial_tauTau,'eff_validity_diTauSingleTau', eff_tautau_validity, eff_tautau_validity_errors, args.wantNEvents)
        tautau_linestyles_validity.append('solid')
        tautau_labels_validity.append("Legacy||singleTau")
        tautau_colors_validity.append("blue")
        #tautau_colors_validity.append("gold")

        ##### preparing dict for application region plot #####

        AddEfficiencyToDict(dataframes_channel[channel].Define(f'total_weight_string',f'{total_weight_string_denum}*weight_ditau_application'), trg_tauTau_list_application_diTau, 'tauTau', nInitial_tauTau,'eff_application_diTau', eff_tautau_application, eff_tautau_application_errors, args.wantNEvents)
        tautau_linestyles_application.append('solid')
        tautau_labels_application.append('Legacy')
        tautau_colors_application.append("red")

        AddEfficiencyToDict(dataframes_channel[channel].Define(f'total_weight_string',f'{total_weight_string_denum}*weight_singleTau_ditau_application'), trg_tauTau_list_application_diTau_singleTau, 'tauTau', nInitial_tauTau,'eff_application_diTauSingleTau', eff_tautau_application, eff_tautau_application_errors, args.wantNEvents)
        tautau_linestyles_application.append('solid')
        tautau_labels_application.append("Legacy||singleTau")
        tautau_colors_application.append("blue")

        AddEfficiencyToDict(dataframes_channel[channel].Define(f'total_weight_string',f'{total_weight_string_denum}*weight_singleTau_ditau_validity'), trg_tauTau_list_validity_diTau_singleTau, 'tauTau', nInitial_tauTau,'eff_validity_diTau_singleTau', eff_tautau_application, eff_tautau_application_errors, args.wantNEvents)
        tautau_labels_application.append("Legacy||singleTau")
        tautau_colors_application.append("dodgerblue")
        tautau_linestyles_application.append('dashed')

        ##### preparing dict for application region plot bigTau = False VS bigTau = True #####


        #AddEfficiencyToDict(dataframes_channel[channel].Define(f'total_weight_string',f'{total_weight_string_denum}*{weight_diTau}'), trg_tauTau_list_application_diTau, 'tauTau', nInitial_tauTau,'eff_application_diTau', eff_tautau_application_bigTau, eff_tautau_application_bigTau_errors, args.wantNEvents)
        #tautau_linestyles_application_bigTau.append('solid')
        #tautau_labels_application_bigTau.append('Legacy')
        #tautau_colors_application_bigTau.append("red")

        AddEfficiencyToDict(dataframes_channel[channel].Define(f'total_weight_string',f'{total_weight_string_denum}*weight_singleTau_ditau_application'), trg_tauTau_list_application_diTau_singleTau, 'tauTau', nInitial_tauTau,'eff_application_diTau_singleTau', eff_tautau_application_bigTau, eff_tautau_application_bigTau_errors, args.wantNEvents)
        tautau_linestyles_application_bigTau.append('solid')
        tautau_labels_application_bigTau.append("Legacy||singleTau - !(bigTau)")
        tautau_colors_application_bigTau.append("cyan")



        #### bigTau = True ####


        #AddEfficiencyToDict(dataframes_channel_bigTau[channel].Define(f'total_weight_string',f'{total_weight_string_denum}*{weight_diTau}'), trg_tauTau_list_application_diTau, 'tauTau', nInitial_tauTau,'eff_application_diTau_bigTau', eff_tautau_application_bigTau, eff_tautau_application_bigTau_errors, args.wantNEvents)
        #tautau_linestyles_application_bigTau.append('solid')
        #tautau_labels_application_bigTau.append('Legacy')
        #tautau_colors_application_bigTau.append("maroon")

        AddEfficiencyToDict(dataframes_channel_bigTau[channel].Define(f'total_weight_string',f'{total_weight_string_denum}*weight_singleTau_ditau_application'), trg_tauTau_list_application_diTau_singleTau, 'tauTau', nInitial_tauTau,'eff_application_diTau_singleTau_bigTau', eff_tautau_application_bigTau, eff_tautau_application_bigTau_errors, args.wantNEvents)
        tautau_linestyles_application_bigTau.append('solid')
        tautau_labels_application_bigTau.append("Legacy||singleTau - bigTau")
        tautau_colors_application_bigTau.append("midnightblue")
        '''
        AddEfficiencyToDict(dataframes_channel[channel].Define(f'total_weight_string',f'{total_weight_string_denum}*weight_singleTau_ditau_validity'), trg_tauTau_list_validity_diTau_singleTau, 'tauTau', nInitial_tauTau,'eff_validity_diTau_singleTau', eff_tautau_application_bigTau, eff_tautau_application_bigTau_errors, args.wantNEvents)
        tautau_labels_application_bigTau.append("Legacy||singleTau - Validity")
        tautau_colors_application_bigTau.append("dodgerblue")
        tautau_linestyles_application_bigTau.append('dashed')
        '''


        ############# bigTau only #############

        AddEfficiencyToDict(dataframes_channel_bigTau[channel].Define(f'total_weight_string',f'{total_weight_string_denum}*weight_ditau_application'), trg_tauTau_list_application_diTau, 'tauTau', nInitial_tauTau,'eff_application_diTau', eff_tautau_application_bigTauOnly, eff_tautau_application_bigTauOnly_errors, args.wantNEvents)
        tautau_linestyles_application_bigTauOnly.append('solid')
        tautau_labels_application_bigTauOnly.append('Legacy')
        tautau_colors_application_bigTauOnly.append("maroon")

        AddEfficiencyToDict(dataframes_channel_bigTau[channel].Define(f'total_weight_string',f'{total_weight_string_denum}*weight_singleTau_ditau_application'), trg_tauTau_list_application_diTau_singleTau, 'tauTau', nInitial_tauTau,'eff_application_diTauSingleTau', eff_tautau_application_bigTauOnly, eff_tautau_application_bigTauOnly_errors, args.wantNEvents)
        tautau_linestyles_application_bigTauOnly.append('solid')
        tautau_labels_application_bigTauOnly.append("Legacy||singleTau")
        tautau_colors_application_bigTauOnly.append("midnightblue")

        AddEfficiencyToDict(dataframes_channel_bigTau[channel].Define(f'total_weight_string',f'{total_weight_string_denum}*weight_singleTau_ditau_application'), trg_tauTau_list_validity_diTau_singleTau, 'tauTau', nInitial_tauTau,'eff_validity_diTau_singleTau', eff_tautau_application_bigTauOnly, eff_tautau_application_bigTauOnly_errors, args.wantNEvents)
        tautau_labels_application_bigTauOnly.append("Legacy||singleTau")
        tautau_colors_application_bigTauOnly.append("dodgerblue")
        tautau_linestyles_application_bigTauOnly.append('dashed')




    x_values = masses
    ratio_ref = eff_tautau_validity['eff_validity_diTau']
    ratio_ref_errors = eff_tautau_validity_errors['eff_validity_diTau_error']

    makeplot(eff_tautau_validity,eff_tautau_validity_errors, tautau_labels_validity, tautau_linestyles_validity, 'tauTau', x_values, args.sample, args.deepTauWP, args.deepTauVersion, ratio_ref,ratio_ref_errors,'_validity_noMET',tautau_colors_validity, args.wantNEvents, args.wantLegend)

    makeplot(eff_tautau_application,eff_tautau_application_errors, tautau_labels_application, tautau_linestyles_application, 'tauTau', x_values, args.sample, args.deepTauWP, args.deepTauVersion, ratio_ref,ratio_ref_errors,'_application_noMET',tautau_colors_application, args.wantNEvents, args.wantLegend)

    makeplot(eff_tautau_application_bigTau,eff_tautau_application_bigTau_errors, tautau_labels_application_bigTau, tautau_linestyles_application_bigTau, 'tauTau', x_values, args.sample,args.deepTauWP, args.deepTauVersion, ratio_ref,ratio_ref_errors,'_application_bigTau_noMET',tautau_colors_application_bigTau, args.wantNEvents, args.wantLegend, True)

    makeplot(eff_tautau_application_bigTauOnly,eff_tautau_application_bigTauOnly_errors, tautau_labels_application_bigTauOnly, tautau_linestyles_application_bigTauOnly, 'tauTau', x_values, args.sample,args.deepTauWP, args.deepTauVersion, ratio_ref,ratio_ref_errors,'_application_bigTauOnly_noMET',tautau_colors_application_bigTauOnly, args.wantNEvents, args.wantLegend)

    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))


