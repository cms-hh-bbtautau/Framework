
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
from Analysis.HistHelper import *
from Analysis.hh_bbtautau import *
from Studies.Triggers.Trig_utilities import *

if __name__ == "__main__":
    import argparse
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=False, type=str, default='GluGluToBulkGraviton')
    parser.add_argument('--deepTauVersion', required=False, type=str, default='v2p1')
    parser.add_argument('--wantBigTau', required=False, type=bool, default=False)
    parser.add_argument('--wantNEvents', required=False, type=bool, default=False)
    parser.add_argument('--wantLegend', required=False, type=bool, default=False)
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

    eff_etau_validity = {}
    eff_etau_application = {}
    eff_etau_application_bigTau={}
    eff_etau_application_bigTauOnly={}

    eff_etau_validity_errors = {}
    eff_etau_application_errors = {}
    eff_etau_application_bigTau_errors={}
    eff_etau_application_bigTauOnly_errors={}

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
        etau_labels_application = []
        etau_labels_validity = []
        etau_labels_application_bigTau = []
        etau_labels_application_bigTauOnly = []

        etau_linestyles_application = []
        etau_linestyles_validity = []
        etau_linestyles_application_bigTau = []
        etau_linestyles_application_bigTauOnly = []

        etau_colors_application = []
        etau_colors_validity = []
        etau_colors_application_bigTau = []
        etau_colors_application_bigTauOnly = []

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
        df_eTau = dfWrapped_central.df.Filter('eTau').Filter("tau1_gen_kind==3 && tau2_gen_kind==5")
        #nInitial_eTau = df_eTau.Count().GetValue()
        total_weight_string = 'weight_total*weight_TauID_Central*weight_tau1_EleidSF_Central*weight_tau1_MuidSF_Central*weight_tau2_EleidSF_Central*weight_tau2_MuidSF_Central*weight_L1PreFiring_Central*weight_L1PreFiring_ECAL_Central*weight_L1PreFiring_Muon_Central'
        nInitial_eTau = df_eTau.Define('total_weight_string', total_weight_string).Sum('total_weight_string').GetValue()


        df_etau_vloose = df_eTau.Filter(os_iso_filtering['VLoose'])
        df_etau_loose = df_eTau.Filter(os_iso_filtering['Loose'])
        df_etau_medium = df_eTau.Filter(os_iso_filtering['Medium'])

        dataframes_channel = {
            #'eTau':df_eTau,
            #'muTau':df_muTau,
            #'eTau': df_eTau,
            'etau_medium': df_etau_medium,
            'etau_loose': df_etau_loose,
            'etau_vloose': df_etau_vloose,
        }
        dataframes_channel_bigTau = {
            #'eTau':df_eTau,
            #'muTau':df_muTau,
            #'tauTau': df_tauTau,
            'etau_medium': df_etau_medium,
            'etau_loose': df_etau_loose,
            'etau_vloose': df_etau_vloose,
        }
        for channel in ['etau_medium','etau_loose', 'etau_vloose']:
            regions_expr= get_regions_dict('eTau', [190,190], 150, False)
            regions_expr_bigTau= get_regions_dict('eTau', [190,190], 150, True)
            for reg_key,reg_exp in regions_expr.items():
                dataframes_channel[channel] = dataframes_channel[channel].Define(reg_key, reg_exp)
            for reg_key,reg_exp in regions_expr_bigTau.items():
                dataframes_channel_bigTau[channel] = dataframes_channel_bigTau[channel].Define(reg_key, reg_exp)

        pass_MET_application = "(HLT_MET && MET_region && !(singleTau_region) && ! (other_trg_region))"
        pass_singleTau_application = "(HLT_singleTau && singleTau_region && !(MET_region) && ! (other_trg_region))"
        pass_other_trg_application= "({} && other_trg_region  && !(MET_region) && ! (singleTau_region) )  "

        pass_singleTau_validity = f"""HLT_singleTau && ({GetTrgValidityRegion("HLT_singleTau")})"""
        pass_MET_validity = f"""HLT_MET && ({GetTrgValidityRegion("HLT_MET")})"""


        #### eTau efficiencies ####
        pass_etau_application = pass_other_trg_application.format("HLT_etau")
        pass_singleEle_application = pass_other_trg_application.format("HLT_singleEle")
        pass_etau_validity = f"""HLT_etau && ({GetTrgValidityRegion("HLT_etau")})"""
        pass_singleEle_validity = f"""HLT_singleEle && ({GetTrgValidityRegion("HLT_singleEle")})"""

        ##### validity triggers #####
        trg_eTau_list_validity_fullOR = [pass_etau_validity,pass_singleEle_validity,pass_singleTau_validity,pass_MET_validity]
        trg_eTau_list_validity_legacy_singleTau = [pass_etau_validity,pass_singleEle_validity,pass_singleTau_validity]
        trg_eTau_list_validity_legacy_MET = [pass_etau_validity,pass_singleEle_validity,pass_MET_validity]
        trg_eTau_list_validity_legacy = [pass_etau_validity,pass_singleEle_validity] # ratio denumerator

        ##### application triggers #####
        trg_eTau_list_application_fullOR = [pass_etau_application,pass_singleEle_application,pass_singleTau_application,pass_MET_application]
        trg_eTau_list_application_legacy_singleTau = [pass_etau_application,pass_singleEle_application,pass_singleTau_application]
        trg_eTau_list_application_legacy_MET = [pass_etau_application,pass_singleEle_application,pass_MET_application]
        trg_eTau_list_application_legacy = [pass_etau_application,pass_singleEle_application] # ratio denumerator

        ##### preparing dict for validity region plot #####
        channel=f'etau_{args.deepTauWP}'

        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_validity_legacy, 'eTau', nInitial_eTau,'eff_validity_legacy', eff_etau_validity,eff_etau_validity_errors, args.wantNEvents)
        etau_linestyles_validity.append('solid')
        etau_labels_validity.append('legacy')
        etau_colors_validity.append("gold")
        #print(mass_string)
        #print(eff_etau_validity['eff_validity_legacy'])

        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_validity_legacy_singleTau, 'eTau', nInitial_eTau,'eff_validity_legacySingleTau', eff_etau_validity,eff_etau_validity_errors, args.wantNEvents)
        etau_linestyles_validity.append('solid')
        etau_labels_validity.append("legacy||singleTau")
        etau_colors_validity.append("red")

        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_validity_legacy_MET, 'eTau', nInitial_eTau,'eff_validity_legacyMET', eff_etau_validity,eff_etau_validity_errors, args.wantNEvents)
        etau_linestyles_validity.append('solid')
        etau_labels_validity.append("legacy||MET")
        etau_colors_validity.append("cyan")

        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_validity_fullOR, 'eTau', nInitial_eTau,'eff_validity_fullOR', eff_etau_validity,eff_etau_validity_errors, args.wantNEvents)
        etau_linestyles_validity.append('solid')
        etau_labels_validity.append("all")
        etau_colors_validity.append("blue")

        ##### preparing dict for application region plot #####
        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_application_legacy, 'eTau', nInitial_eTau,'eff_application_legacy', eff_etau_application,eff_etau_application_errors, args.wantNEvents)
        etau_labels_application.append('legacy')
        etau_colors_application.append("gold")
        etau_linestyles_application.append('solid')

        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_application_legacy_singleTau, 'eTau', nInitial_eTau,'eff_application_legacySingleTau', eff_etau_application,eff_etau_application_errors, args.wantNEvents)
        etau_labels_application.append("legacy||singleTau")
        etau_colors_application.append("red")
        etau_linestyles_application.append('solid')

        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_application_legacy_MET, 'eTau', nInitial_eTau,'eff_application_legacyMET', eff_etau_application,eff_etau_application_errors, args.wantNEvents)
        etau_labels_application.append("legacy||MET")
        etau_colors_application.append("cyan")
        etau_linestyles_application.append('solid')

        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_application_fullOR, 'eTau', nInitial_eTau,'eff_application_fullOR', eff_etau_application,eff_etau_application_errors, args.wantNEvents)
        etau_labels_application.append("all")
        etau_colors_application.append("blue")
        etau_linestyles_application.append('solid')

        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_validity_fullOR, 'eTau', nInitial_eTau,'eff_validity_fullOR', eff_etau_application,eff_etau_application_errors, args.wantNEvents)
        etau_labels_application.append("fullOR")
        etau_colors_application.append("dodgerblue")
        etau_linestyles_application.append('dashed')

        ##### preparing dict for application region plot bigTau = False VS bigTau = True #####

        #### bigTau = False ####
        '''
        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_application_legacy, 'eTau', nInitial_eTau,'eff_application_legacy', eff_etau_application_bigTau,eff_etau_application_bigTau_errors, args.wantNEvents)
        etau_labels_application_bigTau.append('legacy')
        etau_colors_application_bigTau.append("gold")
        etau_linestyles_application_bigTau.append('solid')

        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_application_legacy_singleTau, 'eTau', nInitial_eTau,'eff_application_legacySingleTau', eff_etau_application_bigTau,eff_etau_application_bigTau_errors, args.wantNEvents)
        etau_labels_application_bigTau.append("legacy||singleTau")
        etau_colors_application_bigTau.append("red")
        etau_linestyles_application_bigTau.append('solid')

        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_application_legacy_MET, 'eTau', nInitial_eTau,'eff_application_legacyMET', eff_etau_application_bigTau,eff_etau_application_bigTau_errors, args.wantNEvents)
        etau_labels_application_bigTau.append("legacy||MET")
        etau_colors_application_bigTau.append("cyan")
        etau_linestyles_application_bigTau.append('solid')
        '''

        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_application_fullOR, 'eTau', nInitial_eTau,'eff_application_fullOR', eff_etau_application_bigTau,eff_etau_application_bigTau_errors, args.wantNEvents)
        etau_labels_application_bigTau.append("bigTau=False")
        etau_colors_application_bigTau.append("cyan")
        #etau_colors_application_bigTau.append("blue")
        etau_linestyles_application_bigTau.append('solid')

        #### bigTau = True ####
        '''
        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_eTau_list_application_legacy, 'eTau', nInitial_eTau,'eff_application_legacy_bigTau', eff_etau_application_bigTau,eff_etau_application_bigTau_errors, args.wantNEvents)
        etau_labels_application_bigTau.append('legacy')
        etau_colors_application_bigTau.append("olivedrab")
        etau_linestyles_application_bigTau.append('solid')

        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_eTau_list_application_legacy_singleTau, 'eTau', nInitial_eTau,'eff_application_legacySingleTau_bigTau', eff_etau_application_bigTau,eff_etau_application_bigTau_errors, args.wantNEvents)
        etau_labels_application_bigTau.append("legacy||singleTau")
        etau_colors_application_bigTau.append("maroon")
        etau_linestyles_application_bigTau.append('solid')

        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_eTau_list_application_legacy_MET, 'eTau', nInitial_eTau,'eff_application_legacyMET_bigTau', eff_etau_application_bigTau,eff_etau_application_bigTau_errors, args.wantNEvents)
        etau_labels_application_bigTau.append("legacy||MET")
        etau_colors_application_bigTau.append("paleturquoise")
        etau_linestyles_application_bigTau.append('solid')
        '''
        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_eTau_list_application_fullOR, 'eTau', nInitial_eTau,'eff_application_fullOR_bigTau', eff_etau_application_bigTau,eff_etau_application_bigTau_errors, args.wantNEvents)
        etau_labels_application_bigTau.append("bigTau=True")
        etau_colors_application_bigTau.append("midnightblue")
        etau_linestyles_application_bigTau.append('solid')

        AddEfficiencyToDict(dataframes_channel[channel], trg_eTau_list_validity_fullOR, 'eTau', nInitial_eTau,'eff_validity_fullOR', eff_etau_application_bigTau,eff_etau_application_bigTau_errors, args.wantNEvents)
        etau_labels_application_bigTau.append("fullOR")
        etau_colors_application_bigTau.append("dodgerblue")
        etau_linestyles_application_bigTau.append('dashed')

        ########## bigTau only ###########
        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_eTau_list_application_legacy, 'eTau', nInitial_eTau,'eff_application_legacy_bigTau', eff_etau_application_bigTauOnly,eff_etau_application_bigTauOnly_errors, args.wantNEvents)
        etau_labels_application_bigTauOnly.append('legacy')
        etau_colors_application_bigTauOnly.append("olivedrab")
        etau_linestyles_application_bigTauOnly.append('solid')

        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_eTau_list_application_legacy_singleTau, 'eTau', nInitial_eTau,'eff_application_legacySingleTau_bigTau', eff_etau_application_bigTauOnly,eff_etau_application_bigTauOnly_errors, args.wantNEvents)
        etau_labels_application_bigTauOnly.append("legacy||singleTau")
        etau_colors_application_bigTauOnly.append("maroon")
        etau_linestyles_application_bigTauOnly.append('solid')

        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_eTau_list_application_legacy_MET, 'eTau', nInitial_eTau,'eff_application_legacyMET_bigTau', eff_etau_application_bigTauOnly,eff_etau_application_bigTauOnly_errors, args.wantNEvents)
        etau_labels_application_bigTauOnly.append("legacy||MET")
        etau_colors_application_bigTauOnly.append("paleturquoise")
        etau_linestyles_application_bigTauOnly.append('solid')

        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_eTau_list_application_fullOR, 'eTau', nInitial_eTau,'eff_application_fullOR_bigTau', eff_etau_application_bigTauOnly,eff_etau_application_bigTauOnly_errors, args.wantNEvents)
        etau_labels_application_bigTauOnly.append("all")
        etau_colors_application_bigTauOnly.append("midnightblue")
        etau_linestyles_application_bigTauOnly.append('solid')

        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_eTau_list_validity_fullOR, 'eTau', nInitial_eTau,'eff_validity_fullOR', eff_etau_application_bigTauOnly,eff_etau_application_bigTauOnly_errors, args.wantNEvents)
        etau_labels_application_bigTauOnly.append("fullOR")
        etau_colors_application_bigTauOnly.append("dodgerblue")
        etau_linestyles_application_bigTauOnly.append('dashed')


    #print(eff_etau)
    x_values = masses
    ratio_ref = eff_etau_validity['eff_validity_legacy']
    ratio_ref_errors = eff_etau_validity_errors['eff_validity_legacy_error']

    # arguments: (eff_dict, labels, linestyles, channel, x_values, sample, deepTauWP, deepTauVersion, ratio_ref,suffix='', colors=[])



    makeplot(eff_etau_validity,eff_etau_validity_errors, etau_labels_validity, etau_linestyles_validity, 'eTau', x_values, args.sample, args.deepTauWP, args.deepTauVersion, ratio_ref,ratio_ref_errors,'_validity',etau_colors_validity, args.wantNEvents,args.wantLegend)

    makeplot(eff_etau_application,eff_etau_application_errors, etau_labels_application, etau_linestyles_application, 'eTau', x_values, args.sample, args.deepTauWP, args.deepTauVersion, ratio_ref,ratio_ref_errors,'_application',etau_colors_application, args.wantNEvents,args.wantLegend)

    makeplot(eff_etau_application_bigTau,eff_etau_application_bigTau_errors, etau_labels_application_bigTau, etau_linestyles_application_bigTau, 'eTau', x_values, args.sample,args.deepTauWP, args.deepTauVersion, ratio_ref,ratio_ref_errors,'_application_bigTau',etau_colors_application_bigTau, args.wantNEvents,args.wantLegend, True)

    makeplot(eff_etau_application_bigTauOnly,eff_etau_application_bigTauOnly_errors, etau_labels_application_bigTauOnly, etau_linestyles_application_bigTauOnly, 'eTau', x_values, args.sample,args.deepTauWP, args.deepTauVersion, ratio_ref,ratio_ref_errors,'_application_bigTauOnly',etau_colors_application_bigTauOnly, args.wantNEvents,args.wantLegend)



    plt.show()

    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))


