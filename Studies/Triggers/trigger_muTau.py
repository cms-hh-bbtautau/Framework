
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
    #inDir = "/afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/output/Run2_2018/{}/anaTuples/"
    masses = []

    eff_mutau_validity = {}
    eff_mutau_application = {}
    eff_mutau_application_bigTau={}
    eff_mutau_application_bigTauOnly={}

    eff_mutau_validity_errors  = {}
    eff_mutau_application_errors  = {}
    eff_mutau_application_bigTau_errors ={}
    eff_mutau_application_bigTauOnly_errors ={}

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
        mutau_labels_application = []
        mutau_labels_validity = []
        mutau_labels_application_bigTau = []
        mutau_labels_application_bigTauOnly = []

        mutau_linestyles_application = []
        mutau_linestyles_validity = []
        mutau_linestyles_application_bigTau = []
        mutau_linestyles_application_bigTauOnly = []

        mutau_colors_application = []
        mutau_colors_validity = []
        mutau_colors_application_bigTau = []
        mutau_colors_application_bigTauOnly = []

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
        df_muTau = dfWrapped_central.df.Filter('muTau').Filter("tau1_gen_kind==4 && tau2_gen_kind==5")

        #nInitial_muTau = df_muTau.Count().GetValue()
        total_weight_string = 'weight_total*weight_TauID_Central*weight_tau1_EleidSF_Central*weight_tau1_MuidSF_Central*weight_tau2_EleidSF_Central*weight_tau2_MuidSF_Central*weight_L1PreFiring_Central*weight_L1PreFiring_ECAL_Central*weight_L1PreFiring_Muon_Central'
        nInitial_muTau = df_muTau.Define('total_weight_string', total_weight_string).Sum('total_weight_string').GetValue()


        df_mutau_vloose = df_muTau.Filter(os_iso_filtering['VLoose'])
        df_mutau_loose = df_muTau.Filter(os_iso_filtering['Loose'])
        df_mutau_medium = df_muTau.Filter(os_iso_filtering['Medium'])
        dataframes_channel = {
            #'eTau':df_eTau,
            #'muTau':df_muTau,
            #'eTau': df_eTau,
            'mutau_medium': df_mutau_medium,
            'mutau_loose': df_mutau_loose,
            'mutau_vloose': df_mutau_vloose,
        }
        dataframes_channel_bigTau = {
            #'eTau':df_eTau,
            #'muTau':df_muTau,
            #'tauTau': df_tauTau,
            'mutau_medium': df_mutau_medium,
            'mutau_loose': df_mutau_loose,
            'mutau_vloose': df_mutau_vloose,
        }
        for channel in ['mutau_medium','mutau_loose', 'mutau_vloose']:
            regions_expr= get_regions_dict('muTau', [190,190], 150, False)
            regions_expr_bigTau= get_regions_dict('muTau', [190,190], 150, True)
            for reg_key,reg_exp in regions_expr.items():
                dataframes_channel[channel] = dataframes_channel[channel].Define(reg_key, reg_exp)
            for reg_key,reg_exp in regions_expr_bigTau.items():
                dataframes_channel_bigTau[channel] = dataframes_channel_bigTau[channel].Define(reg_key, reg_exp)

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

        ##### validity triggers #####
        trg_muTau_list_validity_fullOR = [pass_mutau_validity,pass_singleMu_validity,pass_singleTau_validity,pass_MET_validity]
        trg_muTau_list_validity_fullOR_mu50 = [pass_mutau_validity,pass_singleMu_validity,pass_singleTau_validity,pass_MET_validity,pass_singleMu50_validity]
        trg_muTau_list_validity_legacy_singleTau = [pass_mutau_validity,pass_singleMu_validity,pass_singleTau_validity]
        trg_muTau_list_validity_legacy_MET = [pass_mutau_validity,pass_singleMu_validity,pass_MET_validity]
        trg_muTau_list_validity_legacy_MET_singleTau = [pass_mutau_validity,pass_singleMu_validity,pass_MET_validity,pass_singleTau_validity]
        trg_muTau_list_validity_legacy_singleMu50 = [pass_mutau_validity,pass_singleMu_validity,pass_singleMu50_validity]
        trg_muTau_list_validity_legacy = [pass_mutau_validity,pass_singleMu_validity] # ratio denumerator

        ##### application triggers #####
        trg_muTau_list_application_fullOR = [pass_mutau_application,pass_singleMu_application,pass_singleTau_application,pass_MET_application]
        #trg_muTau_list_application_fullOR = [pass_mutau_application,pass_singleMu_application,pass_singleMu50_application,pass_singleTau_application,pass_MET_application]
        trg_muTau_list_application_legacy_singleTau = [pass_mutau_application,pass_singleMu_application,pass_singleTau_application]
        trg_muTau_list_application_legacy_MET_singleTau = [pass_mutau_application,pass_singleMu_application,pass_MET_application, pass_singleTau_application]
        trg_muTau_list_application_legacy_MET = [pass_mutau_application,pass_singleMu_application,pass_MET_application]
        trg_muTau_list_application_legacy_singleMu50 = [pass_mutau_application,pass_singleMu_application,pass_singleMu50_application]
        trg_muTau_list_application_legacy = [pass_mutau_application,pass_singleMu_application] # ratio denumerator

        ######################################################
        ##### preparing dict for validity region plot #####
        ######################################################
        channel=f'mutau_{args.deepTauWP}'

        AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_validity_legacy, 'muTau', nInitial_muTau,'eff_validity_legacy', eff_mutau_validity, eff_mutau_validity_errors, args.wantNEvents)
        mutau_linestyles_validity.append('solid')
        mutau_labels_validity.append('legacy')
        mutau_colors_validity.append("gold")

        AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_validity_legacy_singleTau, 'muTau', nInitial_muTau,'eff_validity_legacySingleTau', eff_mutau_validity, eff_mutau_validity_errors, args.wantNEvents)
        mutau_linestyles_validity.append('solid')
        mutau_labels_validity.append("legacy||singleTau")
        mutau_colors_validity.append("red")

        AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_validity_legacy_MET, 'muTau', nInitial_muTau,'eff_validity_legacyMET', eff_mutau_validity, eff_mutau_validity_errors, args.wantNEvents)
        mutau_linestyles_validity.append('solid')
        mutau_labels_validity.append("legacy||MET")
        mutau_colors_validity.append("cyan")
        '''

        AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_validity_fullOR_mu50, 'muTau', nInitial_muTau,'eff_validity_fullOR_mu50', eff_mutau_validity, eff_mutau_validity_errors, args.wantNEvents)
        mutau_linestyles_validity.append('solid')
        mutau_labels_validity.append("fullOR-mu50")
        mutau_colors_validity.append("cyan")

        '''

        AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_validity_fullOR, 'muTau', nInitial_muTau,'eff_validity_fullOR', eff_mutau_validity, eff_mutau_validity_errors, args.wantNEvents)
        mutau_linestyles_validity.append('solid')
        mutau_labels_validity.append("fullOR")
        mutau_colors_validity.append("blue")

        ######################################################
        ##### preparing dict for application region plot #####
        ######################################################
        AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_application_legacy, 'muTau', nInitial_muTau,'eff_application_legacy', eff_mutau_application, eff_mutau_application_errors, args.wantNEvents)
        mutau_labels_application.append('legacy')
        mutau_colors_application.append("gold")
        mutau_linestyles_application.append('solid')

        AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_application_legacy_singleTau, 'muTau', nInitial_muTau,'eff_application_legacySingleTau', eff_mutau_application, eff_mutau_application_errors, args.wantNEvents)
        mutau_labels_application.append("legacy||singleTau")
        mutau_colors_application.append("red")
        mutau_linestyles_application.append('solid')


        AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_application_legacy_MET, 'muTau', nInitial_muTau,'eff_application_legacyMET', eff_mutau_application, eff_mutau_application_errors, args.wantNEvents)
        mutau_labels_application.append("legacy||MET")
        mutau_colors_application.append("cyan")
        mutau_linestyles_application.append('solid')


        AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_application_fullOR, 'muTau', nInitial_muTau,'eff_application_fullOR', eff_mutau_application, eff_mutau_application_errors, args.wantNEvents)
        mutau_labels_application.append("all")
        mutau_colors_application.append("blue")
        mutau_linestyles_application.append('solid')

        AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_validity_fullOR, 'muTau', nInitial_muTau,'eff_validity_fullOR', eff_mutau_application, eff_mutau_application_errors, args.wantNEvents)
        mutau_labels_application.append("fullOR")
        mutau_colors_application.append("dodgerblue")
        mutau_linestyles_application.append('dashed')

        ############################################################################################################
        ##### preparing dict for application region plot bigTau = False VS bigTau = True #####
        ############################################################################################################

        #### bigTau = False ####
        #AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_application_legacy, 'muTau', nInitial_muTau,'eff_application_legacy', eff_mutau_application_bigTau, eff_mutau_application_bigTau_errors, args.wantNEvents)
        #mutau_labels_application_bigTau.append('legacy')
        #mutau_colors_application_bigTau.append("gold")
        #mutau_linestyles_application_bigTau.append('solid')

        #AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_application_legacy_singleTau, 'muTau', nInitial_muTau,'eff_application_legacySingleTau', eff_mutau_application_bigTau, eff_mutau_application_bigTau_errors, args.wantNEvents)
        #mutau_labels_application_bigTau.append("legacy||singleTau")
        #mutau_colors_application_bigTau.append("red")
        #mutau_linestyles_application_bigTau.append('solid')

        #AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_application_legacy_MET, 'muTau', nInitial_muTau,'eff_application_legacyMET', eff_mutau_application_bigTau, eff_mutau_application_bigTau_errors, args.wantNEvents)
        #mutau_labels_application_bigTau.append("legacy||MET")
        #mutau_colors_application_bigTau.append("cyan")
        #mutau_linestyles_application_bigTau.append('solid')

        AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_application_fullOR, 'muTau', nInitial_muTau,'eff_application_fullOR', eff_mutau_application_bigTau, eff_mutau_application_bigTau_errors,  args.wantNEvents, True)
        mutau_labels_application_bigTau.append("bigTau=False")
        mutau_colors_application_bigTau.append("cyan")
        #mutau_colors_application_bigTau.append("blue")
        mutau_linestyles_application_bigTau.append('solid')

        #### bigTau = True ####
        #AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_muTau_list_application_legacy, 'muTau', nInitial_muTau,'eff_application_legacy_bigTau', eff_mutau_application_bigTau, eff_mutau_application_bigTau_errors, args.wantNEvents)
        #mutau_labels_application_bigTau.append('legacy - bigTau')
        #mutau_colors_application_bigTau.append("olivedrab")
        #mutau_linestyles_application_bigTau.append('solid')

        #AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_muTau_list_application_legacy_singleTau, 'muTau', nInitial_muTau,'eff_application_legacySinglmuTau_bigTau', eff_mutau_application_bigTau, eff_mutau_application_bigTau_errors, args.wantNEvents)
        #mutau_labels_application_bigTau.append("legacy||singleTau - bigTau")
        #mutau_colors_application_bigTau.append("maroon")
        #mutau_linestyles_application_bigTau.append('solid')

        #AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_muTau_list_application_legacy_MET, 'muTau', nInitial_muTau,'eff_application_legacyMET_bigTau', eff_mutau_application_bigTau, eff_mutau_application_bigTau_errors, args.wantNEvents)
        #mutau_labels_application_bigTau.append("legacy||MET - bigTau")
        #mutau_colors_application_bigTau.append("paleturquoise")
        #mutau_linestyles_application_bigTau.append('solid')


        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_muTau_list_application_fullOR, 'muTau', nInitial_muTau,'eff_application_fullOR_bigTau', eff_mutau_application_bigTau, eff_mutau_application_bigTau_errors, args.wantNEvents, True)
        mutau_labels_application_bigTau.append("bigTau=True")
        mutau_colors_application_bigTau.append("midnightblue")
        mutau_linestyles_application_bigTau.append('solid')

        AddEfficiencyToDict(dataframes_channel[channel], trg_muTau_list_validity_fullOR, 'muTau', nInitial_muTau,'eff_validity_fullOR', eff_mutau_application_bigTau, eff_mutau_application_bigTau_errors, args.wantNEvents)
        mutau_labels_application_bigTau.append("fullOR")
        mutau_colors_application_bigTau.append("dodgerblue")
        mutau_linestyles_application_bigTau.append('dashed')

        ###################################
        ########### bigTau only ###########
        ###################################

        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_muTau_list_application_legacy, 'muTau', nInitial_muTau,'eff_application_legacy_bigTau', eff_mutau_application_bigTauOnly, eff_mutau_application_bigTauOnly_errors, args.wantNEvents)
        mutau_labels_application_bigTauOnly.append('legacy')
        mutau_colors_application_bigTauOnly.append("olivedrab")
        mutau_linestyles_application_bigTauOnly.append('solid')

        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_muTau_list_application_legacy_singleTau, 'muTau', nInitial_muTau,'eff_application_legacySinglmuTau_bigTau', eff_mutau_application_bigTauOnly, eff_mutau_application_bigTauOnly_errors, args.wantNEvents)
        mutau_labels_application_bigTauOnly.append("legacy||singleTau")
        mutau_colors_application_bigTauOnly.append("maroon")
        mutau_linestyles_application_bigTauOnly.append('solid')

        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_muTau_list_application_legacy_MET, 'muTau', nInitial_muTau,'eff_application_legacyMET_bigTau', eff_mutau_application_bigTauOnly, eff_mutau_application_bigTauOnly_errors, args.wantNEvents)
        mutau_labels_application_bigTauOnly.append("legacy||MET")
        mutau_colors_application_bigTauOnly.append("paleturquoise")
        mutau_linestyles_application_bigTauOnly.append('solid')

        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_muTau_list_application_fullOR, 'muTau', nInitial_muTau,'eff_application_fullOR_bigTau', eff_mutau_application_bigTauOnly, eff_mutau_application_bigTauOnly_errors, args.wantNEvents)
        mutau_labels_application_bigTauOnly.append("all")
        mutau_colors_application_bigTauOnly.append("midnightblue")
        mutau_linestyles_application_bigTauOnly.append('solid')

        AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_muTau_list_validity_fullOR, 'muTau', nInitial_muTau,'eff_validity_fullOR', eff_mutau_application_bigTauOnly, eff_mutau_application_bigTauOnly_errors, args.wantNEvents)
        mutau_labels_application_bigTauOnly.append("fullOR")
        mutau_colors_application_bigTauOnly.append("dodgerblue")
        mutau_linestyles_application_bigTauOnly.append('dashed')

    #print(eff_mutau)
    x_values = masses
    ratio_ref = eff_mutau_validity['eff_validity_legacy']
    ratio_ref_errors = eff_mutau_validity_errors['eff_validity_legacy_error']

    makeplot(eff_mutau_validity,eff_mutau_validity_errors, mutau_labels_validity, mutau_linestyles_validity, 'muTau', x_values, args.sample, args.deepTauWP, args.deepTauVersion, ratio_ref,ratio_ref_errors,'_validity',mutau_colors_validity, args.wantNEvents, args.wantLegend)

    makeplot(eff_mutau_application,eff_mutau_application_errors, mutau_labels_application, mutau_linestyles_application, 'muTau', x_values, args.sample, args.deepTauWP, args.deepTauVersion, ratio_ref,ratio_ref_errors,'_application',mutau_colors_application, args.wantNEvents, args.wantLegend)

    makeplot(eff_mutau_application_bigTau,eff_mutau_application_bigTau_errors, mutau_labels_application_bigTau, mutau_linestyles_application_bigTau, 'muTau', x_values, args.sample,args.deepTauWP, args.deepTauVersion, ratio_ref,ratio_ref_errors,'_application_bigTau',mutau_colors_application_bigTau, args.wantNEvents, args.wantLegend, True)

    makeplot(eff_mutau_application_bigTauOnly,eff_mutau_application_bigTauOnly_errors, mutau_labels_application_bigTauOnly, mutau_linestyles_application_bigTauOnly, 'muTau', x_values, args.sample,args.deepTauWP, args.deepTauVersion, ratio_ref,ratio_ref_errors,'_application_bigTauOnly',mutau_colors_application_bigTauOnly, args.wantNEvents, args.wantLegend)




    plt.show()

    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))


