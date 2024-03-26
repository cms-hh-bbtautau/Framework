
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

def AddEfficiencyToDict(df,trg_list, channel, n_initial_channel,eff_key, eff_dict,verbose=False):
    filter_expr = f' ('
    filter_expr+= ' || '.join(trg for trg in trg_list)
    filter_expr+= ')'
    n_channel = df.Filter(filter_expr).Count().GetValue()
    eff_channel = n_channel / n_initial_channel
    if verbose:
        print(f"with {filter_expr} : n_initial{channel} = {n_initial_channel}, n_{channel} = {n_channel}, eff_{channel} = {round(eff_channel,2)}")
    if eff_key not in eff_dict.keys():
        eff_dict[eff_key] = []
    eff_dict[eff_key].append(round(eff_channel,2))

def makeplot(eff_dict, labels, linestyles, channel, x_values, sample,wantBigTau,deepTauWP, deepTauVersion,noreg,ratio_ref,suffix=''):
    colors = ['blue', 'green', 'red', 'orange', 'purple', 'pink', 'yellow', 'cyan','black','brown','lime','navy','crimson']
    markers = ['o', '^', 's', 'D', 'x', 'v', 'p', '*','o']
    plt.figure(figsize=(20,15))
    #print(labels)
    key_baseline = ''
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    for i, (key, values) in enumerate(eff_dict.items()):
        print(i, key, key_baseline)
        for channel in ['tautau_medium','tautau_loose','tautau_vloose']:
            if key == f'eff_{channel}_noreg_diTau':
                key_baseline = key
                ratio_ref = values
        if key.split('_')[-1]=='diTau':
            continue

        print(i, key, key_baseline, values)
        label = None
        if i < len(labels):
            label = labels[i]
        #print(label)

        #print(len(values), len(ratio_ref))
        ratio_values = [(val/ref) for val,ref in zip(values,ratio_ref)]
        #print(x_values)
        #print(ratio_values)
        for mass,ratio,value in zip(x_values,ratio_values,values):
            print(f"mass = {mass} and value = {ratio} and value = {value}")
        #print(colors[i % len(colors)],markers[i % len(markers)],linestyles[i % len(markers)], label)
        ax1.plot(x_values, values, color=colors[i % len(colors)],marker=markers[i % len(markers)],linestyle=linestyles[i % len(markers)])#, label=label)
        ax2.plot(x_values, ratio_values, color=colors[i % len(colors)],marker=markers[i % len(markers)],linestyle=linestyles[i % len(markers)],label=label)

    #### Legend + titles + axis adjustment ####
    #plt.title(f'efficiencies for channel {channel} and {sample}')
    plt.xlabel(r'$m_X$')
    ax1.set_ylabel('efficiency')
    ax2.set_ylabel('ratio')
    ax1.set_ylim(0., 1.05)
    #ax2.legend(loc='upper right',  bbox_to_anchor=(0.5, 0.8))
    ax2.legend(loc='upper left',  bbox_to_anchor=(0.8, 0.3))
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    #ax2.legend(loc = 'best')#lines, labels, loc = 'lower center', ncol=5, labelspacing=0.)
    #plt.figlegend(loc='lower left')#lines, labels, loc = 'lower center', ncol=5, labelspacing=0.)
    figName = f'Studies/Triggers/eff_{channel}_{sample}_{deepTauWP}_{deepTauVersion}{suffix}'
    if wantBigTau:
        figName+='_bigTau'
    if noreg:
        figName+='_noReg'
    plt.tight_layout()
    plt.savefig(f'{figName}.png')

'''
def create_dict(dict_to_fill, keys=['key_names','values', 'linestyles','labels']):
    for key in keys:
        if key not in dict_to_fill:
            dict_to_fill[key] = []
'''

if __name__ == "__main__":
    import argparse
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=False, type=str, default='GluGluToBulkGraviton')
    parser.add_argument('--deepTauVersion', required=False, type=str, default='v2p1')
    parser.add_argument('--wantBigTau', required=False, type=bool, default=False)
    parser.add_argument('--wantTightFirstTau', required=False, type=bool, default=False)
    #parser.add_argument('--deepTauWP', required=False, type=str, default='VLoose')


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

    eff_etau = {}
    eff_mutau = {}
    eff_tautau = {}

    eff_etau_noreg = {}
    eff_mutau_noreg = {}
    eff_tautau_noreg = {}

    eff_etau_bigTau={}
    eff_mutau_bigTau={}
    eff_tautau_bigTau={}
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
        #create_dict(eff_etau, keys=['key_names','values', 'linestyles','labels'])
        #create_dict(eff_mutau, keys=['key_names','values', 'linestyles','labels'])
        #create_dict(eff_tautau, keys=['key_names','values', 'linestyles','labels'])
        #create_dict(eff_etau_noreg, keys=['key_names','values', 'linestyles','labels'])
        #create_dict(eff_mutau_noreg, keys=['key_names','values', 'linestyles','labels'])
        #create_dict(eff_tautau_noreg, keys=['key_names','values', 'linestyles','labels'])
        #create_dict(eff_etau_bigTau, keys=['key_names','values', 'linestyles','labels'])
        #create_dict(eff_mutau_bigTau, keys=['key_names','values', 'linestyles','labels'])
        #create_dict(eff_tautau_bigTau, keys=['key_names','values', 'linestyles','labels'])
        etau_labels = []
        mutau_labels = []
        tautau_labels = []
        etau_labels_noreg = []
        mutau_labels_noreg = []
        tautau_labels_noreg = []
        etau_linestyles = []
        mutau_linestyles = []
        tautau_linestyles = []
        tautau_labels_bigTau = []
        tautau_linestyles_bigTau = []
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
        '''

        nInitial_eTau_0 = dfWrapped_central.df.Filter('eTau').Count().GetValue()
        nInitial_muTau_0 = dfWrapped_central.df.Filter('muTau').Count().GetValue()
        nInitial_tauTau_0 = dfWrapped_central.df.Filter('tauTau').Count().GetValue()


        df_eTau = dfWrapped_central.df.Filter('eTau').Filter("tau1_gen_kind==3 && tau2_gen_kind==5")
        nInitial_eTau = df_eTau.Count().GetValue()
        df_eTau = df_eTau.Filter(os_iso_filtering[args.deepTauWP])

        df_muTau = dfWrapped_central.df.Filter('muTau').Filter("tau1_gen_kind==4 && tau2_gen_kind==5")
        nInitial_muTau = df_muTau.Count().GetValue()
        df_muTau = df_muTau.Filter(os_iso_filtering[args.deepTauWP])
        '''
        df_tauTau = dfWrapped_central.df.Filter('tauTau').Filter("tau1_gen_kind==5 && tau2_gen_kind==5")
        nInitial_tauTau = df_tauTau.Count().GetValue()


        df_tautau_vloose = df_tauTau.Filter(os_iso_filtering['VLoose'])
        df_tautau_loose = df_tauTau.Filter(os_iso_filtering['Loose'])
        df_tautau_medium = df_tauTau.Filter(os_iso_filtering['Medium'])
        if args.wantTightFirstTau:
            df_tautau_vloose.Filter("tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}")
            df_tautau_medium.Filter("tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}")
            df_tautau_loose.Filter("tau1_idDeepTau{deepTauYear}{args.deepTauVersion}VSjet >= {Utilities.WorkingPointsTauVSjet.Tight.value}")

        dataframes_channel = {
            #'eTau':df_eTau,
            #'muTau':df_muTau,
            #'tauTau': df_tauTau,
            'tautau_medium': df_tautau_medium,
            'tautau_loose': df_tautau_loose,
            'tautau_vloose': df_tautau_vloose,
        }
        dataframes_channel_bigTau = {
            #'eTau':df_eTau,
            #'muTau':df_muTau,
            #'tauTau': df_tauTau,
            'tautau_medium': df_tautau_medium,
            'tautau_loose': df_tautau_loose,
            'tautau_vloose': df_tautau_vloose,
        }
        for channel in ['tautau_medium','tautau_loose', 'tautau_vloose']:
            regions_expr= get_regions_dict('tauTau', [190,190], 180, False)
            regions_expr_bigTau= get_regions_dict('tauTau', [190,190], 180, True)
            for reg_key,reg_exp in regions_expr.items():
                dataframes_channel[channel] = dataframes_channel[channel].Define(reg_key, reg_exp)
            for reg_key,reg_exp in regions_expr_bigTau.items():
                dataframes_channel_bigTau[channel] = dataframes_channel_bigTau[channel].Define(reg_key, reg_exp)
        pass_met = "(HLT_MET && MET_region && !(singleTau_region) && ! (other_trg_region))"
        pass_singleTau = "(HLT_singleTau && singleTau_region && !(MET_region) && ! (other_trg_region))"
        pass_other_trg = "({} && other_trg_region  && !(MET_region) && ! (singleTau_region) )  "

        pass_singleTau_validity = f"""HLT_singleTau && ({GetTrgValidityRegion("HLT_singleTau")})"""
        pass_MET_validity = f"""HLT_MET && ({GetTrgValidityRegion("HLT_MET")})"""


        #### tauTau efficiencies ####
        pass_diTau = pass_other_trg.format("HLT_ditau")
        pass_ditau_validity = f"""HLT_ditau && ({GetTrgValidityRegion("HLT_ditau")})"""
        trg_tauTau_list_noreg_diTauOnly = [pass_ditau_validity] # denumerator
        trg_tauTau_list_fullOR = [pass_ditau_validity,pass_singleTau_validity,pass_MET_validity]
        trg_tauTau_list_reg_all = [pass_diTau, pass_singleTau, pass_met]



        for channel in ['tautau_medium','tautau_loose','tautau_vloose']:
            AddEfficiencyToDict(dataframes_channel[channel], trg_tauTau_list_noreg_diTauOnly, 'tauTau', nInitial_tauTau,f'eff_{channel}_noreg_diTau', eff_tautau_bigTau)
            tautau_labels_bigTau.append(channel)
            tautau_linestyles_bigTau.append('solid')
            AddEfficiencyToDict(dataframes_channel[channel], trg_tauTau_list_reg_all, 'tauTau', nInitial_tauTau,f'eff_{channel}_regAll', eff_tautau_bigTau)
            tautau_labels_bigTau.append(channel)
            tautau_linestyles_bigTau.append('solid')
            AddEfficiencyToDict(dataframes_channel_bigTau[channel], trg_tauTau_list_reg_all, 'tauTau', nInitial_tauTau,f'eff_{channel}_regAll_bigTau', eff_tautau_bigTau)
            tautau_labels_bigTau.append(f"{channel} - BigTau")
            tautau_linestyles_bigTau.append('solid')
            #AddEfficiencyToDict(dataframes_channel[channel], trg_tauTau_list_fullOR, 'tauTau', nInitial_tauTau,f'eff_{channel}_fullOR', eff_tautau_bigTau)
            #tautau_labels_bigTau.append(f"{channel} - fullOR")
            #tautau_linestyles_bigTau.append('dashed')
        for a,b in zip (tautau_labels_bigTau,tautau_linestyles_bigTau):
            print(a, b)


    #print(eff_etau)
    x_values = masses
    makeplot(eff_tautau_bigTau, tautau_labels_bigTau,tautau_linestyles_bigTau, 'tauTau', x_values, args.sample, args.wantBigTau, 'allWPs', args.deepTauVersion,False, eff_tautau_bigTau['eff_tautau_medium_noreg_diTau'],'_bigTauComp')
    plt.show()

    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))


