import ROOT
import sys
import os
import math
import shutil
import time
ROOT.EnableThreadSafety()

from RunKit.run_tools import ps_call
if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])

import Common.Utilities as Utilities
from Analysis.HistHelper import *
from Analysis.hh_bbtautau import *

def createCacheQuantities(dfWrapped_cache, cache_map_name):
    df_cache = dfWrapped_cache.df
    map_creator_cache = ROOT.analysis.CacheCreator(*dfWrapped_cache.colTypes)()
    df_cache = map_creator_cache.processCache(ROOT.RDF.AsRNode(df_cache), Utilities.ListToVector(dfWrapped_cache.colNames), cache_map_name)
    return df_cache

def createVoidHist(outFileName, hist_cfg_dict):
    x_bins = hist_cfg_dict[args.var]['x_bins']
    if type(hist_cfg_dict[args.var]['x_bins'])==list:
        x_bins_vec = Utilities.ListToVector(x_bins, "double")
        hvoid = ROOT.TH1F("", "", x_bins_vec.size()-1, x_bins_vec.data())
    else:
        n_bins, bin_range = x_bins.split('|')
        start,stop = bin_range.split(':')
        hvoid = ROOT.TH1F("", "",int(n_bins), float(start), float(stop))
    outFile = ROOT.TFile(outFileName, "RECREATE")
    hvoid.Write()
    outFile.Close()

def AddCacheColumnsInDf(dfWrapped_central, dfWrapped_cache,cache_map_name='cache_map_placeholder'):
    col_names_cache =  dfWrapped_cache.colNames
    col_tpyes_cache =  dfWrapped_cache.colTypes
    #print(col_names_cache)
    #if "kinFit_result" in col_names_cache:
    #    col_names_cache.remove("kinFit_result")
    dfWrapped_cache.df = createCacheQuantities(dfWrapped_cache, cache_map_name)
    if dfWrapped_cache.df.Filter(f"{cache_map_name} > 0").Count().GetValue() <= 0 : raise RuntimeError("no events passed map placeolder")
    dfWrapped_central.AddCacheColumns(col_names_cache,col_tpyes_cache)

    #ROOT.gInterpreter.ProcessLine("""delete analysis::GetEntriesMap();""")

def createCentralQuantities(df_central, central_col_types, central_columns):
    map_creator = ROOT.analysis.MapCreator(*central_col_types)()
    df_central = map_creator.processCentral(ROOT.RDF.AsRNode(df_central), Utilities.ListToVector(central_columns), 1)
    #df_central = map_creator.getEventIdxFromShifted(ROOT.RDF.AsRNode(df_central))
    return df_central

def SaveHists(histograms, out_file):
    #print(histograms.items())
    for key_tuple,hist_list in histograms.items():
        (key_1, key_2) = key_tuple
        ch, reg, cat = key_1
        sample_type,uncName,scale = key_2
        #if cat == 'btag_shape': continue
        dir_name = '/'.join(key_1)
        dir_ptr = mkdir(out_file,dir_name)
        merged_hist = hist_list[0].GetValue()
        for hist in hist_list[1:] :
            merged_hist.Add(hist.GetValue())
            #print(key_tuple, merged_hist.Integral(0, merged_hist.GetNbinsX()+1))
        isCentral = 'Central' in key_2
        hist_name =  sample_type
        if not isCentral:
            hist_name+=f"_{uncName}{scale}"
        #print(key_tuple, merged_hist.Integral(0, merged_hist.GetNbinsX()+1))
        dir_ptr.WriteTObject(merged_hist, hist_name, "Overwrite")


def GetHistogramDictFromDataframes(var, all_dataframes, key_2 , key_filter_dict, unc_cfg_dict,hist_cfg_dict, furtherCut=''):
    dataframes = all_dataframes[key_2]
    sample_type,uncName,scale = key_2
    isCentral = 'Central' in key_2
    histograms = {}

    for key_1,key_cut in key_filter_dict.items():
        ch, reg, cat = key_1
        if (key_1, key_2) in histograms.keys(): continue
        if cat == 'boosted' and var in bjet_vars: continue
        if cat == 'boosted' and uncName in unc_to_not_consider_boosted: continue
        if cat != 'boosted' and var in var_to_add_boosted: continue

        total_weight_expression = GetWeight(ch,cat) if sample_type!='data' else "1"
        weight_name = "final_weight"
        if not isCentral:
            if type(unc_cfg_dict)==dict:
                if uncName in unc_cfg_dict.keys() and 'expression' in unc_cfg_dict[uncName].keys():
                    weight_name = unc_cfg_dict[uncName]['expression'].format(scale=scale)
        if (key_1, key_2) not in histograms.keys():
            histograms[(key_1, key_2)] = []
        for dataframe in dataframes:
            if "weight_tau1_TrgSF_singleTau_Central" not in dataframe.GetColumnNames():
                dataframe=dataframe.Define("weight_tau1_TrgSF_singleTau_Central","1.f")
            if "weight_tau2_TrgSF_singleTau_Central" not in dataframe.GetColumnNames():
                dataframe=dataframe.Define("weight_tau2_TrgSF_singleTau_Central","1.f")
            if furtherCut != '' : key_cut += f' && {furtherCut}'
            dataframe_new = dataframe.Filter(key_cut)
            dataframe_new = dataframe_new.Define(f"final_weight_0_{ch}_{cat}_{reg}", f"{total_weight_expression}")
            final_string_weight = ApplyBTagWeight(cat,applyBtag=False, finalWeight_name = f"final_weight_0_{ch}_{cat}_{reg}") if sample_type!='data' else "1"
            dataframe_new = dataframe_new.Filter(f"{cat}")
            if cat == 'btag_shape':
                final_string_weight = f"final_weight_0_{ch}_{cat}_{reg}"
            histograms[(key_1, key_2)].append(dataframe_new.Define("final_weight", final_string_weight).Define("weight_for_hists", f"{weight_name}").Histo1D(GetModel(hist_cfg_dict, var), var, "weight_for_hists"))

    return histograms


def GetShapeDataFrameDict(all_dataframes, key, key_central, inFile, inFileCache, compute_variations, deepTauVersion, colNames, colTypes, hasCache=True ):
    sample_type,uncName,scale=key
    if compute_variations and key!=key_central and sample_type!='data':
        if key not in all_dataframes.keys():
            all_dataframes[key] = []

        fileToOpen = ROOT.TFile(inFile, 'READ')
        file_keys= []
        for keyFile in fileToOpen.GetListOfKeys():
            if keyFile.GetName() == 'Events' : continue
            obj = keyFile.ReadObj()
            if not obj.IsA().InheritsFrom(ROOT.TTree.Class()):
                continue
            file_keys.append(keyFile.GetName())
        fileToOpen.Close()
        #print(file_keys)
        treeName = f"Events_{uncName}{scale}"
        #treeName = f"Events_nanoHTT_{uncName}{scale}"
        #print(treeName)
        treeName_noDiff = f"{treeName}_noDiff"
        if treeName_noDiff in file_keys:
            #print(treeName_noDiff)
            dfWrapped_noDiff = DataFrameBuilder(ROOT.RDataFrame(treeName_noDiff, inFile))
            dfWrapped_noDiff.CreateFromDelta(colNames, colTypes)
            if hasCache:
                dfWrapped_cache_noDiff = DataFrameBuilder(ROOT.RDataFrame(treeName_noDiff,inFileCache), args.deepTauVersion)
                AddCacheColumnsInDf(dfWrapped_noDiff, dfWrapped_cache_noDiff,f"cache_map_{uncName}{scale}_noDiff")
            all_dataframes[key].append(PrepareDfWrapped(dfWrapped_noDiff).df)


        treeName_Valid = f"{treeName}_Valid"
        if treeName_Valid in file_keys:
            #print(treeName_Valid)
            dfWrapped_Valid = DataFrameBuilder(ROOT.RDataFrame(treeName_Valid, inFile))
            dfWrapped_Valid.CreateFromDelta(colNames, colTypes)
            if hasCache:
                dfWrapped_cache_Valid = DataFrameBuilder(ROOT.RDataFrame(treeName_Valid,inFileCache), args.deepTauVersion)
                AddCacheColumnsInDf(dfWrapped_Valid, dfWrapped_cache_Valid,f"cache_map_{uncName}{scale}_Valid")
            all_dataframes[key].append(PrepareDfWrapped(dfWrapped_Valid).df)


        treeName_nonValid = f"{treeName}_nonValid"
        if treeName_nonValid in file_keys:
            #print(treeName_nonValid)
            dfWrapped_nonValid = DataFrameBuilder(ROOT.RDataFrame(treeName_nonValid, inFile))
            if hasCache:
                dfWrapped_cache_nonValid = DataFrameBuilder(ROOT.RDataFrame(treeName_nonValid,inFileCache), args.deepTauVersion)
                AddCacheColumnsInDf(dfWrapped_nonValid, dfWrapped_cache_nonValid, f"cache_map_{uncName}{scale}_nonValid")
            all_dataframes[key].append(PrepareDfWrapped(dfWrapped_nonValid).df)

        if not all_dataframes[key]:
            all_dataframes.pop(key)



if __name__ == "__main__":
    import argparse
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('--inFile', required=True, type=str)
    parser.add_argument('--cacheFile', required=False, type=str, default = '')
    parser.add_argument('--outFileName', required=True, type=str)
    parser.add_argument('--dataset', required=True, type=str)
    parser.add_argument('--deepTauVersion', required=False, type=str, default='v2p1')
    parser.add_argument('--compute_unc_variations', type=bool, default=False)
    parser.add_argument('--compute_rel_weights', type=bool, default=False)
    parser.add_argument('--histConfig', required=True, type=str)
    parser.add_argument('--globalConfig', required=True, type=str)
    parser.add_argument('--uncConfig', required=True, type=str)
    parser.add_argument('--var', required=True, type=str)
    parser.add_argument('--sampleConfig', required=True, type=str)
    parser.add_argument('--furtherCut', required=False, type=str, default = "")
    args = parser.parse_args()


    startTime = time.time()
    headers_dir = os.path.dirname(os.path.abspath(__file__))
    ROOT.gROOT.ProcessLine(f".include {os.environ['ANALYSIS_PATH']}")
    ROOT.gInterpreter.Declare(f'#include "include/KinFitNamespace.h"')
    ROOT.gInterpreter.Declare(f'#include "include/HistHelper.h"')
    ROOT.gInterpreter.Declare(f'#include "include/Utilities.h"')
    ROOT.gROOT.ProcessLine('#include "include/AnalysisTools.h"')
    #if not os.path.isdir(args.outDir):
    #    os.makedirs(args.outDir)
    if args.furtherCut:
        print(f"further cut = {args.furtherCut}")
    hist_cfg_dict = {}
    with open(args.histConfig, 'r') as f:
        hist_cfg_dict = yaml.safe_load(f)
    sample_cfg_dict = {}
    with open(args.sampleConfig, 'r') as f:
        sample_cfg_dict = yaml.safe_load(f)
    unc_cfg_dict = {}
    with open(args.uncConfig, 'r') as f:
        unc_cfg_dict = yaml.safe_load(f)
    global_cfg_dict = {}
    with open(args.globalConfig, 'r') as f:
        global_cfg_dict = yaml.safe_load(f)

    # central hist definition
    create_new_hist = False
    key_not_exist = False
    df_empty = False
    inFile_root = ROOT.TFile.Open(args.inFile,"READ")
    inFile_keys = [k.GetName() for k in inFile_root.GetListOfKeys()]
    if 'Events' not in inFile_keys:
        key_not_exist = True
    inFile_root.Close()
    if not key_not_exist and ROOT.RDataFrame('Events',args.inFile).Count().GetValue() == 0:
        df_empty = True

    create_new_hist = key_not_exist or df_empty

    if not create_new_hist:
        dfWrapped_central = DataFrameBuilder(ROOT.RDataFrame('Events',args.inFile), args.deepTauVersion)
        all_dataframes = {}
        all_histograms = {}
        sample_type = sample_cfg_dict[args.dataset]['sampleType'] if args.dataset != 'data' else 'data'
        key_central = (sample_type, "Central", "Central")
        key_filter_dict = createKeyFilterDict(global_cfg_dict)
        outfile  = ROOT.TFile(args.outFileName,'RECREATE')
        col_names_central =  dfWrapped_central.colNames
        col_tpyes_central =  dfWrapped_central.colTypes

        hasCache= args.cacheFile != ''
        if hasCache:
            dfWrapped_cache_central = DataFrameBuilder(ROOT.RDataFrame('Events',args.cacheFile), args.deepTauVersion)
            AddCacheColumnsInDf(dfWrapped_central, dfWrapped_cache_central, "cache_map_Central")

        if key_central not in all_dataframes:
            all_dataframes[key_central] = [PrepareDfWrapped(dfWrapped_central).df]
        central_histograms = GetHistogramDictFromDataframes(args.var, all_dataframes,  key_central , key_filter_dict, unc_cfg_dict['norm'],hist_cfg_dict,args.furtherCut)
        #print(central_histograms)
        # central quantities definition
        compute_variations = ( args.compute_unc_variations or args.compute_rel_weights ) and args.dataset != 'data'
        if compute_variations:
            all_dataframes[key_central][0] = createCentralQuantities(all_dataframes[key_central][0], col_tpyes_central, col_names_central)
            if all_dataframes[key_central][0].Filter("map_placeholder > 0").Count().GetValue() <= 0 : raise RuntimeError("no events passed map placeolder")

        # norm weight histograms
        if args.compute_rel_weights and args.dataset!='data':
            for uncName in unc_cfg_dict['norm'].keys():
                for scale in scales:
                    #print(uncName, scale)
                    key_2 = (sample_type, uncName, scale)
                    if key_2 not in all_dataframes.keys():
                        all_dataframes[key_2] = []
                    all_dataframes[key_2] = [all_dataframes[key_central][0]]
                    norm_histograms =  GetHistogramDictFromDataframes(args.var,all_dataframes, key_2, key_filter_dict,unc_cfg_dict['norm'], hist_cfg_dict, args.furtherCut)
                    central_histograms.update(norm_histograms)

        # save histograms
        SaveHists(central_histograms, outfile)

        # shape weight  histograms
        if args.compute_unc_variations and args.dataset!='data':
            for uncName in unc_cfg_dict['shape']:
                #print(uncName)
                for scale in scales:
                    key_2 = (sample_type, uncName, scale)
                    #print(key_2)
                    GetShapeDataFrameDict(all_dataframes, key_2, key_central, args.inFile, args.cacheFile, compute_variations, args.deepTauVersion, col_names_central, col_tpyes_central, hasCache)
                    if key_2 not in all_dataframes.keys(): continue
                    if not all_dataframes[key_2] : continue
                    shape_histograms=  GetHistogramDictFromDataframes(args.var, all_dataframes, key_2 , key_filter_dict,unc_cfg_dict['shape'], hist_cfg_dict, args.furtherCut)
                    SaveHists(shape_histograms, outfile)


        outfile.Close()
    else:
        print(f"NO HISTOGRAM CREATED!!!! dataset: {args.dataset} ")
        createVoidHist(args.outFileName, hist_cfg_dict)

    #finally:
    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))
