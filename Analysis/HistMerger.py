import ROOT
import sys
import os
import math
import shutil
import json
import time
from RunKit.run_tools import ps_call
if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])

import Common.Utilities as Utilities
from Analysis.HistHelper import *
from Analysis.hh_bbtautau import *



def checkLists(list1, list2):
    if len(list1) != len(list2):
        print(f"lists have different length: {list1} and {list2}")
        return False
    for item in list1:
        if item not in list2:
            print(f"{item} in {list1} but not in {list2}")
            return False
    return True

def checkFile(inFileRoot, channels, qcdRegions, categories, var):
    keys_channel = [str(key.GetName()) for key in inFileRoot.GetListOfKeys()]
    if not (checkLists(keys_channel, channels)):
        print("check list not worked for channels")
        return False
    for channel in channels:
        dir_0 = inFileRoot.Get(channel)
        keys_qcdRegions = [str(key.GetName()) for key in dir_0.GetListOfKeys()]
        if not checkLists(keys_qcdRegions, QCDregions):
            print("check list not worked for QCDregions")
            return False
        for qcdRegion in QCDregions:
            dir_1 = dir_0.Get(qcdRegion)
            keys_categories = [str(key.GetName()) for key in dir_1.GetListOfKeys()]
            if var in bjet_vars and 'boosted' not in keys_categories:
                keys_categories.extend(['boosted'])
            if not checkLists(keys_categories, categories):
                    print("check list not worked for categories")
                    return False
            for cat in categories:
                if cat == 'boosted' and var in bjet_vars: continue
                if cat != 'boosted' and var in var_to_add_boosted: continue
                #print(cat, var)
                dir_2 = dir_1.Get(cat)
                keys_histograms = [str(key.GetName()) for key in dir_2.GetListOfKeys()]
                #print(keys_histograms)
                if not keys_histograms: return False
    return True


def getHistDict(var, all_histograms, inFileRoot,channels, QCDregions, categories, uncSource,sample_name):
    #print("STARTING GET HIST DICT")
    #print(sample_name)
    for channel in channels:
        #if sample_name=='data' and uncSource!='Central': break
        dir_0 = inFileRoot.Get(channel)
        #print(dir_0.GetListOfKeys())
        for qcdRegion in QCDregions:
            dir_1 = dir_0.Get(qcdRegion)
            #print(dir_1.GetListOfKeys())
            for cat in categories:
                if cat == 'boosted' and var in bjet_vars: continue
                if cat != 'boosted' and var in var_to_add_boosted: continue
                if cat == 'boosted' and uncSource in unc_to_not_consider_boosted: continue
                #print(cat, var)
                dir_2 = dir_1.Get(cat)
                for key in dir_2.GetListOfKeys():
                    obj = key.ReadObj()
                    key_name = key.GetName()
                    key_name_split = key_name.split('_')
                    obj.SetDirectory(0)
                    if not obj.IsA().InheritsFrom(ROOT.TH1.Class()): continue

                    sample_type = key_name_split[0]
                    #if sample_name == 'data': print(key_name_split)
                    name_to_use = sample_name
                    if sample_type in sample_types_to_merge:
                        name_to_use = sample_type
                    #if sample_name == 'data': print(f"NAME TO USE IS {name_to_use}")
                    if name_to_use not in all_histograms.keys():
                        all_histograms[name_to_use] = {}
                    key_total = ((channel, qcdRegion, cat), ('Central', 'Central'))
                    #if sample_name == 'data': print(f"UNC SOURCE IS {uncSource}")
                    if uncSource=='Central':
                        key_total = ((channel, qcdRegion, cat), ('Central', 'Central'))
                        if len(key_name_split)>1:continue
                        if key_total not in all_histograms[name_to_use].keys():
                            all_histograms[name_to_use][key_total] = []
                        all_histograms[name_to_use][key_total].append(obj)
                    else:
                        uncName_scale ='' if name_to_use == 'data' else '_'.join(n for n in key_name_split[1:])
                        for scale in ['Up','Down']:
                            if name_to_use != 'data' and uncSource+scale != uncName_scale: continue
                            #print(uncName_scale)
                            key_total = ((channel, qcdRegion, cat), (uncSource, scale))
                            if key_total not in all_histograms[name_to_use].keys():
                                all_histograms[name_to_use][key_total] = []
                            all_histograms[name_to_use][key_total].append(obj)
                    #if sample_name == 'data': print(f"KEY TOTAL FOR {name_to_use} IS {key_total}")
                    #if sample_name == 'data': print(f"LEN OF ALL HISTOGRAMS IS {len(all_histograms[name_to_use][key_total])}")

def MergeHistogramsPerType(all_histograms):
    for sample_type in all_histograms.keys():
        if sample_type == 'data': print(f"DURING MERGE HISTOGRAMS, sample_type is {sample_type}")
        for key_name,histlist in all_histograms[sample_type].items():
            final_hist =  histlist[0]
            objsToMerge = ROOT.TList()
            for hist in histlist[1:]:
                objsToMerge.Add(hist)
            final_hist.Merge(objsToMerge)
            all_histograms[sample_type][key_name] = final_hist
            if len(histlist)!=1:
                print(f"for {sample_type} the lenght of histlist is {len(histlist)}")
            #if sample_type=='data' :
            #    print(key_name, len(histlist))
            #    print("Integral")
            #    print(key_name, final_hist.Integral(0,final_hist.GetNbinsX()+1))



def GetBTagWeightDict(var, all_histograms):
    all_histograms_1D = {}
    final_json_dict = {}
    for sample_type in all_histograms.keys():
        #print(sample_type)
        all_histograms_1D[sample_type] = {}
        final_json_dict[sample_type]={}
        for key_name,histogram in all_histograms[sample_type].items():
            (key_1, key_2) = key_name
            ch, reg, cat = key_1
            uncName,scale = key_2
            key_tuple_num = ((ch, reg, 'btag_shape'), key_2)
            key_tuple_den = ((ch, reg, 'inclusive'), key_2)
            ratio_num_hist = all_histograms[sample_type][key_tuple_num] if key_tuple_num in all_histograms[sample_type].keys() else None
            ratio_den_hist = all_histograms[sample_type][key_tuple_den] if key_tuple_den in all_histograms[sample_type].keys() else None
            histlist =[]
            second_key = f'{ch}_{reg}_{cat}_{uncName}{scale}'
            if second_key not in final_json_dict[sample_type].keys():
                final_json_dict[sample_type][second_key]={}
            if 'yBin' not in final_json_dict[sample_type][second_key].keys():
                final_json_dict[sample_type][second_key]['yBin']=[]
            if 'nJets' not in final_json_dict[sample_type][second_key].keys():
                final_json_dict[sample_type][second_key]['nJets']=[]
            if 'num' not in final_json_dict[sample_type][second_key].keys():
                final_json_dict[sample_type][second_key]['num']=[]
            if 'den' not in final_json_dict[sample_type][second_key].keys():
                final_json_dict[sample_type][second_key]['den']=[]
            if 'ratio' not in final_json_dict[sample_type][second_key].keys():
                final_json_dict[sample_type][second_key]['ratio']=[]

            num = ratio_num_hist.Integral(0,ratio_num_hist.GetNbinsX()+1)
            den = ratio_den_hist.Integral(0,ratio_den_hist.GetNbinsX()+1)
            ratio = 0.
            if ratio_den_hist.Integral(0,ratio_den_hist.GetNbinsX()+1) != 0 :
                ratio = ratio_num_hist.Integral(0,ratio_num_hist.GetNbinsX()+1)/ratio_den_hist.Integral(0,ratio_den_hist.GetNbinsX()+1)
            final_json_dict[sample_type][second_key]['num'].append(num)
            final_json_dict[sample_type][second_key]['den'].append(den)
            final_json_dict[sample_type][second_key]['ratio'].append(ratio)
            #if ratio == 0 and hist1D.Integral(0, hist1D.GetNbinsX()+1) ==0 :
            #    continue
            histogram.Scale(ratio)
            all_histograms_1D[sample_type][key_name] = histogram
    return all_histograms_1D,final_json_dict


if __name__ == "__main__":
    import argparse
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('inputFile', nargs='+', type=str)
    #parser.add_argument('datasetFile', nargs='+', type=str)
    parser.add_argument('--outFile', required=True, type=str)
    parser.add_argument('--jsonFile', required=True, type=str)
    parser.add_argument('--year', required=True, type=str)
    parser.add_argument('--datasetFile', required=True, type=str)
    parser.add_argument('--var', required=True, type=str)
    parser.add_argument('--sampleConfig', required=True, type=str)
    parser.add_argument('--globalConfig', required=True, type=str)
    parser.add_argument('--uncConfig', required=True, type=str)
    parser.add_argument('--uncSource', required=False, type=str,default='Central')

    args = parser.parse_args()
    startTime = time.time()
    with open(args.uncConfig, 'r') as f:
        unc_cfg_dict = yaml.safe_load(f)
    with open(args.sampleConfig, 'r') as f:
        sample_cfg_dict = yaml.safe_load(f)
    with open(args.globalConfig, 'r') as f:
        global_cfg_dict = yaml.safe_load(f)

    all_samples_list = args.datasetFile.split(',')
    all_samples_types = {}
    #print(all_samples_list)
    #print(all_samples_types)
    uncNameTypes = GetUncNameTypes(unc_cfg_dict)

    if args.uncSource != 'Central' and args.uncSource not in uncNameTypes:
        print("unknown unc source {args.uncSource}")

    categories = list(global_cfg_dict['categories'])
    QCDregions = list(global_cfg_dict['QCDRegions'])
    channels = list(global_cfg_dict['channelSelection'])
    signals = list(global_cfg_dict['signal_types'])
    files_separated = {}
    all_histograms ={}
    #all_histograms_1D ={}
    #final_json_dict = {}
    all_infiles = [ fileName for fileName in args.inputFile ]
    #print(all_samples_list)
    #print(len(all_infiles), len(all_samples_list))
    if len(all_infiles) != len(all_samples_list):
        raise RuntimeError(f"all_infiles have len {len(all_infiles)} and all_samples_list have len {len(all_samples_list)}")
    #print()
    #print("STARTING LOOP ON FILES/SAMPLES")
    #print()
    ignore_samples = []
    for (inFileName, sample_name) in zip(all_infiles, all_samples_list):
        #print("BEGINNING OF LOOP")
        #print("INFILE NAME")
        #print(inFileName)
        #print("SAMPLE NAME")
        #print(sample_name)
        #print()
        if not os.path.exists(inFileName):
            print(f"{inFileName} does not exist")
            continue
            #raise RuntimeError(f"{inFileName} removed")
        #print(sample_name)
        inFileRoot = ROOT.TFile.Open(inFileName, "READ")
        if inFileRoot.IsZombie():
            inFileRoot.Close()
            os.remove(inFileName)
            ignore_samples.append(sample_name)
            raise RuntimeError(f"{inFileName} is Zombie")
        if  not checkFile(inFileRoot, channels, QCDregions, categories,args.var):
            print(f"{sample_name} has void file")
            ignore_samples.append(sample_name)
            #all_samples_list.remove(sample_name)
            inFileRoot.Close()
            continue
            #os.remove(inFileName)
            #print(all_infiles.index(fileName))
            #raise RuntimeError(f"{inFileName} has problems")
        #print("AFTER CHECK FILES")
        #print(inFileName)
        #print(sample_name)
        if sample_name == 'data':
            print("BEFORE GET HIST DICT DEFINING SAMPLES_TYPES")
            print(sample_name)
        getHistDict(args.var,all_histograms, inFileRoot,channels, QCDregions, categories, args.uncSource,sample_name)
        #print(all_histograms)
        if sample_name == 'data':
            all_samples_types['data'] = ['data']
        else:
            sample_type=sample_cfg_dict[sample_name]['sampleType']
            sample_key = sample_type if sample_type in sample_types_to_merge else sample_name
            if sample_name not in ignore_samples:
                if sample_key not in all_samples_types.keys(): all_samples_types[sample_key] = []
                all_samples_types[sample_key].append(sample_name)
        inFileRoot.Close()
        #if sample_name == 'data':
        #    print("AFTER DEFINING SAMPLES_TYPES")
        print(sample_name)
    #print()
    #print("STARTING MERGE HISTOGRAMS PER TYPE")
    #print()
    #print('data' in all_histograms.keys())
    MergeHistogramsPerType(all_histograms)
    all_histograms_1D,final_json_dict=GetBTagWeightDict(args.var,all_histograms)
    fixNegativeContributions=True
    if args.var == 'kinFit_m':
        fixNegativeContributions=False

    AddQCDInHistDict(args.var,all_histograms_1D, channels, categories, args.uncSource, all_samples_types.keys(), scales,fixNegativeContributions)

    with open(f"{args.jsonFile}.json", "w") as write_file:
        json.dump(final_json_dict, write_file, indent=4)

    outFile = ROOT.TFile(args.outFile, "RECREATE")

    for sample_type in all_histograms_1D.keys():
        #print(sample_type)
        for key in all_histograms_1D[sample_type]:
            (channel, qcdRegion, cat), (uncNameType, uncScale) = key
            if qcdRegion != 'OS_Iso': continue
            dirStruct = (channel, cat)
            dir_name = '/'.join(dirStruct)
            dir_ptr = mkdir(outFile,dir_name)
            hist = all_histograms_1D[sample_type][key]
            #print(sample_type, key, hist.GetEntries())
            hist_name =  sample_type
            if uncNameType!=args.uncSource: continue
            if uncNameType != 'Central':
                if sample_type == 'data' : continue
                if uncScale == 'Central': continue
                hist_name+=f"_{uncNameType}_{uncScale}"
            else:
                if uncScale!='Central':continue
            hist.SetTitle(hist_name)
            hist.SetName(hist_name)
            dir_ptr.WriteTObject(hist, hist_name, "Overwrite")
    outFile.Close()
    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))