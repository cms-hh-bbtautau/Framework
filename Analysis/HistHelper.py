import sys
import math
import ROOT
import os
if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])

import Common.Utilities as Utilities

scales = ['Up','Down']
unc_2018 = ['JES_BBEC1_2018', 'JES_Absolute_2018', 'JES_EC2_2018', 'JES_HF_2018', 'JES_RelativeSample_2018' ]
unc_2017 = ['JES_BBEC1_2017', 'JES_Absolute_2017', 'JES_EC2_2017', 'JES_HF_2017', 'JES_RelativeSample_2017' ]
unc_2016preVFP = ['JES_BBEC1_2016preVFP', 'JES_Absolute_2016preVFP', 'JES_EC2_2016preVFP', 'JES_HF_2016preVFP', 'JES_RelativeSample_2016preVFP' ]
unc_2016postVFP = ['JES_BBEC1_2016postVFP', 'JES_Absolute_2016postVFP', 'JES_EC2_2016postVFP', 'JES_HF_2016postVFP', 'JES_RelativeSample_2016postVFP' ]
sample_types_to_merge = ['DY','TT','ZJNuNu','ZQQ','W', 'ttH']

uncs_to_exclude = {
    'Run2_2018': unc_2017+ unc_2016preVFP + unc_2016postVFP,
    'Run2_2017': unc_2018+ unc_2016preVFP + unc_2016postVFP,
    'Run2_2016': unc_2017+ unc_2018 + unc_2016preVFP,
    'Run2_2016_HIPM':unc_2017+ unc_2018 + unc_2016postVFP,
    }

def GetUncNameTypes(unc_cfg_dict):
    uncNames = []
    uncNames.extend(list(unc_cfg_dict['norm'].keys()))
    uncNames.extend([unc for unc in unc_cfg_dict['shape']])
    return uncNames

def GetSamplesStuff(sample_cfg_dict,global_cfg_dict,wantSignals=True,wantAllMasses=True,wantOneMass=True,mass=500):
    all_samples_list = []
    all_samples_types = {'data':['data'],}
    signals = list(global_cfg_dict['signal_types'])
    for sample in sample_cfg_dict.keys():
        if sample == 'GLOBAL' : continue
        if 'sampleType' not in sample_cfg_dict[sample].keys(): continue
        sample_type = sample_cfg_dict[sample]['sampleType']
        if sample_type in ['QCD', 'VBFToRadion','VBFToBulkGraviton', 'data']: continue
        if sample != 'GluGluToHHTo2B2Tau_node_SM' and sample_type == 'HHnonRes': continue
        if sample == "W0JetsToLNu-amcatnloFXFX": continue
        if sample == "W1JetsToLNu-amcatnloFXFX": continue
        if sample == "W2JetsToLNu-amcatnloFXFX": continue
        sample_name = sample_type if sample_type in ['DY','TT','ZJNuNu','ZQQ','W'] else sample
        if wantOneMass:
            if 'mass' in sample_cfg_dict[sample].keys():
                if sample_type in signals and sample_cfg_dict[sample]['mass']!=mass : continue
        if not wantAllMasses and not wantOneMass:
            if 'mass' in sample_cfg_dict[sample].keys(): continue
        isSignal = False
        if sample_type in signals:
            isSignal = True
            if not wantSignals: continue
            sample_name=sample
        if sample_name not in all_samples_types.keys() :
            all_samples_types[sample_name] = []
        all_samples_types[sample_name].append(sample)
        if not wantAllMasses and isSignal: continue
        if sample_type in all_samples_list: continue
        all_samples_list.append(sample_name)
        #print(sample_type, sample_name, all_samples_types)
    print(all_samples_types)
    return all_samples_list, all_samples_types


def CreateNamesDict(histNamesDict, sample_types, uncName, scales, sample_cfg_dict,global_cfg_dict):
    signals = list(global_cfg_dict['signal_types'])
    for sample_key in sample_types.keys():
        final_sampleKey=f"{sample_key}"
        if sample_key == 'data':
            histNamesDict[final_sampleKey] = (sample_key, 'Central','Central')
            continue
        else:
            if uncName == 'Central':
                histNamesDict[final_sampleKey] = (sample_key, 'Central','Central')
                continue
            else:
                for scale in scales:
                    histName = f"{final_sampleKey}_{uncName}{scale}"
                    histKey = (sample_key,  uncName, scale)
                    histNamesDict[histName] = histKey



def defineP4(df, name):
    df = df.Define(f"{name}_p4", f"ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>({name}_pt,{name}_eta,{name}_phi,{name}_mass)")
    return df

def defineAllP4(df):
    df = df.Define(f"SelectedFatJet_idx", f"CreateIndexes(SelectedFatJet_pt.size())")
    df = df.Define(f"SelectedFatJet_p4", f"GetP4(SelectedFatJet_pt, SelectedFatJet_eta, SelectedFatJet_phi, SelectedFatJet_mass, SelectedFatJet_idx)")
    for idx in [0,1]:
        df = defineP4(df, f"tau{idx+1}")
        df = defineP4(df, f"b{idx+1}")
    return df

def createInvMass(df):
    particleNet_mass = 'particleNet_mass' if 'SelectedFatJet_particleNet_mass_boosted' in df.GetColumnNames() else 'particleNetLegacy_mass'
    df = df.Define("tautau_m_vis", "static_cast<float>((tau1_p4+tau2_p4).M())")
    df = df.Define("bb_m_vis", f"""
                   if (!boosted){{
                       return static_cast<float>((b1_p4+b2_p4).M());
                       }}
                    return static_cast<float>(SelectedFatJet_{particleNet_mass}_boosted);""")
    df = df.Define("bbtautau_mass", """
                   if (!boosted){
                       return static_cast<float>((b1_p4+b2_p4+tau1_p4+tau2_p4).M());
                       }
                    return static_cast<float>((SelectedFatJet_p4_boosted+tau1_p4+tau2_p4).M());""")
    df = df.Define("dR_tautau", 'ROOT::Math::VectorUtil::DeltaR(tau1_p4, tau2_p4)')
    return df

def RenormalizeHistogram(histogram, norm, include_overflows=True):
    integral = histogram.Integral(0, histogram.GetNbinsX()+1) if include_overflows else histogram.Integral()
    if integral!=0:
        histogram.Scale(norm / integral)

def FixNegativeContributions(histogram):
    correction_factor = 0.

    ss_debug = ""
    ss_negative = ""

    original_Integral = histogram.Integral(0, histogram.GetNbinsX()+1)
    ss_debug += "\nSubtracted hist for '{}'.\n".format(histogram.GetName())
    ss_debug += "Integral after bkg subtraction: {}.\n".format(original_Integral)
    if original_Integral < 0:
        print(ss_debug)
        print("Integral after bkg subtraction is negative for histogram '{}'".format(histogram.GetName()))
        return False,ss_debug, ss_negative

    for n in range(1, histogram.GetNbinsX()+1):
        if histogram.GetBinContent(n) >= 0:
            continue
        prefix = "WARNING" if histogram.GetBinContent(n) + histogram.GetBinError(n) >= 0 else "ERROR"

        ss_negative += "{}: {} Bin {}, content = {}, error = {}, bin limits=[{},{}].\n".format(
            prefix, histogram.GetName(), n, histogram.GetBinContent(n), histogram.GetBinError(n),
            histogram.GetBinLowEdge(n), histogram.GetBinLowEdge(n+1))

        error = correction_factor - histogram.GetBinContent(n)
        new_error = math.sqrt(math.pow(error, 2) + math.pow(histogram.GetBinError(n), 2))
        histogram.SetBinContent(n, correction_factor)
        histogram.SetBinError(n, new_error)

    RenormalizeHistogram(histogram, original_Integral, True)
    return True, ss_debug, ss_negative

class DataFrameBuilderBase:
    def CreateColumnTypes(self):
        #colNames = [str(c) for c in self.df.GetColumnNames() if 'kinFit_result' not in str(c)]
        colNames = [str(c) for c in self.df.GetColumnNames()]#if 'kinFit_result' not in str(c)]
        entryIndexIdx = colNames.index("entryIndex")
        runIdx = colNames.index("run")
        eventIdx = colNames.index("event")
        lumiIdx = colNames.index("luminosityBlock")
        colNames[entryIndexIdx], colNames[0] = colNames[0], colNames[entryIndexIdx]
        colNames[runIdx], colNames[1] = colNames[1], colNames[runIdx]
        colNames[eventIdx], colNames[2] = colNames[2], colNames[eventIdx]
        colNames[lumiIdx], colNames[3] = colNames[3], colNames[lumiIdx]
        self.colNames = colNames
        #if "kinFit_result" in self.colNames:
        #    self.colNames.remove("kinFit_result")
        self.colTypes = [str(self.df.GetColumnType(c)) for c in self.colNames]

    def __init__(self, df):
        self.df = df
        self.colNames=[]
        self.colTypes=[]
        self.var_list = []
        self.CreateColumnTypes()

    def CreateFromDelta(self,central_columns,central_col_types):
        var_list =[]
        for var_idx,var_name in enumerate(self.colNames):
            if not var_name.endswith("Diff"):
                continue
            var_name_forDelta = var_name.removesuffix("Diff")
            central_col_idx = central_columns.index(var_name_forDelta)
            if central_columns[central_col_idx]!=var_name_forDelta:
                raise RuntimeError(f"CreateFromDelta: {central_columns[central_col_idx]} != {var_name_forDelta}")
            self.df = self.df.Define(f"{var_name_forDelta}", f"""analysis::FromDelta({var_name},
                                     analysis::GetEntriesMap()[std::make_tuple(entryIndex, run, event, luminosityBlock)]->GetValue<{self.colTypes[var_idx]}>({central_col_idx}) )""")
            var_list.append(f"{var_name_forDelta}")
        for central_col_idx,central_col in enumerate(central_columns):
            if central_col in var_list or central_col in self.colNames: continue
            self.df = self.df.Define(central_col, f"""analysis::GetEntriesMap()[std::make_tuple(entryIndex, run, event, luminosityBlock)]->GetValue<{central_col_types[central_col_idx]}>({central_col_idx})""")

    def AddCacheColumns(self,cache_cols,cache_col_types):
        for cache_col_idx,cache_col in enumerate(cache_cols):
            if  cache_col in self.df.GetColumnNames(): continue
            if cache_col.replace('.','_') in self.df.GetColumnNames(): continue
            self.df = self.df.Define(cache_col.replace('.','_'), f"""analysis::GetCacheEntriesMap().at(std::make_tuple(entryIndex, run, event, luminosityBlock))->GetValue<{cache_col_types[cache_col_idx]}>({cache_col_idx})""")

def GetModel(hist_cfg, var):
    x_bins = hist_cfg[var]['x_bins']
    if type(hist_cfg[var]['x_bins'])==list:
        x_bins_vec = Utilities.ListToVector(x_bins, "double")
        model = ROOT.RDF.TH1DModel("", "", x_bins_vec.size()-1, x_bins_vec.data())
    else:
        n_bins, bin_range = x_bins.split('|')
        start,stop = bin_range.split(':')
        model = ROOT.RDF.TH1DModel("", "",int(n_bins), float(start), float(stop))
    return model


def Get2DModel(hist_cfg, var):
    x_bins = hist_cfg[var]['x_bins']
    if type(hist_cfg[var]['x_bins'])==list:
        x_bins_vec = Utilities.ListToVector(x_bins, "double")
        model = ROOT.RDF.TH2DModel("", "", x_bins_vec.size()-1, x_bins_vec.data(), 13, -0.5, 12.5)
    else:
        n_bins, bin_range = x_bins.split('|')
        start,stop = bin_range.split(':')
        model = ROOT.RDF.TH2DModel("", "",int(n_bins), float(start), float(stop), 13, -0.5, 12.5)
    return model


def mkdir(file, path):
    dir_names = path.split('/')
    current_dir = file
    for n, dir_name in enumerate(dir_names):
        dir_obj = current_dir.Get(dir_name)
        full_name = f'{file.GetPath()}' + '/'.join(dir_names[:n])
        if dir_obj:
            if not dir_obj.IsA().InheritsFrom(ROOT.TDirectory.Class()):
                raise RuntimeError(f'{dir_name} already exists in {full_name} and it is not a directory')
        else:
            dir_obj = current_dir.mkdir(dir_name)
            if not dir_obj:

                raise RuntimeError(f'Failed to create {dir_name} in {full_name}')
        current_dir = dir_obj
    return current_dir




