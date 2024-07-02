import ROOT
if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])

from Analysis.HistHelper import *
from Common.Utilities import *

deepTauYears = {'v2p1':'2017','v2p5':'2018'}
QCDregions = ['OS_Iso', 'SS_Iso', 'OS_AntiIso', 'SS_AntiIso']

#categories = ['res2b', 'res1b', 'inclusive', 'btag_shape']
categories = ['res2b', 'res1b', 'inclusive', 'boosted', 'btag_shape','baseline']
gen_channels = {'eTau':[3,5], 'muTau':[4,5], 'tauTau':[5,5]}
channels = {'eTau':13, 'muTau':23, 'tauTau':33}
triggers = {'eTau':  'HLT_singleTau || HLT_singleEle || HLT_etau || HLT_MET' , 'muTau': 'HLT_singleTau || HLT_singleMu || HLT_mutau || HLT_MET', 'tauTau' : 'HLT_singleTau || HLT_ditau || HLT_MET' }
channels = {'eTau':13, 'muTau':23, 'tauTau':33}
'''
triggers = {
    'eTau':  '(HLT_singleTau || HLT_singleEle || HLT_etau || HLT_MET)' ,
    'muTau': '(HLT_singleTau || HLT_singleMu || HLT_mutau || HLT_MET)',
    'tauTau' : '(HLT_singleTau || HLT_ditau || HLT_MET)'
}

trigger_list = {
    'eTau':  ['HLT_singleTau', 'HLT_singleEle', 'HLT_etau', 'HLT_MET'] ,
    'muTau': ['HLT_singleTau', 'HLT_singleMu', 'HLT_mutau', 'HLT_MET'],
    'tauTau' : ['HLT_singleTau' , 'HLT_ditau', 'HLT_MET']
 }
'''
triggers = {
    'eTau':  '( HLT_singleEle || HLT_etau )' ,
    'muTau': '( HLT_singleMu || HLT_mutau )',
    'tauTau' : '( HLT_ditau )'
}
trigger_list = {
    'eTau':  [ 'HLT_singleEle', 'HLT_etau'] ,
    'muTau': [ 'HLT_singleMu', 'HLT_mutau'],
    'tauTau' : [ 'HLT_ditau']
 }
btag_wps = {'res2b':'Medium', 'res1b':'Medium', 'boosted':"Loose", 'inclusive':'','btag_shape':'','baseline':''}

mass_cut_limits = {'bb_m_vis':[50,270],'tautau_m_vis':[20,130]}

scales = ['Up', 'Down']
bjet_vars = ["b1_pt","b2_pt","b1_eta","b2_eta"]
var_to_add_boosted= ["SelectedFatJet_pt_boosted","SelectedFatJet_eta_boosted"]
unc_to_not_consider_boosted = ["PUJetID", "JER","JES_FlavorQCD","JES_RelativeBal","JES_HF","JES_BBEC1","JES_EC2","JES_Absolute","JES_Total","JES_BBEC1_2018","JES_Absolute_2018","JES_EC2_2018","JES_HF_2018","JES_RelativeSample_2018","bTagSF_Loose_btagSFbc_correlated",  "bTagSF_Loose_btagSFbc_uncorrelated",  "bTagSF_Loose_btagSFlight_correlated",  "bTagSF_Loose_btagSFlight_uncorrelated",  "bTagSF_Medium_btagSFbc_correlated",  "bTagSF_Medium_btagSFbc_uncorrelated",  "bTagSF_Medium_btagSFlight_correlated",  "bTagSF_Medium_btagSFlight_uncorrelated",  "bTagSF_Tight_btagSFbc_correlated",  "bTagSF_Tight_btagSFbc_uncorrelated",  "bTagSF_Tight_btagSFlight_correlated",  "bTagSF_Tight_btagSFlight_uncorrelated","bTagShapeSF_lf","bTagShapeSF_hf","bTagShapeSF_lfstats1","bTagShapeSF_lfstats2","bTagShapeSF_hfstats1","bTagShapeSF_hfstats2","bTagShapeSF_cferr1","bTagShapeSF_cferr2"]



filters = {
        'channels':[('eTau','eTau && HLT_singleEle'), ('muTau','muTau && HLT_singleMu'),('tauTau','tauTau && HLT_ditau')],
        'QCD_regions':[('OS_Iso','OS_Iso'),('SS_Iso','SS_Iso'),('OS_AntiIso','OS_AntiIso'),('SS_AntiIso','SS_AntiIso')] ,
        'categories': [('res2b', 'res2b'), ('res1b', 'res1b'), ('inclusive', 'return true;'),('btag_shape', 'return true;')],
        }


def createKeyFilterDict(global_cfg_dict):
    reg_dict = {}
    filter_str = ""
    channels_to_consider = global_cfg_dict['channelSelection']
    qcd_regions_to_consider = global_cfg_dict['QCDRegions']
    categories_to_consider = global_cfg_dict["categories"]
    for ch in channels_to_consider:
        for reg in qcd_regions_to_consider:
            for cat in categories_to_consider:
                filter_base = f"{ch} && {triggers[ch]} && {reg} && {cat}"
                if cat =='boosted' :
                    filter_str =  filter_base
                elif cat == 'baseline':
                    filter_str = f"{filter_base}"
                else:
                    filter_str = f"b1_pt>0 && b2_pt>0 && {filter_base}"
                    for mass_name,mass_limits in mass_cut_limits.items():
                        filter_str+=f" && {mass_name} >= {mass_limits[0]} && {mass_name} <= {mass_limits[1]}"
                key = (ch, reg, cat)
                reg_dict[key] = filter_str
                #print(key, filter_str)
                #print()
    return reg_dict

def QCD_Estimation(histograms, all_samples_list, channel, category, uncName, scale, wantNegativeContributions):
    #print(channel, category)
    #print(histograms.keys())
    #key_B_data = ((channel, 'OS_AntiIso', category), ('Central', 'Central'))
    key_B = ((channel, 'OS_AntiIso', category), (uncName, scale))
    #key_C_data = ((channel, 'SS_Iso', category), ('Central', 'Central'))
    key_C = ((channel, 'SS_Iso', category), (uncName, scale))
    #key_D_data = ((channel, 'SS_AntiIso', category), ('Central', 'Central'))
    key_D = ((channel, 'SS_AntiIso', category), (uncName, scale))
    hist_data = histograms['data']
    #print(hist_data.keys())
    hist_data_B = hist_data[key_B].Clone()
    #if channel != 'tauTau' and category != 'inclusive': return hist_data_B
    hist_data_C = hist_data[key_C].Clone()
    hist_data_D = hist_data[key_D].Clone()
    n_data_B = hist_data_B.Integral(0, hist_data_B.GetNbinsX()+1)
    n_data_C = hist_data_C.Integral(0, hist_data_C.GetNbinsX()+1)
    n_data_D = hist_data_D.Integral(0, hist_data_D.GetNbinsX()+1)
    print(f"Initially Yield for data in OS AntiIso region is {key_B} is {n_data_B}")
    print(f"Initially Yield for data in SS Iso region is{key_C} is {n_data_C}")
    print(f"Initially Yield for data in SS AntiIso region is{key_D} is {n_data_D}")
    for sample in all_samples_list:
        if sample=='data' or 'GluGluToBulkGraviton' in sample or 'GluGluToRadion' in sample or 'VBFToBulkGraviton' in sample or 'VBFToRadion' in sample or sample=='QCD':
            print(f"sample {sample} is not considered")
            continue
        print(sample)
        # find kappa value
        hist_sample = histograms[sample]
        #print(histograms[sample].keys())
        hist_sample_B = hist_sample[key_B].Clone()
        hist_sample_C = hist_sample[key_C].Clone()
        hist_sample_D = hist_sample[key_D].Clone()
        n_sample_B= hist_sample_B.Integral(0, hist_sample_B.GetNbinsX()+1)
        n_data_B-=n_sample_B
        n_sample_C = hist_sample_C.Integral(0, hist_sample_C.GetNbinsX()+1)
        n_data_C-=n_sample_C
        n_sample_D = hist_sample_D.Integral(0, hist_sample_D.GetNbinsX()+1)
        n_data_D-=n_sample_D
        print(f"Yield for data in OS AntiIso region {key_B} after removing {sample} with yield {n_sample_B} is {n_data_B}")
        print(f"Yield for data in SS Iso region {key_C} after removing {sample} with yield {n_sample_C} is {n_data_C}")
        print(f"Yield for data in SS AntiIso region {key_D} after removing {sample} with yield {n_sample_D} is {n_data_D}")
        hist_data_B.Add(hist_sample_B, -1)
    #if n_data_C <= 0 or n_data_D <= 0:
        #print(f"n_data_C = {n_data_C}")
        #print(f"n_data_D = {n_data_D}")
    kappa = n_data_C/n_data_D if n_data_D != 0 else 0
    if kappa<0:
        print(f"transfer factor <0, {category}, {channel}, {uncName}, {scale}")
        return ROOT.TH1D("","",hist_data_B.GetNbinsX(), hist_data_B.GetXaxis().GetBinLowEdge(1), hist_data_B.GetXaxis().GetBinUpEdge(hist_data_B.GetNbinsX()))
        #raise  RuntimeError(f"transfer factor <=0 ! {kappa}")
    hist_data_B.Scale(kappa)
    if wantNegativeContributions:
        fix_negative_contributions,debug_info,negative_bins_info = FixNegativeContributions(hist_data_B)
        if not fix_negative_contributions:
            #return hist_data_B
            print(debug_info)
            print(negative_bins_info)
            print("Unable to estimate QCD")
            return ROOT.TH1D("","",hist_data_B.GetNbinsX(), hist_data_B.GetXaxis().GetBinLowEdge(1), hist_data_B.GetXaxis().GetBinUpEdge(hist_data_B.GetNbinsX()))

            #raise RuntimeError("Unable to estimate QCD")
    return hist_data_B


def CompareYields(histograms, all_samples_list, channel, category, uncName, scale):
    #print(channel, category)
    #print(histograms.keys())key_B_data = ((channel, 'OS_AntiIso', category), ('Central', 'Central'))
    key_A_data = ((channel, 'OS_Iso', category), ('Central', 'Central'))
    key_A = ((channel, 'OS_Iso', category), (uncName, scale))
    key_B_data = ((channel, 'OS_AntiIso', category), ('Central', 'Central'))
    key_B = ((channel, 'OS_AntiIso', category), (uncName, scale))
    key_C_data = ((channel, 'SS_Iso', category), ('Central', 'Central'))
    key_C = ((channel, 'SS_Iso', category), (uncName, scale))
    key_D_data = ((channel, 'SS_AntiIso', category), ('Central', 'Central'))
    key_D = ((channel, 'SS_AntiIso', category), (uncName, scale))
    hist_data = histograms['data']
    #print(hist_data.keys())
    hist_data_A = hist_data[key_A_data]
    hist_data_B = hist_data[key_B_data]
    #if channel != 'tauTau' and category != 'inclusive': return hist_data_B
    hist_data_C = hist_data[key_C_data]
    hist_data_D = hist_data[key_D_data]
    n_data_A = hist_data_A.Integral(0, hist_data_A.GetNbinsX()+1)
    n_data_B = hist_data_B.Integral(0, hist_data_B.GetNbinsX()+1)
    n_data_C = hist_data_C.Integral(0, hist_data_C.GetNbinsX()+1)
    n_data_D = hist_data_D.Integral(0, hist_data_D.GetNbinsX()+1)
    print(f"data || {key_A_data} || {n_data_A}")
    print(f"data || {key_B_data} || {n_data_B}")
    print(f"data || {key_C_data} || {n_data_C}")
    print(f"data || {key_D_data} || {n_data_D}")
    for sample in all_samples_list:
        #print(sample)
        # find kappa value
        hist_sample = histograms[sample]
        #print(histograms[sample].keys())
        hist_sample_A = hist_sample[key_A]
        hist_sample_B = hist_sample[key_B]
        hist_sample_C = hist_sample[key_C]
        hist_sample_D = hist_sample[key_D]
        n_sample_A = hist_sample_A.Integral(0, hist_sample_A.GetNbinsX()+1)
        n_sample_B = hist_sample_B.Integral(0, hist_sample_B.GetNbinsX()+1)
        n_sample_C = hist_sample_C.Integral(0, hist_sample_C.GetNbinsX()+1)
        n_sample_D = hist_sample_D.Integral(0, hist_sample_D.GetNbinsX()+1)

        print(f"{sample} || {key_A} || {n_sample_A}")
        print(f"{sample} || {key_B} || {n_sample_B}")
        print(f"{sample} || {key_C} || {n_sample_C}")
        print(f"{sample} || {key_D} || {n_sample_D}")

def AddQCDInHistDict(var, all_histograms, channels, categories, uncName, all_samples_list, scales, wantNegativeContributions=False):
    if 'QCD' not in all_histograms.keys():
            all_histograms['QCD'] = {}
    for channel in channels:
        for cat in categories:
            #if cat == 'btag_shape': continue
            #key =( (channel, 'OS_Iso', cat), ('Central', 'Central'))
            #all_histograms['QCD'][key] = QCD_Estimation(all_histograms, all_samples_list, channel, cat, 'Central', 'Central')
            for scale in scales + ['Central']:
                if uncName=='Central' and scale != 'Central': continue
                if uncName!='Central' and scale == 'Central': continue
                if cat == 'boosted' and var in bjet_vars: continue
                if cat == 'boosted' and uncName in unc_to_not_consider_boosted: continue
                if cat != 'boosted' and var in var_to_add_boosted: continue
                key =( (channel, 'OS_Iso', cat), (uncName, scale))
                all_histograms['QCD'][key] = QCD_Estimation(all_histograms, all_samples_list, channel, cat, uncName, scale,wantNegativeContributions)

def ApplyBTagWeight(cat,applyBtag=False, finalWeight_name = 'final_weight_0'):
    btag_weight = "1"
    btagshape_weight = "1"
    if applyBtag:
        if btag_wps[cat]!='' : btag_weight = f"weight_bTagSF_{btag_wps[cat]}_Central"
    else:
        if cat !='btag_shape' and cat !='boosted': btagshape_weight = "weight_bTagShapeSF"
    return f'{finalWeight_name}*{btag_weight}*{btagshape_weight}'


def GetWeight(channel, cat):
    weights_to_apply = ["weight_MC_Lumi_pu", "weight_L1PreFiring_Central","weight_L1PreFiring_ECAL_Central","weight_L1PreFiring_Muon_Central"]
    trg_weights_dict = {
        'eTau':["weight_tau1_TrgSF_singleEle_Central","weight_tau2_TrgSF_singleEle_Central", "weight_tau1_TrgSF_etau_Central", "weight_tau2_TrgSF_etau_Central"],
        'muTau':["weight_tau1_TrgSF_singleMu_Central","weight_tau2_TrgSF_singleMu_Central", "weight_tau1_TrgSF_mutau_Central", "weight_tau2_TrgSF_mutau_Central"],
        'tauTau':["weight_tau1_TrgSF_ditau_Central","weight_tau2_TrgSF_ditau_Central","weight_tau1_TrgSF_singleTau_Central","weight_tau2_TrgSF_singleTau_Central"]
        }
    tau_weights =["EleSF_EleIDCentral", "HighPt_MuonID_SF_HighPtIDCentral", "HighPt_MuonID_SF_HighPtIdRelTkIsoCentral", "HighPt_MuonID_SF_RecoCentral", "HighPt_MuonID_SF_TightIDCentral", "MuonID_SF_HighPtID_TrkCentral", "MuonID_SF_HighPtIdRelTkIsoCentral", "MuonID_SF_RecoCentral", "MuonID_SF_TightID_TrkCentral", "MuonID_SF_TightRelIsoCentral", "TauID_SF_Medium_Central"]
    for tau_suffix in tau_weights:
        for tau_idx in [1,2]:
            weights_to_apply.append(f"weight_tau{tau_idx}_{tau_suffix}")
    if cat != 'boosted':
         weights_to_apply.extend(["weight_Jet_PUJetID_Central_b1", "weight_Jet_PUJetID_Central_b2"])
    weights_to_apply.extend(trg_weights_dict[channel])
    total_weight = '*'.join(weights_to_apply)
    return total_weight


class DataFrameBuilder(DataFrameBuilderBase):

    def defineBoostedVariables(self):
        #FatJetObservables = ["area", "btagCSVV2", "btagDDBvLV2", "btagDeepB", "btagHbb", "deepTagMD_HbbvsQCD",
        #             "deepTagMD_ZHbbvsQCD", "deepTagMD_ZbbvsQCD", "deepTagMD_bbvsLight", "deepTag_H",
        #             "jetId", "msoftdrop", "nBHadrons", "nCHadrons",
        #             "nConstituents", "particleNetMD_QCD", "particleNetMD_Xbb", "particleNet_HbbvsQCD",
        #             "particleNet_mass", "rawFactor", "p4","pt","eta","phi","mass" ]
        FatJetObservables = ["area", "btagCSVV2", "btagDDBvLV2", "btagDeepB", "btagHbb", "deepTagMD_HbbvsQCD",
                     "deepTagMD_ZHbbvsQCD", "deepTagMD_ZbbvsQCD", "deepTagMD_bbvsLight", "deepTag_H",
                     "jetId", "msoftdrop", "nBHadrons", "nCHadrons", "nConstituents","rawFactor",
                      "particleNetMD_QCD", "particleNetMD_Xbb", "particleNet_HbbvsQCD", "particleNet_mass", # 2018
                     "particleNet_QCD","particleNet_XbbVsQCD", # 2016
                     "particleNetLegacy_QCD", "particleNetLegacy_Xbb", "particleNetLegacy_mass", # 2016
                     "particleNetWithMass_QCD", "particleNetWithMass_HbbvsQCD", "particleNet_massCorr", "p4","pt","eta","phi","mass"  # 2016
                     ]
        self.df = self.df.Define("SelectedFatJet_size_boosted","SelectedFatJet_p4[fatJet_sel].size()")
        # def the correct discriminator
        particleNet_Xbb= 'particleNetMD_Xbb' if 'SelectedFatJet_particleNetMD_Xbb' in self.df.GetColumnNames() else 'particleNetLegacy_Xbb'
        self.df = self.df.Define(f"SelectedFatJet_{particleNet_Xbb}_boosted_vec",f"SelectedFatJet_{particleNet_Xbb}[fatJet_sel]")
        self.df = self.df.Define("SelectedFatJet_idxUnordered", "CreateIndexes(SelectedFatJet_p4[fatJet_sel].size())")
        self.df = self.df.Define("SelectedFatJet_idxOrdered", f"ReorderObjects(SelectedFatJet_{particleNet_Xbb}_boosted_vec, SelectedFatJet_idxUnordered)")
        for fatJetVar in FatJetObservables:
            if f'SelectedFatJet_{fatJetVar}' in self.df.GetColumnNames() and f'SelectedFatJet_{fatJetVar}_boosted_vec' not in self.df.GetColumnNames():
                self.df = self.df.Define(f'SelectedFatJet_{fatJetVar}_boosted_vec',f""" SelectedFatJet_{fatJetVar}[fatJet_sel];""")
                self.df = self.df.Define(f'SelectedFatJet_{fatJetVar}_boosted',f"""
                                    SelectedFatJet_{fatJetVar}_boosted_vec[SelectedFatJet_idxOrdered[0]];
                                   """)
        #self.df.Display({"SelectedFatJet_p4_boosted", "SelectedFatJet_size_boosted"}).Print()

    def defineTriggers(self):
        for ch in channels:
            for trg in trigger_list[ch]:
                #print(trg)
                if trg not in self.df.GetColumnNames():
                    self.df = self.df.Define(trg, "1")

    def defineSelectionRegions(self):
        self.df = self.df.Define("nSelBtag", f"int(b1_idbtagDeepFlavB >= {self.bTagWP}) + int(b2_idbtagDeepFlavB >= {self.bTagWP})")
        particleNet_HbbvsQCD = 'particleNet_HbbvsQCD' if 'SelectedFatJet_particleNet_HbbvsQCD' in self.df.GetColumnNames() else 'particleNetWithMass_HbbvsQCD'

        self.df = self.df.Define("fatJet_presel", f"SelectedFatJet_pt>250 && SelectedFatJet_{particleNet_HbbvsQCD}>={self.pNetWP}").Define("fatJet_sel"," RemoveOverlaps(SelectedFatJet_p4, fatJet_presel, { {tau1_p4, tau2_p4},}, 2, 0.8)").Define("boosted", "SelectedFatJet_p4[fatJet_sel].size()>0")
        self.df = self.df.Define("res1b", f"!boosted && nSelBtag == 1")
        self.df = self.df.Define("res2b", f"!boosted && nSelBtag == 2")
        self.df = self.df.Define("inclusive", f"!boosted")
        self.df = self.df.Define("btag_shape", f"!boosted")
        self.df = self.df.Define("baseline",f"return true;")

    def defineChannels(self):
        for channel,ch_value in channels.items():
            self.df = self.df.Define(f"{channel}", f"channelId=={ch_value}")
            if 'genchannelId' in self.df.GetColumnNames():
                self.df = self.df.Define(f"gen_{channel}", f"genchannelId=={ch_value}")
        '''
        for gen_channel,gen_ch_value in gen_channels.items():
            if f"tau1_gen_kind" in self.df.GetColumnNames() and f"tau2_gen_kind" in self.df.GetColumnNames():
                self.df = self.df.Define(f"gen_{gen_channel}", f"tau1_gen_kind == {gen_ch_value[0]} && tau2_gen_kind == {gen_ch_value[1]}")
                print(f"gen_{gen_channel} && {gen_channel}",gen_channel, gen_ch_value, self.df.Filter(f"gen_{gen_channel} && {gen_channel}").Count().GetValue() )
            else:
                self.df = self.df.Define(f"gen_{gen_channel}", f"{gen_channel}")
        '''

    def defineL1PrefiringRelativeWeights(self):
        if "weight_L1PreFiringDown_rel" not in self.df.GetColumnNames():
            self.df = self.df.Define("weight_L1PreFiringDown_rel","weight_L1PreFiring_Down/weight_L1PreFiring_Central")
        if "weight_L1PreFiringUp_rel" not in self.df.GetColumnNames():
            self.df = self.df.Define("weight_L1PreFiringUp_rel","weight_L1PreFiringUp/weight_L1PreFiring_Central")
        if "weight_L1PreFiring_ECALDown_rel" not in self.df.GetColumnNames():
            self.df = self.df.Define("weight_L1PreFiring_ECALDown_rel","weight_L1PreFiring_ECALDown/weight_L1PreFiring_ECAL_Central")
        if "weight_L1PreFiring_Muon_StatUp_rel" not in self.df.GetColumnNames():
            self.df = self.df.Define("weight_L1PreFiring_Muon_StatUp_rel","weight_L1PreFiring_Muon_StatUp/weight_L1PreFiring_Muon_Central")
        if "weight_L1PreFiring_Muon_StatDown_rel" not in self.df.GetColumnNames():
            self.df = self.df.Define("weight_L1PreFiring_Muon_StatDown_rel","weight_L1PreFiring_Muon_StatDown/weight_L1PreFiring_Muon_Central")
        if "weight_L1PreFiring_Muon_SystUp_rel" not in self.df.GetColumnNames():
            self.df = self.df.Define("weight_L1PreFiring_Muon_SystUp_rel","weight_L1PreFiring_Muon_SystUp/weight_L1PreFiring_Muon_Central")
        if "weight_L1PreFiring_Muon_SystDown_rel" not in self.df.GetColumnNames():
            self.df = self.df.Define("weight_L1PreFiring_Muon_SystDown_rel","weight_L1PreFiring_Muon_SystDown/weight_L1PreFiring_Muon_Central")

    def defineQCDRegions(self):
        #print(self.deepTauVersion)
        tau2_iso_var = f"tau2_idDeepTau{self.deepTauYear()}{self.deepTauVersion}VSjet"
        tau1_iso_var = f"tau1_idDeepTau{self.deepTauYear()}{self.deepTauVersion}VSjet"
        self.df = self.df.Define("tau1_iso_medium", f"{tau1_iso_var} >= {Utilities.WorkingPointsTauVSjet.Medium.value}")
        self.df = self.df.Define("OS", "tau1_charge*tau2_charge < 0")
        self.df = self.df.Define("Iso", f"{tau2_iso_var} >= {Utilities.WorkingPointsTauVSjet.Medium.value}")
        self.df = self.df.Define("AntiIso", f"{tau2_iso_var} >= {Utilities.WorkingPointsTauVSjet.VVVLoose.value} && !Iso")
        self.df = self.df.Define("OS_Iso", f"tau1_iso_medium && OS && Iso")
        self.df = self.df.Define("SS_Iso", f"tau1_iso_medium && !OS && Iso")
        self.df = self.df.Define("OS_AntiIso", f"tau1_iso_medium && OS && AntiIso")
        self.df = self.df.Define("SS_AntiIso", f"tau1_iso_medium && !OS && AntiIso")

    def deepTauYear(self):
        return deepTauYears[self.deepTauVersion]

    def selectTrigger(self, trigger):
        self.df = self.df.Filter(trigger)

    def addCut (self, cut=""):
        if cut!="":
            self.df = self.df.Filter(cut)

    def ApplyMassCut(self):
        for mass_name,mass_limits in mass_cut_limits.items():
            self.df = self.df.Filter(f"{mass_name} >= {mass_limits[0]} && {mass_name} <= {mass_limits[1]}")


    def __init__(self, df, deepTauVersion='v2p1', bTagWP = 2, pNetWP = 0.9172):
        super(DataFrameBuilder, self).__init__(df)
        self.deepTauVersion = deepTauVersion
        self.bTagWP = bTagWP
        self.pNetWP = pNetWP

def PrepareDfWrapped(dfWrapped):
    dfWrapped.df = defineAllP4(dfWrapped.df)
    #dfWrapped.defineL1PrefiringRelativeWeights()
    dfWrapped.defineQCDRegions()
    dfWrapped.defineSelectionRegions()
    dfWrapped.defineBoostedVariables()
    dfWrapped.defineChannels()
    dfWrapped.defineTriggers()
    dfWrapped.df = createInvMass(dfWrapped.df)
    return dfWrapped