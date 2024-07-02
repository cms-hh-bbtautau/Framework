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


processes= ["DY","EWK_WminusToLNu","EWK_WplusToLNu","EWK_ZTo2L","EWK_ZTo2Nu","GluGluHToTauTau_M125","GluGluHToWWTo2L2Nu_M125","GluGluToBulkGravitonToHHTo2B2Tau_M-1000","GluGluToBulkGravitonToHHTo2B2Tau_M-1250","GluGluToBulkGravitonToHHTo2B2Tau_M-1500","GluGluToBulkGravitonToHHTo2B2Tau_M-1750","GluGluToBulkGravitonToHHTo2B2Tau_M-2000","GluGluToBulkGravitonToHHTo2B2Tau_M-250","GluGluToBulkGravitonToHHTo2B2Tau_M-2500","GluGluToBulkGravitonToHHTo2B2Tau_M-260","GluGluToBulkGravitonToHHTo2B2Tau_M-270","GluGluToBulkGravitonToHHTo2B2Tau_M-280","GluGluToBulkGravitonToHHTo2B2Tau_M-300","GluGluToBulkGravitonToHHTo2B2Tau_M-3000","GluGluToBulkGravitonToHHTo2B2Tau_M-320","GluGluToBulkGravitonToHHTo2B2Tau_M-350","GluGluToBulkGravitonToHHTo2B2Tau_M-400","GluGluToBulkGravitonToHHTo2B2Tau_M-450","GluGluToBulkGravitonToHHTo2B2Tau_M-500","GluGluToBulkGravitonToHHTo2B2Tau_M-550","GluGluToBulkGravitonToHHTo2B2Tau_M-600","GluGluToBulkGravitonToHHTo2B2Tau_M-650","GluGluToBulkGravitonToHHTo2B2Tau_M-700","GluGluToBulkGravitonToHHTo2B2Tau_M-750","GluGluToBulkGravitonToHHTo2B2Tau_M-800","GluGluToBulkGravitonToHHTo2B2Tau_M-850","GluGluToBulkGravitonToHHTo2B2Tau_M-900","GluGluToRadionToHHTo2B2Tau_M-1000","GluGluToRadionToHHTo2B2Tau_M-1250","GluGluToRadionToHHTo2B2Tau_M-1500","GluGluToRadionToHHTo2B2Tau_M-1750","GluGluToRadionToHHTo2B2Tau_M-2000","GluGluToRadionToHHTo2B2Tau_M-250","GluGluToRadionToHHTo2B2Tau_M-2500","GluGluToRadionToHHTo2B2Tau_M-260","GluGluToRadionToHHTo2B2Tau_M-270","GluGluToRadionToHHTo2B2Tau_M-280","GluGluToRadionToHHTo2B2Tau_M-300","GluGluToRadionToHHTo2B2Tau_M-3000","GluGluToRadionToHHTo2B2Tau_M-320","GluGluToRadionToHHTo2B2Tau_M-350","GluGluToRadionToHHTo2B2Tau_M-400","GluGluToRadionToHHTo2B2Tau_M-450","GluGluToRadionToHHTo2B2Tau_M-500","GluGluToRadionToHHTo2B2Tau_M-550","GluGluToRadionToHHTo2B2Tau_M-600","GluGluToRadionToHHTo2B2Tau_M-650","GluGluToRadionToHHTo2B2Tau_M-700","GluGluToRadionToHHTo2B2Tau_M-750","GluGluToRadionToHHTo2B2Tau_M-800","GluGluToRadionToHHTo2B2Tau_M-850","GluGluToRadionToHHTo2B2Tau_M-900","GluGluZH_HToWW_ZTo2L_M125","HWminusJ_HToWW_M125","HWplusJ_HToWW_M125","HZJ_HToWW_M125","ST_t-channel_antitop_4f_InclusiveDecays","ST_t-channel_top_4f_InclusiveDecays","ST_tW_antitop_5f_InclusiveDecays","ST_tW_top_5f_InclusiveDecays","TT4b","TTGG","TTGamma_Dilept","TTGamma_Hadronic","TTGamma_SingleLept","TTTJ","TTTT","TTTW","TT","TTWH","TTWJetsToLNu","TTWJetsToQQ","TTWW","TTWZ","TTZH","TTZHTo4b","TTZHToNon4b","TTZToLLNuNu_M-10","TTZToNuNu","TTZToQQ","TTZZ","TTZZTo4b","TTZZToNon4b","VBFHToTauTau_M125","VBFHToWWTo2L2Nu_M125","W","WW","WWW_4F","WWZ_4F","WZ","WZZ","WminusHToTauTau_M125","WplusHToTauTau_M125","ZHToTauTau_M125","ZH_Hbb_Zll","ZH_Hbb_Zqq","ZJNuNu","ZZ","ZZZ","ggHToZZTo2L2Q","ttH","data","QCD"]

uncReNames = {
    "bTagShapeSF_lf":"CMS_btag_LF",
    "bTagShapeSF_hf":"CMS_btag_HF",
    "bTagShapeSF_lfstats1":"CMS_btag_lfstats1_{}",
    "bTagShapeSF_lfstats2":"CMS_btag_lfstats2_{}",
    "bTagShapeSF_hfstats1":"CMS_btag_hfstats1_{}",
    "bTagShapeSF_hfstats2":"CMS_btag_hfstats2_{}",
    "bTagShapeSF_cferr1":"CMS_btag_cferr1",
    "bTagShapeSF_cferr2":"CMS_btag_cferr2",
    "EleID":"CMS_eff_e_{}",
    "PUJetID":"CMS_eff_j_PUJET_id_{}",
    "HighPtMuon_HighPtID":"CMS_eff_m_highpt_id_HighPtID_{}",
    "HighPtMuon_HighPtIDIso":"CMS_eff_m_highpt_id_HighPtIDIso_{}",
    "HighPtMuon_Reco":"CMS_eff_m_highpt_id_Reco_{}",
    "HighPtMuon_TightID":"CMS_eff_m_highpt_id_TightID_{}",
    "Muon_HighPtID":"CMS_eff_m_id_highptid{}",
    "MuonID_HighPtIso":"CMS_eff_m_id_HighPtIso_{}",
    "MuonID_Reco":"CMS_eff_m_id_Reco_{}",
    "MuonID_TightID":"CMS_eff_m_id_TightID_{}",
    "MuonID_TightIDIso":"CMS_eff_m_id_TightIDIso_{}",
    "TauID_total":"CMS_eff_t_id_total_{}", # da capire
    "TauID_stat1_DM0":"CMS_eff_t_id_stat1_DM0_{}",
    "TauID_stat1_DM10":"CMS_eff_t_id_stat1_DM10_{}",
    "TauID_stat1_DM11":"CMS_eff_t_id_stat1_DM11_{}",
    "TauID_stat1_DM1":"CMS_eff_t_id_stat1_DM1_{}",
    "TauID_stat2_DM0":"CMS_eff_t_id_stat2_DM0_{}",
    "TauID_stat2_DM10":"CMS_eff_t_id_stat2_DM10_{}",
    "TauID_stat2_DM11":"CMS_eff_t_id_stat2_DM11_{}",
    "TauID_stat2_DM1":"CMS_eff_t_id_stat2_DM1_{}",
    "TauID_syst_alleras":"CMS_eff_t_id_syst_alleras",
    "TauID_syst_year":"CMS_eff_t_id_syst_{}",
    "TauID_syst_year_DM0":"CMS_eff_t_id_syst_{}_DM0",
    "TauID_syst_year_DM10":"CMS_eff_t_id_syst_{}_DM10",
    "TauID_syst_year_DM11":"CMS_eff_t_id_syst_{}_DM11",
    "TauID_syst_year_DM1":"CMS_eff_t_id_syst_{}_DM1",
    "TauID_genuineElectron_barrel":"CMS_eff_t_id_etauFR_barrel_{}",
    "TauID_genuineElectron_endcaps":"CMS_eff_t_id_etauFR_endcaps_{}",
    "TauID_genuineMuon_eta0p4to0p8":"CMS_eff_t_id_mutauFR_eta0p4to0p8_{}",
    "TauID_genuineMuon_eta0p8to1p2":"CMS_eff_t_id_mutauFR_eta0p8to1p2_{}",
    "TauID_genuineMuon_eta1p2to1p7":"CMS_eff_t_id_mutauFR_eta1p2to1p7_{}",
    "TauID_genuineMuon_etaGt1p7":"CMS_eff_t_id_mutauFR_etaGt1p7_{}",
    "TauID_genuineMuon_etaLt0p4":"CMS_eff_t_id_mutauFR_etaLt0p4_{}",
    "TauID_stat_highpT_bin1":"CMS_eff_t_id_stat_highpT_bin1_{}",
    "TauID_stat_highpT_bin2":"CMS_eff_t_id_stat_highpT_bin2_{}",
    "TauID_syst_highpT_bin1":"CMS_eff_t_id_syst_highpT_bin1",
    "TauID_syst_highpT_bin2":"CMS_eff_t_id_syst_highpT_bin2",
    "TauID_syst_highpT_extrap":"CMS_eff_t_id_syst_highpT_extrap",
    "L1Prefiring":"CMS_l1_prefiring_{}",
    "L1Prefiring_ECAL":"CMS_l1_ecal_prefiring_{}",
    "L1Prefiring_Muon_Stat":"CMS_l1_muon_stat_prefiring_{}",
    "L1Prefiring_Muon_Syst":"CMS_l1_muon_syst_prefiring_{}",
    "PileUp_Lumi_MC":"CMS_pileup_{}",
    "TauES_DM0": "CMS_scale_t_DM0_{}",
    "TauES_DM1": "CMS_scale_t_DM1_{}",
    "TauES_3prong": "CMS_scale_t_3prong_{}",
    "EleFakingTauES_DM0": "CMS_scale_t_eFake_DM0_{}",
    "EleFakingTauES_DM1": "CMS_scale_t_eFake_DM1_{}",
    "MuFakingTauES": "CMS_scale_t_muFake_{}",
    "JER": "CMS_res_j_{}",
    "JES_Absolute": "CMS_scale_j_Abs",
    "JES_Absolute_2018": "CMS_scale_j_Abs_2018",
    "JES_Absolute_2017": "CMS_scale_j_Abs_2017",
    "JES_Absolute_2016preVFP": "CMS_scale_j_Abs_2016_HIPM" ,
    "JES_Absolute_2016postVFP": "CMS_scale_j_Abs_2016" ,
    "JES_BBEC1": "CMS_scale_j_BBEC1",
    "JES_BBEC1_2018": "CMS_scale_j_BBEC1_2018",
    "JES_BBEC1_2017": "CMS_scale_j_BBEC1_2017",
    "JES_BBEC1_2016preVFP": "CMS_scale_j_BBEC1_2016_HIPM" ,
    "JES_BBEC1_2016postVFP": "CMS_scale_j_BBEC1_2016" ,
    "JES_EC2": "CMS_scale_j_EC2",
    "JES_EC2_2018": "CMS_scale_j_EC2_2018",
    "JES_EC2_2017": "CMS_scale_j_EC2_2017",
    "JES_EC2_2016preVFP": "CMS_scale_j_EC2_2016_HIPM" ,
    "JES_EC2_2016postVFP": "CMS_scale_j_EC2_2016" ,
    "JES_FlavorQCD": "CMS_scale_j_FlavQCD",
    "JES_HF": "CMS_scale_j_HF",
    "JES_HF_2018": "CMS_scale_j_HF_2018",
    "JES_HF_2017": "CMS_scale_j_HF_2017",
    "JES_HF_2016preVFP": "CMS_scale_j_HF_2016_HIPM" ,
    "JES_HF_2016postVFP": "CMS_scale_j_HF_2016" ,
    "JES_RelativeBal": "CMS_scale_j_RelBal",
    "JES_RelativeSample_2018": "CMS_scale_j_RelSample_2018",
    "JES_RelativeSample_2017": "CMS_scale_j_RelSample_2017",
    "JES_RelativeSample_2016preVFP": "CMS_scale_j_RelSample_2016_HIPM" ,
    "JES_RelativeSample_2016postVFP": "CMS_scale_j_RelSample_2016" ,
    "TrigSF_diTau_DM0":"CMS_bbtt_trigdiTau_DM0_{}",
    "TrigSF_diTau_DM1":"CMS_bbtt_trigdiTau_DM1_{}",
    "TrigSF_diTau_3Prong":"CMS_bbtt_trigdiTau_3Prong_{}",
    "TrigSF_singleEle":"CMS_bbtt_trigsingleEle_{}",
    "TrigSF_singleMu":"CMS_bbtt_trigsingleMu_{}",
    "TrigSF_eTau_Ele":"CMS_bbtt_trigeTau_Ele_{}",
    "TrigSF_eTau_DM0":"CMS_bbtt_trigeTau_DM0_{}",
    "TrigSF_eTau_DM1":"CMS_bbtt_trigeTau_DM1_{}",
    "TrigSF_eTau_3Prong":"CMS_bbtt_trigeTau_3Prong_{}",
    "TrigSF_muTau_Mu":"CMS_bbtt_trigmuTau_Mu_{}",
    "TrigSF_muTau_DM0":"CMS_bbtt_trigmuTau_DM0_{}",
    "TrigSF_muTau_DM1":"CMS_bbtt_trigmuTau_DM1_{}",
    "TrigSF_muTau_3Prong":"CMS_bbtt_trigmuTau_3Prong_{}",
    "TrgSF_etau_DM0":"CMS_bbtt_trgetau_DM0_{}",
    "TrgSF_etau_DM1":"CMS_bbtt_trgetau_DM1_{}",
    "TrgSF_etau_ele":"CMS_bbtt_trgetau_ele_{}",
    "TrgSF_mutau_3Prong":"CMS_bbtt_trgmutau_3Prong_{}",
    "TrgSF_mutau_DM0":"CMS_bbtt_trgmutau_DM0_{}",
    "TrgSF_mutau_DM1":"CMS_bbtt_trgmutau_DM1_{}",
    "TrgSF_mutau_mu":"CMS_bbtt_trgmutau_mu_{}",
    "TrgSF_singleEle":"CMS_bbtt_trgsingleEle_{}",
    "TrgSF_singleMu24":"CMS_bbtt_trgsingleMu24_{}",
}


if __name__ == "__main__":
    import argparse
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('--inFile', required=True)
    parser.add_argument('--outFile', required=True)
    parser.add_argument('--year', required=True)
    parser.add_argument('--var', required=False, type=str, default='tau1_pt')
    #parser.add_argument('--remove-files', required=False, type=bool, default=False)
    args = parser.parse_args()


all_histnames = {}

for process in processes:
    all_histnames[process] = process
    for unc_old in uncReNames.keys():
        new_unc = uncReNames[unc_old].format(args.year)
        #print(unc_old.split('_'))
        #print(args.year)
        for scale in ['Up','Down']:
            #if args.year in unc_old.split('_'):
            all_histnames[f"{process}_{unc_old}_{args.year}{scale}"] = f"{process}_{new_unc}{scale}"
            all_histnames[f"{process}_{unc_old}_{args.year}_{scale}"] = f"{process}_{new_unc}{scale}"
            #else:
            all_histnames[f"{process}_{unc_old}{scale}"] = f"{process}_{new_unc}{scale}"
            all_histnames[f"{process}_{unc_old}_{scale}"] = f"{process}_{new_unc}{scale}"

#print(all_histnames.keys())
inFile = ROOT.TFile.Open(args.inFile, "READ")
outFile = ROOT.TFile.Open(args.outFile, "RECREATE")
channels =[str(key.GetName()) for key in inFile.GetListOfKeys()]
for channel in channels:
    dir_0 = inFile.Get(channel)
    keys_categories = [str(key.GetName()) for key in dir_0.GetListOfKeys()]
    for cat in keys_categories:
        dir_1= dir_0.Get(cat)
        for key_hist in dir_1.GetListOfKeys():
            key_name = key_hist.GetName()
            if key_name not in all_histnames.keys():
                print(f"{key_name} not in all_histnames keys")
                continue
            #key_hist = dir_1.Get(key_name)
            obj = key_hist.ReadObj()
            obj.SetDirectory(0)
            dirStruct = (channel, cat)
            dir_name = '/'.join(dirStruct)
            dir_ptr = mkdir(outFile,dir_name)
            obj.SetTitle(all_histnames[key_name])
            obj.SetName(all_histnames[key_name])
            dir_ptr.WriteTObject(obj, all_histnames[key_name], "Overwrite")
outFile.Close()
inFile.Close()

