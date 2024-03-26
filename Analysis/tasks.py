import law
import luigi
import os
import shutil
import time
import yaml
from RunKit.sh_tools import sh_call
from RunKit.checkRootFile import checkRootFileSafe

from run_tools.law_customizations import Task, HTCondorWorkflow, copy_param, get_param_value
from AnaProd.tasks import AnaTupleTask, DataMergeTask
unc_2018 = ['JES_BBEC1_2018', 'JES_Absolute_2018', 'JES_EC2_2018', 'JES_HF_2018', 'JES_RelativeSample_2018' ]
unc_2017 = ['JES_BBEC1_2017', 'JES_Absolute_2017', 'JES_EC2_2017', 'JES_HF_2017', 'JES_RelativeSample_2017' ]
unc_2016preVFP = ['JES_BBEC1_2016preVFP', 'JES_Absolute_2016preVFP', 'JES_EC2_2016preVFP', 'JES_HF_2016preVFP', 'JES_RelativeSample_2016preVFP' ]
unc_2016postVFP = ['JES_BBEC1_2016postVFP', 'JES_Absolute_2016postVFP', 'JES_EC2_2016postVFP', 'JES_HF_2016postVFP', 'JES_RelativeSample_2016postVFP' ]

uncs_to_exclude = {
    'Run2_2018': unc_2017+ unc_2016preVFP + unc_2016postVFP,
    'Run2_2017': unc_2018+ unc_2016preVFP + unc_2016postVFP,
    'Run2_2016': unc_2017+ unc_2018 + unc_2016preVFP,
    'Run2_2016_HIPM':unc_2017+ unc_2018 + unc_2016postVFP,
    }
vars_to_plot = None
def load_vars_to_plot(hists):
    global vars_to_plot
    #vars_to_plot = list(hists.keys())# ['bbtautau_mass']
    vars_to_plot = ['kinFit_m']
    #vars_to_plot = ["tau1_pt", "tau1_eta", "tau1_phi", "tau1_mass", "tau1_idDeepTau2017v2p1VSe", "tau1_idDeepTau2017v2p1VSmu", "tau1_idDeepTau2017v2p1VSjet", "tau1_idDeepTau2018v2p5VSe", "tau1_idDeepTau2018v2p5VSmu", "tau1_idDeepTau2018v2p5VSjet", "tau1_charge", "tau1_iso"]
    #vars_to_plot += ["tau2_pt", "tau2_eta", "tau2_phi", "tau2_mass", "tau2_idDeepTau2017v2p1VSe", "tau2_idDeepTau2017v2p1VSmu", "tau2_idDeepTau2017v2p1VSjet", "tau2_idDeepTau2018v2p5VSe", "tau2_idDeepTau2018v2p5VSmu", "tau2_idDeepTau2018v2p5VSjet", "tau2_charge", "tau2_iso"]
    #vars_to_plot += ["b1_pt", "b1_eta", "b1_phi", "b1_mass", "b1_btagDeepFlavB", "b1_btagDeepFlavCvB", "b1_btagDeepFlavCvL", "b1_particleNetAK4_B", "b1_particleNetAK4_CvsB", "b1_particleNetAK4_CvsL", "b1_HHbtag"]#, "b1_hadronFlavour"]
    #vars_to_plot += ["b2_pt", "b2_eta", "b2_phi", "b2_mass", "b2_btagDeepFlavB", "b2_btagDeepFlavCvB", "b2_btagDeepFlavCvL", "b2_particleNetAK4_B", "b2_particleNetAK4_CvsB", "b2_particleNetAK4_CvsL", "b2_HHbtag"]#, "b2_hadronFlavour"]
    #vars_to_plot += ["tautau_m_vis", "bb_m_vis", "bbtautau_mass", "dR_tautau", "nBJets", "met_pt", "met_phi", "kinFit_convergence", "kinFit_result_probability", "kinFit_m", "kinFit_chi2"]
    #vars_to_plot += ["SVfit_valid", "SVfit_pt", "SVfit_pt_error", "SVfit_eta", "SVfit_eta_error", "SVfit_phi", "SVfit_phi_error", "SVfit_m", "SVfit_m_error", "SVfit_mt", "MT2"]#, "SVfit_mt_error", "MT2"]
    return vars_to_plot

hists = None
def load_hist_config(hist_config):
    global hists
    with open(hist_config, 'r') as f:
        hists = yaml.safe_load(f)
    return hists
unc_cfg_dict = None
def load_unc_config(unc_cfg):
    global unc_cfg_dict
    with open(unc_cfg, 'r') as f:
        unc_cfg_dict = yaml.safe_load(f)
    return unc_cfg_dict

def get2DOutFileName(var, sample_name, central_Histograms_path, btag_dir, suffix=''):
    outDir = os.path.join(central_Histograms_path, sample_name,var, btag_dir)
    fileName = f'{var}2D{suffix}.root'
    outFile = os.path.join(outDir,fileName)
    return outFile

def getOutFileName(var, sample_name, central_Histograms_path, btag_dir, suffix=''):
    outDir = os.path.join(central_Histograms_path, sample_name,var, btag_dir)
    fileName = f'{var}{suffix}.root'
    outFile = os.path.join(outDir,fileName)
    return outFile

# **************** 1D HISTOGRAMS *******************
class HistProducerFileTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 5.0)
    n_cpus = copy_param(HTCondorWorkflow.n_cpus, 1)
    hist_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'plot','histograms.yaml')
    hists = load_hist_config(hist_config)

    def GetBTagDir(self):
        return "bTag_WP" if self.wantBTag else "bTag_shape"

    def workflow_requires(self):
        branch_map = self.create_branch_map()
        workflow_dict = {}
        workflow_dict["anaTuple"] = {
            idx: AnaTupleTask.req(self, branch=br, branches=())
            for idx, (sample, br,var) in branch_map.items() if sample !='data' and var==vars_to_plot[0]
            for idx, (sample, br,var) in branch_map.items() if sample !='data' and var==vars_to_plot[0]
        }
        workflow_dict["dataMergeTuple"] = {
            idx: DataMergeTask.req(self, branch=0, branches=())
            for idx, (sample, br,var) in branch_map.items() if sample =='data'
        }
        return workflow_dict

    def requires(self):
        sample_name, prod_br,var = self.branch_data
        if sample_name =='data':
            return [ DataMergeTask.req(self,max_runtime=DataMergeTask.max_runtime._default, branch=prod_br, branches=())]
        return [ AnaTupleTask.req(self, branch=prod_br, max_runtime=AnaTupleTask.max_runtime._default, branches=())]

    def create_branch_map(self):
        n = 0
        k = 0
        branches = {}
        #codeprint(vars_to_plot)
        anaProd_branch_map = AnaTupleTask.req(self, branch=-1, branches=()).create_branch_map()
        #vars_to_plot = ['bbtautau_mass']#list(hists.keys())# ['tau1_pt','tau2_pt','tau1_eta','tau2_eta','SelectedFatJet_pt_boosted','SelectedFatJet_eta_boosted'] #
        #vars_to_plot = ['bbtautau_mass']#['tau1_pt','tau2_pt','tau1_eta','tau2_eta','SelectedFatJet_pt_boosted','SelectedFatJet_eta_boosted'] #
        for var in vars_to_plot:
            for prod_br, (sample_id, sample_name, sample_type, input_file) in anaProd_branch_map.items():
                if sample_type=='data' or sample_type=='VBFToRadion' or sample_type=='VBFToBulkGraviton' or 'QCD' in sample_type : continue
                branches[n] = (sample_name, prod_br,var)
                n+=1
            branches[n] = ('data', 0, var)
            n+=1
        return branches


    def output(self):
        sample_name, prod_br,var = self.branch_data
        input_file = self.input()[0].path
        fileName  = os.path.basename(input_file)
        #vars_to_plot = list(hists.keys())
        #local_files_target = []
        #outFile_histProdSample = get2DOutFileName(var, sample_name, self.central_Histograms_path(), self.GetBTagDir(), '_onlyCentral')
        #if os.path.exists(outFile_histProdSample):
        #    return law.LocalFileTarget(outFile_histProdSample)
        outDir = os.path.join(self.central_Histograms_path(), sample_name, 'tmp', var, self.GetBTagDir())
        if not os.path.isdir(outDir):
            os.makedirs(outDir)
        outFile = os.path.join(outDir,fileName)
        return law.LocalFileTarget(outFile)

    def run(self):
        sample_name, prod_br,var = self.branch_data
        if len(self.input()) > 1:
            raise RuntimeError(f"multple input files!! {' '.join(f.path for f in self.input())}")
        input_file = self.input()[0].path
        outFileName = os.path.basename(input_file)
        outFileName_split = os.path.basename(input_file).split('.')[0]
        hist_config = self.hist_config
        unc_config = os.path.join(self.ana_path(), 'config', 'weight_definition.yaml')
        sample_config = self.sample_config
        outDir = os.path.join(self.central_Histograms_path(), sample_name, 'tmp')
        cacheFile = f'/eos/home-k/kandroso/cms-hh-bbtautau/anaCache/Run2_2018/{sample_name}/{self.version}/{outFileName}'
        HistProducerFile = os.path.join(self.ana_path(), 'Analysis', 'HistProducerFile.py')
        HistProducerFile_cmd = ['python3', HistProducerFile,'--inFile', input_file, '--cacheFile', cacheFile, '--outFileName',outFileName_split, '--dataset', sample_name, '--outDir', outDir, '--compute_unc_variations', 'True', '--compute_rel_weights', 'True', '--uncConfig', unc_config, '--histConfig', hist_config,'--sampleConfig', sample_config, '--var', var]
        if self.version.split('_')[1] == 'deepTau2p5':
            HistProducerFile_cmd.extend([ '--deepTauVersion', 'v2p5'])
        if self.wantBTag:
            HistProducerFile_cmd.extend([ '--wantBTag', f'{self.wantBTag}'])
        #print(HistProducerFile_cmd)
        sh_call(HistProducerFile_cmd,verbose=1)


class HistProducerSampleTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 1.0)
    hist_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'plot','histograms.yaml')
    hists = load_hist_config(hist_config)

    def GetBTagDir(self):
        return "bTag_WP" if self.wantBTag else "bTag_shape"


    def workflow_requires(self):
        self_map = self.create_branch_map()
        workflow_dict = {}
        workflow_dict["histProducerFile"] = {
            n: HistProducerFileTask.req(self, branches=tuple((idx,) for idx in idx_list))
            for n, (sample_name, idx_list,var) in self_map.items() if var==vars_to_plot[0]
        }
        return workflow_dict

    def requires(self):
        sample_name, idx_list,var  = self.branch_data
        deps = [HistProducerFileTask.req(self, max_runtime=HistProducerFileTask.max_runtime._default, branch=idx) for idx in idx_list ]
        return deps


    def create_branch_map(self):
        branches = {}
        histProducerFile_map = HistProducerFileTask.req(self,branch=-1, branches=()).create_branch_map()
        all_samples = {}
        for br_idx, (sample_name, prod_br,var) in histProducerFile_map.items():
            if sample_name not in all_samples:
                all_samples[sample_name] = {}
            if var not in all_samples[sample_name]:
                all_samples[sample_name][var]=[]
            all_samples[sample_name][var].append(br_idx)
        k=0
        for n, key in enumerate(all_samples.items()):
            (smpl_name, dictionary)=key
            if ('QCD' in smpl_name) or ('VBF' in smpl_name and ('Radion' in smpl_name or 'BulkGraviton' in smpl_name)): continue
            for m, (var, idx_list) in enumerate(dictionary.items()):
                if var == vars_to_plot[0]:
                    branches[k] = (smpl_name, idx_list, var)
                    k+=1
        #print(branches)
        return branches

    def output(self):
        sample_name, idx_list,var  = self.branch_data
        local_files_target = []
        #vars_to_plot = ["tau1_pt", "tau1_eta", "tau1_phi", "tau1_mass", "tau1_idDeepTau2017v2p1VSe", "tau1_idDeepTau2017v2p1VSmu", "tau1_idDeepTau2017v2p1VSjet", "tau1_idDeepTau2018v2p5VSe", "tau1_idDeepTau2018v2p5VSmu", "tau1_idDeepTau2018v2p5VSjet", "tau1_charge", "tau1_iso"]
        for obs in vars_to_plot:
            outFile = getOutFileName(obs, sample_name, self.central_Histograms_path(), self.GetBTagDir())
            local_files_target.append(law.LocalFileTarget(outFile))
        return local_files_target #law.LocalFileTarget(outFile)

    def run(self):
        sample_name, idx_list,var  = self.branch_data
        #print(histProducerFile_map)
        hists_str = ','.join(var for var in vars_to_plot)
        files_idx = []
        file_ids_str = ''
        file_name_pattern = 'nano'
        if(len(idx_list)>1):
            file_name_pattern +="_{id}"
            file_ids_str = f"0-{len(idx_list)}"
        file_name_pattern += ".root"
        HistProducerSample = os.path.join(self.ana_path(), 'Analysis', 'HistProducerSample.py')
        outDir = os.path.join(self.central_Histograms_path(), sample_name)
        histDir = os.path.join(self.central_Histograms_path(), sample_name, 'tmp')
        HistProducerSample_cmd = ['python3', HistProducerSample,'--histDir', histDir, '--outDir', outDir, '--hists', hists_str, '--file-name-pattern', file_name_pattern, '--remove-files', 'True']
        if self.wantBTag:
            HistProducerSample_cmd.extend([ '--wantBTag', f'{self.wantBTag}'])
        if(len(idx_list)>1):
            HistProducerSample_cmd.extend(['--file-ids', file_ids_str])
        #print(HistProducerSample_cmd)
        #print(f"outFile is {self.output().path}")
        sh_call(HistProducerSample_cmd,verbose=1)


class MergeTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 30.0)
    hist_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'plot','histograms.yaml')
    hists = load_hist_config(hist_config)
    unc_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'weight_definition.yaml')
    unc_cfg_dict = load_unc_config(unc_config)

    def GetBTagDir(self):
        return "bTag_WP" if self.wantBTag else "bTag_shape"


    def workflow_requires(self):
        prod_branches = HistProducerSampleTask.req(self, branch=-1, branches=()).create_branch_map()
        workflow_dict = {}
        workflow_dict["histProducerSample"] = {
            0: HistProducerSampleTask.req(self, branches=tuple((br,) for br in prod_branches))
        }
        return workflow_dict

    def requires(self):
        prod_branches = HistProducerSampleTask.req(self, branch=-1, branches=()).create_branch_map()
        deps = [HistProducerSampleTask.req(self, max_runtime=HistProducerSampleTask.max_runtime._default, branch=prod_br) for prod_br in prod_branches ]
        return deps


    def create_branch_map(self):
        uncNames = ['Central']
        uncNames.extend(list(unc_cfg_dict['norm'].keys()))
        uncNames.extend([unc for unc in unc_cfg_dict['shape']])
        #vars_to_plot = list(hists.keys())
        n = 0
        branches = {}
        #vars_to_plot = ["tau1_pt", "tau1_eta", "tau1_phi", "tau1_mass", "tau1_idDeepTau2017v2p1VSe", "tau1_idDeepTau2017v2p1VSmu", "tau1_idDeepTau2017v2p1VSjet", "tau1_charge", "tau1_iso"]
        #vars_to_plot += ["tau2_pt", "tau2_eta", "tau2_phi", "tau2_mass", "tau2_idDeepTau2017v2p1VSe", "tau2_idDeepTau2017v2p1VSmu", "tau2_idDeepTau2017v2p1VSjet", "tau2_charge", "tau2_iso"]
        for var in vars_to_plot :
            for uncName in uncNames:
                if uncName in uncs_to_exclude[self.period]: continue
                branches[n] = (var, uncName)
                n += 1
        return branches

    def output(self):
        var, uncName = self.branch_data
        #outFile_haddMergedFiles = os.path.join(self.central_Histograms_path(), 'all_histograms',var,self.GetBTagDir(),f'all_histograms_{var}.root')
        #if os.path.exists(outFile_haddMergedFiles):
        #    return law.LocalFileTarget((outFile_haddMergedFiles))
        local_file_target = os.path.join(self.central_Histograms_path(), 'all_histograms',var,self.GetBTagDir(),f'all_histograms_{var}_{uncName}.root')
        return law.LocalFileTarget(local_file_target)

    def run(self):
        var, uncName = self.branch_data
        #vars_to_plot = list(hists.keys())
        print(f"outFile is {self.output().path}")
        sample_config = self.sample_config
        unc_config = os.path.join(self.ana_path(), 'config', 'weight_definition.yaml')
        json_dir= os.path.join(self.valeos_path(), 'jsonFiles', self.period, self.version)
        MergerProducer = os.path.join(self.ana_path(), 'Analysis', 'HistMerger.py')
        MergerProducer_cmd = ['python3', MergerProducer,'--histDir', self.central_Histograms_path(), '--sampleConfig', sample_config, '--var', var,'--uncConfig', unc_config, '--uncSource', uncName,'--jsonDir', json_dir]
        if self.wantBTag:
            MergerProducer_cmd.extend([ '--wantBTag', f'{self.wantBTag}'])
        sh_call(MergerProducer_cmd,verbose=1)
        #print(MergerProducer_cmd)




class HistRebinnerTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 10.0)
    hist_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'plot','histograms.yaml')
    hists = load_hist_config(hist_config)
    unc_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'weight_definition.yaml')
    unc_cfg_dict = load_unc_config(unc_config)

    def GetBTagDir(self):
        return "bTag_WP" if self.wantBTag else "bTag_shape"


    def workflow_requires(self):
        workflow_dict = {}
        workflow_dict["MergeTask"] = {
            0: MergeTask.req(self,  branch=-1, branches=())
        }
        return workflow_dict

    def requires(self):
        deps = [MergeTask.req(self, branch=-1, branches=()) ]
        return deps


    def create_branch_map(self):
        uncNames = ['Central']
        uncNames.extend(list(unc_cfg_dict['norm'].keys()))
        uncNames.extend([unc for unc in unc_cfg_dict['shape']])
        #vars_to_plot = list(hists.keys())
        n = 0
        branches = {}
        #vars_to_plot = ["tau1_pt", "tau1_eta", "tau1_phi", "tau1_mass", "tau1_idDeepTau2017v2p1VSe", "tau1_idDeepTau2017v2p1VSmu", "tau1_idDeepTau2017v2p1VSjet", "tau1_charge", "tau1_iso"]
        #vars_to_plot += ["tau2_pt", "tau2_eta", "tau2_phi", "tau2_mass", "tau2_idDeepTau2017v2p1VSe", "tau2_idDeepTau2017v2p1VSmu", "tau2_idDeepTau2017v2p1VSjet", "tau2_charge", "tau2_iso"]
        for var in vars_to_plot :
            for uncName in uncNames:
                if uncName in uncs_to_exclude[self.period]: continue
                branches[n] = (var, uncName)
                n += 1
        return branches

    def output(self):
        var, uncName = self.branch_data
        local_file_target = os.path.join(self.central_Histograms_path(), 'all_histograms',var,self.GetBTagDir(),f'all_histograms_{var}_{uncName}_Rebinned.root')
        return law.LocalFileTarget(local_file_target)

    #python3 /afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/Analysis/HistRebinner.py --histDir /eos/home-v/vdamante/HH_bbtautau_resonant_Run2/histograms/Run2_2018/v9_deepTau2p1 --inFileName all_histograms --sampleConfig /afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/config/samples_Run2_2018.yaml --var kinFit_m --uncConfig /afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/config/weight_definition.yaml --histConfig /afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/config/plot/histograms.yaml --uncSource TauID_stat1_DM0

    def run(self):
        var, uncName = self.branch_data
        print(f"outFile is {self.output().path}")
        sample_config = self.sample_config
        unc_config = os.path.join(self.ana_path(), 'config', 'weight_definition.yaml')
        json_dir= os.path.join(self.valeos_path(), 'jsonFiles', self.period, self.version)
        RebinnerProducer = os.path.join(self.ana_path(), 'Analysis', 'HistRebinner.py')
        RebinnerProducer_cmd = ['python3', RebinnerProducer,'--histDir', self.central_Histograms_path(), '--inFileName', 'all_histograms', '--sampleConfig', sample_config, '--var', var,'--uncConfig', unc_config, '--histConfig', self.hist_config, '--uncSource', uncName]
        if self.wantBTag:
            RebinnerProducer_cmd.extend([ '--wantBTag', f'{self.wantBTag}'])
        sh_call(RebinnerProducer_cmd,verbose=1)
        #print(MergerProducer_cmd)


# ************* CENTRAL ONLY **********************

class HistSampleTaskCentral(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 6.0)
    n_cpus = copy_param(HTCondorWorkflow.n_cpus, 1)
    hist_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'plot','histograms.yaml')
    hists = load_hist_config(hist_config)
    vars_to_plot = load_vars_to_plot(hists)
    vars_to_plot = load_vars_to_plot(hists)

    def GetBTagDir(self):
        return "bTag_WP" if self.wantBTag else "bTag_shape"

    def workflow_requires(self):
        branch_map = self.create_branch_map()
        workflow_dict = {}
        workflow_dict["anaTuple"] = {
            idx: AnaTupleTask.req(self, branches=tuple((br,) for br in branches))
            for idx, (sample, branches,var) in branch_map.items() if sample !='data' and var==vars_to_plot[0]
            for idx, (sample, branches,var) in branch_map.items() if sample !='data' and var==vars_to_plot[0]
        }
        workflow_dict["dataMergeTuple"] = {
            idx: DataMergeTask.req(self, branch=0, branches=())
            for idx, (sample, branches,var) in branch_map.items() if sample =='data' and var==vars_to_plot[0]
        }
        return workflow_dict

    def requires(self):
        sample, prod_branches,var = self.branch_data
        if sample =='data':
            return [ DataMergeTask.req(self,max_runtime=DataMergeTask.max_runtime._default, branch=0, branches=())]
        return [AnaTupleTask.req(self, max_runtime=AnaTupleTask.max_runtime._default, branch=prod_br) for prod_br in prod_branches ]

    def create_branch_map(self):
        n=0
        branches = {}
        anaProd_branch_map = AnaTupleTask.req(self, branch=-1, branches=()).create_branch_map()
        for var in vars_to_plot:
            all_samples = {}
            for prod_br, (sample_id, sample_name, sample_type, input_file) in anaProd_branch_map.items():
                if sample_type=='data' or sample_type=='VBFToRadion' or sample_type=='VBFToBulkGraviton' or 'QCD' in sample_type : continue
                if sample_name not in all_samples:
                    all_samples[sample_name] = {}
                if var not in all_samples[sample_name]:
                    all_samples[sample_name][var] = []
                all_samples[sample_name][var].append(prod_br)
            for sample_name,all_samples_sample in all_samples.items():
                for var, idx_list in all_samples_sample.items():
                    branches[n] = (sample_name, idx_list, var)
                n+=1
            branches[n] = ('data',[0], var)
            n+=1
        return branches


    def output(self):
        sample_name, idx_list,var = self.branch_data
        #print(sample_name)
        #print(sample_name)
        outDir = os.path.join(self.central_Histograms_path(), sample_name, var, self.GetBTagDir())
        fileName = f'{var}_onlyCentral.root'
        outFile = os.path.join(outDir, fileName)
        return law.LocalFileTarget(outFile)

    def run(self):
        sample_name, prod_branches,var = self.branch_data
        fileName = f'{var}_onlyCentral.root'
        outFile = os.path.join(self.central_Histograms_path(), sample_name, var, self.GetBTagDir(), fileName)
        sample_config = self.sample_config
        CheckFile = os.path.join(self.ana_path(), 'Analysis', 'checkFile.py')
        CheckFile_cmd = ['python3', CheckFile,'--inFile', outFile, '--sampleConfig', sample_config]

        #print(os.path.exists(outFile))
        if (os.path.exists(outFile)):
            print(f"{outFile} exists")
            sh_call(CheckFile_cmd,verbose=1)
        while(not os.path.exists(outFile)):
            print(f"{outFile} does not exist")

            inDir = os.path.join(self.central_anaTuples_path(), sample_name)
            hist_config = self.hist_config
            outDir = os.path.join(self.central_Histograms_path(), sample_name)
            cacheDir = f'/eos/home-k/kandroso/cms-hh-bbtautau/anaCache/Run2_2018/{sample_name}/{self.version}'
            HistProducerFile = os.path.join(self.ana_path(), 'Analysis', 'HistSampleCentral.py')
            HistProducerFile_cmd = ['python3', HistProducerFile,'--inDir', inDir, '--cacheDir', cacheDir, '--outDir', outDir, '--dataset', sample_name, '--histConfig', hist_config,'--sampleConfig', sample_config, '--var', var]
            if self.version.split('_')[1] == 'deepTau2p5':
                HistProducerFile_cmd.extend([ '--deepTauVersion', 'v2p5'])
            if self.wantBTag:
                HistProducerFile_cmd.extend([ '--wantBTag', f'{self.wantBTag}'])
            sh_call(HistProducerFile_cmd,verbose=1)
            sh_call(CheckFile_cmd,verbose=1)

class MergeTaskCentral(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 30.0)
    hist_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'plot','histograms.yaml')
    hists = load_hist_config(hist_config)
    unc_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'weight_definition.yaml')
    unc_cfg_dict = load_unc_config(unc_config)
    vars_to_plot = load_vars_to_plot(hists)
    vars_to_plot = load_vars_to_plot(hists)

    def GetBTagDir(self):
        return "bTag_WP" if self.wantBTag else "bTag_shape"

    def workflow_requires(self):
        workflow_dict = {}
        workflow_dict["histProducerSample"] = {
            0: HistSampleTaskCentral.req(self, branches=())
        }
        return workflow_dict

    def requires(self):
        deps = [HistSampleTaskCentral.req(self, branches=()) ]
        return deps

    def create_branch_map(self):
        n = 0
        branches = {}
        for var in vars_to_plot :
            branches[n] = var
            n += 1
        return branches

    def output(self):
        var = self.branch_data
        outDir = os.path.join(self.central_Histograms_path(), 'all_histograms',var,self.GetBTagDir())
        outFile = f'all_histograms_{var}_onlyCentral.root'
        json_dir= os.path.join(self.valeos_path(), 'jsonFiles', self.period, self.version, 'all_ratios')
        json_file = f'all_ratios_{var}_onlyCentral.json'
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        if not os.path.exists(json_dir):
            os.makedirs(json_dir)
        local_files_target = []
        local_files_target.append(law.LocalFileTarget(os.path.join(outDir,outFile)))
        local_files_target.append(law.LocalFileTarget(os.path.join(json_dir,json_file)))
        return local_files_target

    def run(self):
        var = self.branch_data
        unc_config = os.path.join(self.ana_path(), 'config', 'weight_definition.yaml')
        sample_config = self.sample_config
        unc_config = os.path.join(self.ana_path(), 'config', 'weight_definition.yaml')
        MergerProducer = os.path.join(self.ana_path(), 'Analysis', 'HistMerger.py')
        json_dir= os.path.join(self.valeos_path(), 'jsonFiles', self.period, self.version)
        MergerProducer_cmd = ['python3', MergerProducer,'--histDir', self.central_Histograms_path(), '--sampleConfig', sample_config, '--var', var,'--uncConfig', unc_config, '--jsonDir', json_dir, '--suffix','_onlyCentral']
        if self.wantBTag:
            MergerProducer_cmd.extend(['--wantBTag','True' ])
        sh_call(MergerProducer_cmd,verbose=1)


class PlotterTaskCentral(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 10.0)
    hist_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'plot','histograms.yaml')
    hists = load_hist_config(hist_config)
    vars_to_plot = load_vars_to_plot(hists)

    def GetBTagDir(self):
        return "bTag_WP" if self.wantBTag else "bTag_shape"


    def workflow_requires(self):
        return { "MergeTaskCentral" : MergeTaskCentral.req(self, branches=()) }

    def requires(self):
        return [MergeTaskCentral.req(self, branches=())]

    def create_branch_map(self):
        sample_cfg_dict = self.sample_config
        categories = ['res1b', 'res2b', 'inclusive', 'baseline']#list(self.global_params['categories'])
        QCD_region = 'OS_Iso'
        uncName = 'Central'
        channels = list(self.global_params['channelSelection'])
        n = 0
        branches = {}
        for var in vars_to_plot:
            for channel in channels:
                for category in categories:
                    branches[n] = (var,channel, QCD_region, category, uncName)
                    n+=1
        #print(branches)
        print(self.version)
        print(self.version.split('_')[1] == 'deepTau2p5')
        return branches

    def output(self):
        var,channel, QCD_region, category, uncName = self.branch_data

        plotdirname = 'plots_2p5' if self.version.split('_')[1] == 'deepTau2p5' else 'plots_2p1'
        eosdir = '/eos/home-v/vdamante/www'
        final_directory = os.path.join(eosdir, plotdirname, var, self.GetBTagDir(), 'from1D')
        #final_directory = os.path.join(self.central_Histograms_path(), 'plots', var, self.GetBTagDir())
        final_fileName = f'{var}_{channel}_{category}_XMass{self.mass}.pdf'
        #final_fileName = f'{var}_{channel}_{QCD_region}_{category}_XMass{self.mass}.pdf'
        final_file= os.path.join(final_directory, final_fileName)
        return law.LocalFileTarget(final_file)

    def run(self):
        var,channel, QCD_region, category, uncName = self.branch_data
        unc_config = os.path.join(self.ana_path(), 'config', 'weight_definition.yaml')
        sample_config = self.sample_config
        PlotterProducer = os.path.join(self.ana_path(), 'Analysis', 'HistPlotter.py')
        plotdirname = 'plots_2p5' if self.version.split('_')[1] == 'deepTau2p5' else 'plots_2p1'
        eosdir = '/eos/home-v/vdamante/www/'
        outDir = os.path.join(eosdir, plotdirname)
        PlotterProducer_cmd = ['python3', PlotterProducer,'--mass', self.mass, '--histDir', self.central_Histograms_path(), '--outDir',outDir, '--inFileName', 'all_histograms', '--var', var, '--sampleConfig',sample_config, '--channel', channel, '--category', category, '--uncSource', uncName, '--suffix','_onlyCentral']
        sh_call(PlotterProducer_cmd,verbose=1)

# ***************************************

class HaddMergedTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 1.0)
    hist_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'plot','histograms.yaml')
    hists = load_hist_config(hist_config)
    vars_to_plot=load_vars_to_plot(hists)

    def GetBTagDir(self):
        return "bTag_WP" if self.wantBTag else "bTag_shape"


    def workflow_requires(self):
        return { "Merge" : MergeTask.req(self, branches=()) }

    def requires(self):
        return [MergeTask.req(self, branches=())]

    def create_branch_map(self):
        #vars_to_plot = list(hists.keys())
        n = 0
        branches = {}
        for var in vars_to_plot :
            branches[n] = var
            n += 1
        return branches

    def output(self):
        var = self.branch_data
        #local_file_target = os.path.join(self.central_Histograms_path(), 'all_histograms',var,self.GetBTagDir(),f'all_histograms_{var}.root')
        local_file_target = os.path.join(self.central_Histograms_path(), 'all_histograms',var,self.GetBTagDir(),f'all_histograms_{var}_Rebinned.root')
        return law.LocalFileTarget(local_file_target)

    def run(self):
        print(f"outFile is {self.output().path}")
        var = self.branch_data
        unc_config = os.path.join(self.ana_path(), 'config', 'weight_definition.yaml')
        outDir_json = os.path.join(self.valeos_path(), 'jsonFiles', self.period, self.version,'fitResults')
        histDir = self.central_Histograms_path()
        ShapeOrLogNormalProducer = os.path.join(self.ana_path(), 'Analysis', 'ShapeOrLogNormal.py')
        ShapeOrLogNormal_cmd = ['python3', ShapeOrLogNormalProducer,'--mass', self.mass, '--histDir', histDir, '--inFileName', 'all_histograms', '--outDir', outDir_json, '--sampleConfig',self.sample_config,'--uncConfig', unc_config,'--var', var, '--suffix','_Rebinned']
        sh_call(ShapeOrLogNormal_cmd,verbose=1)
        final_json_file = os.path.join(outDir_json, var, self.GetBTagDir(), 'slopeInfo.json')
        HaddMergedHistsProducer = os.path.join(self.ana_path(), 'Analysis', 'hadd_merged_hists.py')
        HaddMergedHistsProducer_cmd = ['python3', HaddMergedHistsProducer,'--histDir', histDir, '--file-name-pattern', 'all_histograms', '--hists', var, '--uncConfig', unc_config, '--suffix','_Rebinned', '--period', self.period]
        #if os.path.exists(final_json_file):
            #HaddMergedHistsProducer_cmd.extend(['--remove-files', 'True'])
        sh_call(HaddMergedHistsProducer_cmd,verbose=1)

