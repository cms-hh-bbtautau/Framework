import law
import luigi
import os
import shutil
import time
import yaml
import contextlib
import ROOT


from RunKit.run_tools import ps_call, natural_sort
from RunKit.crabLaw import cond as kInit_cond, update_kinit_thread
from run_tools.law_customizations import Task, HTCondorWorkflow, copy_param
from Common.Utilities import SerializeObjectToString
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
    #### Kinematics #####
    vars_to_plot = ["kinFit_m"]
    #vars_to_plot += ["tau1_pt", "tau1_eta","tau2_pt", "tau2_eta"]
    #vars_to_plot += ["b1_pt", "b1_eta","b2_pt","b2_eta"]
    #vars_to_plot += ["b1_phi", "b1_mass","b2_phi","b2_mass"]
    #vars_to_plot += ["tau1_phi", "tau1_mass","tau2_phi","tau2_mass"]
    #vars_to_plot += ["tau1_charge","tau1_iso","tau2_charge","tau2_iso"]


    #### Global Obs #####
    #vars_to_plot += ["tautau_m_vis", "bb_m_vis", "bbtautau_mass", "dR_tautau", "nBJets"]
    #vars_to_plot += ["met_pt", "met_phi"]

    #### Vars ANACACHE ####
    # vars_to_plot = ["kinFit_m","kinFit_convergence", "kinFit_result_probability", "kinFit_chi2"]
    # vars_to_plot += ["SVfit_valid", "SVfit_pt", "SVfit_eta", "SVfit_phi", "SVfit_m", "SVfit_mt", "MT2"]
    # vars_to_plot += ["SVfit_pt_error",  "SVfit_eta_error",  "SVfit_phi_error", "SVfit_m_error", "SVfit_mt_error"]

    #### Taggers #####
    #vars_to_plot += ["tau1_idDeepTau2017v2p1VSe", "tau1_idDeepTau2017v2p1VSmu", "tau1_idDeepTau2017v2p1VSjet", "tau1_idDeepTau2018v2p5VSe", "tau1_idDeepTau2018v2p5VSmu", "tau1_idDeepTau2018v2p5VSjet"]
    #vars_to_plot += ["tau2_idDeepTau2017v2p1VSe", "tau2_idDeepTau2017v2p1VSmu", "tau2_idDeepTau2017v2p1VSjet", "tau2_idDeepTau2018v2p5VSe", "tau2_idDeepTau2018v2p5VSmu", "tau2_idDeepTau2018v2p5VSjet"]
    #vars_to_plot += ["b1_btagDeepFlavB", "b1_btagDeepFlavCvB", "b1_btagDeepFlavCvL", "b1_particleNetAK4_B", "b1_particleNetAK4_CvsB", "b1_particleNetAK4_CvsL", "b1_HHbtag"]#, "b1_hadronFlavour"]
    #vars_to_plot += ["b2_btagDeepFlavB", "b2_btagDeepFlavCvB", "b2_btagDeepFlavCvL", "b2_particleNetAK4_B", "b2_particleNetAK4_CvsB", "b2_particleNetAK4_CvsL", "b2_HHbtag"]#, "b2_hadronFlavour"]

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

def getYear(period):
    year_dict = {
        'Run2_2016_HIPM':'2016_HIPM',
        'Run2_2016':'2016',
        'Run2_2017':'2017',
        'Run2_2018':'2018',
    }
    return year_dict[period]


def GetSamples(samples, backgrounds, wantExt=True, signals=['GluGluToRadion','GluGluToBulkGraviton']):
    global samples_to_consider
    samples_to_consider = ['data']

    for sample_name in samples.keys():
        if sample_name == 'GLOBAL': continue
        sample_type = samples[sample_name]['sampleType']
        if sample_type in signals:
            samples_to_consider.append(sample_name)
            continue
        if sample_type not in signals and sample_name not in backgrounds.keys(): continue
        #print(sample_name)
        if 'hasExt1' in backgrounds[sample_name].keys() and backgrounds[sample_name]['hasExt1']==True and wantExt == True:
            samples_to_consider.append("{sample_name}_ext1")
        samples_to_consider.append(sample_name)
    return samples_to_consider

backgrounds = None
def load_background_samples():
    global backgrounds
    background_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config','background_samples.yaml')
    with open(background_config,'r') as f:
        backgrounds= yaml.safe_load(f)
    return backgrounds

# **************** 1D HISTOGRAMS *******************
class HistProducerFileTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 5.0)
    n_cpus = copy_param(HTCondorWorkflow.n_cpus, 1)
    hist_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'plot','histograms.yaml')
    hists = load_hist_config(hist_config)
    vars_to_plot = load_vars_to_plot(hists)
    backgrounds = load_background_samples()


    def workflow_requires(self):
        branch_map = self.create_branch_map()
        workflow_dict = {}
        workflow_dict["anaTuple"] = {
            idx: AnaTupleTask.req(self, branch=br, branches=())
            for idx, (sample,sample_name_ext1, br,var, n_branch) in branch_map.items() if sample !='data' #and var==vars_to_plot[0]
        }
        workflow_dict["dataMergeTuple"] = {
            idx: DataMergeTask.req(self, branch=0, branches=())
            for idx, (sample, sample_name_ext1,br,var, n_branch) in branch_map.items() if sample =='data'
        }
        return workflow_dict


    def requires(self):
        sample_name,sample_name_ext1, prod_br,var,n_branch = self.branch_data
        if sample_name =='data':
            return [ DataMergeTask.req(self,max_runtime=DataMergeTask.max_runtime._default, branch=prod_br, branches=())]
        return [ AnaTupleTask.req(self, branch=prod_br, max_runtime=AnaTupleTask.max_runtime._default, branches=())]

    def create_branch_map(self):
        n = 0
        branches = {}
        anaProd_branch_map = AnaTupleTask.req(self, branch=-1, branches=()).create_branch_map()
        samples_to_consider = GetSamples(self.samples, backgrounds, True)
        for var in vars_to_plot:
            for prod_br,(sample_id, sample_name, sample_type, input_file) in anaProd_branch_map.items():
                sample_name_split = sample_name.split('_')
                sample_name_ext1 = sample_name
                if sample_name_split[-1]=='ext1':
                    sample_name = '_'.join(sample_name_split[:-1])
                if sample_name not in samples_to_consider or sample_name=='data': continue
                #print(sample_name)
                branches[n] = (sample_name,sample_name_ext1, prod_br,var, n)
                n+=1
            branches[n] = ('data','data', 0, var, n)
            n+=1
        return branches

    def output(self):
        if len(self.branch_data) == 0:
            return self.local_target('dummy.txt')

        sample_name,sample_name_ext1, prod_br,var,n_branch = self.branch_data
        input_file = self.input()[0].path
        outFileName = os.path.basename(self.input()[0].path)
        #print(outFileName)
        if sample_name_ext1 != sample_name:
            outFileName_split = outFileName.split('.')
            outFileName = f'{outFileName_split[0]}_{prod_br}.root'
        output_path = os.path.join(self.period, sample_name, self.version,'tmp', var, outFileName)
        return self.remote_target(output_path,  fs=self.fs_histograms)

    def run(self):
        sample_name,sample_name_ext1, prod_br,var,n_branch = self.branch_data
        if len(self.input()) > 1:
            raise RuntimeError(f"multple input files!! {' '.join(f.path for f in self.input())}")
        input_file = self.input()[0]
        inputFileName = os.path.basename(self.input()[0].path)
        print(f'input file is {input_file.path}')
        hist_config = self.hist_config
        global_config = os.path.join(self.ana_path(), 'config','HH_bbtautau', f'global.yaml')
        unc_config = os.path.join(self.ana_path(), 'config',self.period, f'weights.yaml')
        unc_cfg_dict = load_unc_config(unc_config)
        sample_config = os.path.join(self.ana_path(), 'config',self.period, f'samples.yaml')
        anaCache_path = os.path.join('anaCache',self.period, sample_name_ext1, self.version, inputFileName)
        HistProducerFile = os.path.join(self.ana_path(), 'Analysis', 'HistProducerFile.py')
        print(f'output file is {self.output().path}')
        with input_file.localize("r") as local_input, self.remote_target(anaCache_path, fs=self.fs_anaCache).localize("r") as local_anacache, self.output().localize("w") as local_output:
            #print(outFile.path)
            HistProducerFile_cmd = ['python3', HistProducerFile,'--inFile', local_input.path, '--cacheFile', local_anacache.path, '--outFileName',local_output.path, '--dataset', sample_name, '--compute_unc_variations', 'True', '--compute_rel_weights', 'True', '--uncConfig', unc_config, '--histConfig', hist_config,'--globalConfig', global_config, '--sampleConfig', sample_config, '--var', var]
            if 'deepTau2p5' in self.version.split('_'):
                print("deepTau2p5 in use")
                HistProducerFile_cmd.extend([ '--deepTauVersion', 'v2p5'])
            ps_call(HistProducerFile_cmd,verbose=1)


class HistProducerSampleTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 1.0)
    n_cpus = copy_param(HTCondorWorkflow.n_cpus, 1)
    hist_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'plot','histograms.yaml')
    hists = load_hist_config(hist_config)
    vars_to_plot = load_vars_to_plot(hists)
    backgrounds = load_background_samples()


    def workflow_requires(self):
        self_map = self.create_branch_map()
        all_samples = {}
        branches = {}
        for br_idx, (smpl_name, idx_list, var) in self_map.items():
            if var not in all_samples:
                all_samples[var] = []
            for idx_list_idx in idx_list:
                all_samples[var].append(idx_list_idx)
        workflow_dict = {}
        n=0
        for var in all_samples.keys():
            workflow_dict[var] = {
                n: HistProducerFileTask.req(self, branches=tuple((idx,) for idx in all_samples[var]))
            }
            n+=1
        return workflow_dict

    def requires(self):
        sample_name, idx_list,var  = self.branch_data
        deps = [HistProducerFileTask.req(self, max_runtime=HistProducerFileTask.max_runtime._default,branch=idx) for idx in idx_list ]
        return deps

    def create_branch_map(self):
        branches = {}
        histProducerFile_map = HistProducerFileTask.req(self,branch=-1, branches=()).create_branch_map()
        all_samples = {}
        samples_to_consider = GetSamples(self.samples, backgrounds, True)
        for j,(sample_name,sample_name_ext1,prod_br,var,n_branch)  in histProducerFile_map.items():
            if sample_name not in all_samples:
                all_samples[sample_name] = {}
            if var not in all_samples[sample_name]:
                all_samples[sample_name][var]=[]
            all_samples[sample_name][var].append(n_branch)
        k=0
        for n, key in enumerate(all_samples.items()):
            (smpl_name, dictionary)=key
            if smpl_name not in samples_to_consider: continue
            for m, (var, idx_list) in enumerate(dictionary.items()):
                branches[k] = (smpl_name, idx_list, var)
                k+=1
        return branches

    def output(self):
        sample_name, idx_list,var  = self.branch_data
        local_files_target = []
        outFileName_histProdSample = f'{var}.root'
        output_path = os.path.join(self.period, sample_name, self.version, var, outFileName_histProdSample)
        return self.remote_target(output_path,  fs=self.fs_histograms)

    def run(self):
        sample_name, idx_list,var  = self.branch_data
        HistProducerSample = os.path.join(self.ana_path(), 'Analysis', 'HistProducerSample.py')
        with contextlib.ExitStack() as stack:
            #print(self.input())
            local_inputs = [stack.enter_context(inp.localize('r')).path for inp in self.input()]
            with self.output().localize("w") as tmp_local_file:
                HistProducerSample_cmd = ['python3', HistProducerSample,'--outFile', tmp_local_file.path]#, '--remove-files', 'True']
                HistProducerSample_cmd.extend(local_inputs)
                ps_call(HistProducerSample_cmd,verbose=1)


class MergeTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 30.0)
    hist_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'plot','histograms.yaml')
    hists = load_hist_config(hist_config)
    backgrounds = load_background_samples()


    def workflow_requires(self):
        histProducerSample_map = HistProducerSampleTask.req(self,branch=-1, branches=()).create_branch_map()
        all_samples = {}
        branches = {}
        for br_idx, (smpl_name, idx_list, var) in histProducerSample_map.items():
            if var not in all_samples:
                all_samples[var] = []
            all_samples[var].append(br_idx)
        workflow_dict = {}
        n=0
        for var in all_samples.keys():
            workflow_dict[var] = {
                n: HistProducerSampleTask.req(self, branches=tuple((idx,) for idx in all_samples[var]))
            }
            n+=1
        return workflow_dict

    def requires(self):
        uncName, var, branches_idx = self.branch_data
        deps = [HistProducerSampleTask.req(self, max_runtime=HistProducerSampleTask.max_runtime._default, branch=prod_br) for prod_br in branches_idx ]
        return deps

    def create_branch_map(self):
        histProducerSample_map = HistProducerSampleTask.req(self,branch=-1, branches=()).create_branch_map()
        all_samples = {}
        branches = {}
        for br_idx, (smpl_name, idx_list, var) in histProducerSample_map.items():
            if var not in all_samples:
                all_samples[var] = []
            all_samples[var].append(br_idx)
        k=0
        uncNames = ['Central']
        unc_config = os.path.join(self.ana_path(), 'config',self.period, f'weights.yaml')
        unc_cfg_dict = load_unc_config(unc_config)
        uncNames.extend(list(unc_cfg_dict['norm'].keys()))
        uncNames.extend([unc for unc in unc_cfg_dict['shape']])
        for uncName in uncNames:
            if uncName in uncs_to_exclude[self.period]: continue
            for n, key in enumerate(all_samples.items()):
                var, branches_idx = key
                branches[k] = (uncName, var, branches_idx)
                k+=1
        return branches


    def output(self):
        if len(self.branch_data) == 0:
            return self.local_target('dummy.txt')
        uncName, var, branches_idx = self.branch_data
        outDir_MergeHists = os.path.join(self.period, 'all_histograms', self.version, var)
        outFileName_MergeHists = f'all_histograms_{var}_{uncName}.root'
        output_path = os.path.join(outDir_MergeHists, outFileName_MergeHists)
        return self.remote_target(output_path,  fs=self.fs_histograms)

    def run(self):
        uncName, var, branches_idx = self.branch_data
        sample_config = os.path.join(self.ana_path(), 'config',self.period, f'samples.yaml')
        unc_config = os.path.join(self.ana_path(), 'config',self.period, f'weights.yaml')
        global_config = os.path.join(self.ana_path(), 'config','HH_bbtautau', f'global.yaml')
        jsonFile_name = os.path.join('jsonFiles', self.period, self.version, 'all_ratios',f'all_ratios_{var}_{uncName}.yaml')
        MergerProducer = os.path.join(self.ana_path(), 'Analysis', 'HistMerger.py')
        all_inputs = []
        samples_to_consider = GetSamples(self.samples, backgrounds, True)
        for sample_name in self.samples.keys():
            if sample_name not in samples_to_consider: continue
            output_path_hist_prod_sample = os.path.join(self.period, sample_name, self.version, var, f'{var}.root')
            all_inputs.append((self.remote_target(output_path_hist_prod_sample, fs=self.fs_histograms),sample_name))
        output_path_hist_prod_sample_data = os.path.join(self.period, 'data', self.version, var, f'{var}.root')
        all_inputs.append((self.remote_target(output_path_hist_prod_sample_data, fs=self.fs_histograms),'data'))
        all_datasets=[]
        with contextlib.ExitStack() as stack:
            local_inputs = []
            for inp, smpl in all_inputs:
                local_inputs.append(stack.enter_context(inp.localize('r')).path)
                all_datasets.append(smpl)
            with self.output().localize("w") as tmp_local_file, self.remote_target(jsonFile_name, fs=self.fs_json).localize("w") as json_file :
                tmpFile = tmp_local_file.path
                jsonFile = json_file.path
                dataset_names = ','.join(smpl for smpl in all_datasets)
                MergerProducer_cmd = ['python3', MergerProducer,'--outFile', tmpFile, '--jsonFile', jsonFile, '--var', var, '--uncSource', uncName, '--uncConfig', unc_config, '--sampleConfig', sample_config, '--datasetFile', dataset_names,  '--year', getYear(self.period) , '--globalConfig', global_config]#, '--remove-files', 'True']
                MergerProducer_cmd.extend(local_inputs)
                ps_call(MergerProducer_cmd,verbose=1)


class HaddMergedTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 1.0)
    hist_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'plot','histograms.yaml')
    hists = load_hist_config(hist_config)
    vars_to_plot=load_vars_to_plot(hists)
    backgrounds = load_background_samples()

    def workflow_requires(self):
        return { "Merge" : MergeTask.req(self, branches=()) }

    def requires(self):
        return [MergeTask.req(self, branches=())]

    def create_branch_map(self):
        n = 0
        branches = {}
        for var in vars_to_plot :
            branches[n] = var
            n += 1
        return branches

    def output(self):
        var = self.branch_data
        output_path = os.path.join(self.period, 'all_histograms', self.version, var, f'all_histograms_{var}_HAdded.root')
        return self.remote_target(output_path, fs=self.fs_histograms)


    def run(self):
        var = self.branch_data
        sample_config = os.path.join(self.ana_path(), 'config',self.period, f'samples.yaml')
        unc_config = os.path.join(self.ana_path(), 'config',self.period, f'weights.yaml')
        unc_cfg_dict = load_unc_config(unc_config)
        inDir_allHistograms = os.path.join(self.period, 'all_histograms', self.version, var)
        uncNames = []
        uncNames.extend(list(unc_cfg_dict['norm'].keys()))
        uncNames.extend([unc for unc in unc_cfg_dict['shape']])
        inFileNameCentral = f'all_histograms_{var}_Central.root'
        inFileCentral = os.path.join(inDir_allHistograms,inFileNameCentral)
        all_inputs = []
        all_uncertainties = ['Central']
        for uncName in ['Central']+uncNames :
            if uncName in uncs_to_exclude[self.period]: continue
            inFileName =  f'all_histograms_{var}_{uncName}.root'
            #print(inFileName)
            all_inputs.append(self.remote_target(os.path.join(inDir_allHistograms,inFileName), fs=self.fs_histograms))
            all_uncertainties.append(uncName)
        tmp_outFileName =  f'all_histograms_{var}_tmp.root'

        tmp_outFile = self.remote_target(os.path.join(inDir_allHistograms,tmp_outFileName), fs=self.fs_histograms)
        all_uncertainties_string = ','.join(unc for unc in all_uncertainties)
        HaddMergedHistsProducer = os.path.join(self.ana_path(), 'Analysis', 'hadd_merged_hists.py')
        RenameHistsProducer = os.path.join(self.ana_path(), 'Analysis', 'renameHists.py')
        ShapeOrLogNormalProducer = os.path.join(self.ana_path(), 'Analysis', 'ShapeOrLogNormal.py')
        jsonFile_name = os.path.join('jsonFiles', self.period, self.version, 'fitResults',f'slopeInfo_{var}.yaml')
        with contextlib.ExitStack() as stack:
            local_inputs = []
            for inp in all_inputs:
                local_inputs.append(stack.enter_context(inp.localize('r')).path)
            with tmp_outFile.localize("w") as tmpFile:
                HaddMergedHistsProducer_cmd = ['python3', HaddMergedHistsProducer,'--outFile', tmpFile.path, '--var', var]
                HaddMergedHistsProducer_cmd.extend(local_inputs)
                ps_call(HaddMergedHistsProducer_cmd,verbose=1)

        with tmp_outFile.localize("r") as tmpFile, self.output().localize("w") as outFile:
            RenameHistsProducer_cmd = ['python3', RenameHistsProducer,'--inFile', tmpFile.path, '--outFile', outFile.path, '--var', var, '--year', getYear(self.period)]
            ps_call(RenameHistsProducer_cmd,verbose=1)
