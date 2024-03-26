import law
import luigi
import os
import shutil
import time
import threading
import yaml

from RunKit.sh_tools import sh_call
from RunKit.checkRootFile import checkRootFileSafe
from RunKit.crabLaw import cond as kInit_cond,update_kinit_thread
from run_tools.law_customizations import Task, HTCondorWorkflow, copy_param, get_param_value


unc_cfg_dict = None
def load_unc_config(unc_cfg):
    global unc_cfg_dict
    with open(unc_cfg, 'r') as f:
        unc_cfg_dict = yaml.safe_load(f)
    return unc_cfg_dict

class AnaCacheTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 10.0)

    def create_branch_map(self):
        n = 0
        branches = {}
        for sample_name in sorted(self.samples.keys()):
            #print(sample_name)
            isData = self.samples[sample_name]['sampleType'] == 'data'
            branches[n] = (sample_name, isData)
            n += 1
        return branches

    def output(self):
        sample_name, isData = self.branch_data
        outDir = os.path.join(self.central_anaCache_path(), sample_name)
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        sample_out = os.path.join(outDir, 'anaCache.yaml')
        return law.LocalFileTarget(sample_out)

    def run(self):
        sample_name, isData = self.branch_data
        print(f'Creating anaCache for sample {sample_name} into {self.output().path}')
        if isData:
            self.output().touch()
        else:
            producer = os.path.join(self.ana_path(), 'AnaProd', 'anaCacheProducer.py')
            inDir = os.path.join(self.central_nanoAOD_path(), sample_name)
            #if self.period!='Run2_2018':
            if self.version.split('_')[-1]=='HTT':
                inDir = os.path.join(self.central_nanoAOD_path_HLepRare(), sample_name)
            os.makedirs(os.path.dirname(self.output().path), exist_ok=True)
            sh_call(['python3', producer, '--config', self.sample_config, '--inDir', inDir, '--sample', sample_name,
                    '--outFile', self.output().path, '--customisations', self.customisations ], env=self.cmssw_env())
        print(f'anaCache for sample {sample_name} is created in {self.output().path}')



class InputFileTask(Task, law.LocalWorkflow):
    max_files_per_job = luigi.IntParameter(default=1, description="maximum number of input files per job")

    def create_branch_map(self):
        branches = {}
        for n, sample_name in enumerate(sorted(self.samples.keys())):
            branches[n] = sample_name
        return branches

    def output(self):
        sample_name = self.branch_data
        sample_out = os.path.join(self.local_path(sample_name, "input_files.txt"))
        return law.LocalFileTarget(sample_out)

    def run(self):
        sample_name = self.branch_data
        print(f'Creating inputFile for sample {sample_name} into {self.output().path}')
        os.makedirs(os.path.join(self.local_path(),sample_name), exist_ok=True)
        txtFile_tmp = os.path.join(self.local_path(), sample_name, "tmp.txt")
        inDir = os.path.join(self.central_nanoAOD_path(), sample_name)
        #if self.period!='Run2_2018':
        if self.version.split('_')[-1]=='HTT':
            inDir = os.path.join(self.central_nanoAOD_path_HLepRare(), sample_name)
        #print(f"inDir is {inDir}")
        input_files = []
        for root, dirs, files in os.walk(inDir):
            for file in files:
                if file.endswith('.root') and not file.startswith('.'):
                    if os.path.join(root, file) not in input_files:
                        input_files.append(os.path.join(root, file))
        #print(f"inDir is {input_files}")
        with open(txtFile_tmp, 'w') as inputFileTxt:
            for input_line in input_files:
                inputFileTxt.write(input_line+'\n')
        finalFile = self.output().path
        shutil.move(txtFile_tmp,finalFile)
        print(f'inputFile for sample {sample_name} is created in {self.output().path}')



class AnaTupleTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 30.0)

    def workflow_requires(self):
        return { "anaCache" : AnaCacheTask.req(self, branches=()), "inputFile": InputFileTask.req(self,workflow='local', branches=()) }

    def requires(self):
        sample_id, sample_name, sample_type, input_file = self.branch_data
        return [ AnaCacheTask.req(self, branch=sample_id, max_runtime=AnaCacheTask.max_runtime._default, branches=()), InputFileTask.req(self, branch=sample_id, workflow='local', branches=()) ]

    def create_branch_map(self):
        n = 0
        branches = {}
        for sample_id, sample_name in enumerate(sorted(self.samples.keys())):
            inputFileTxt = InputFileTask.req(self, branch=sample_id,workflow='local', branches=(sample_id,)).output().path
            with open(inputFileTxt, 'r') as inputtxtFile:
                input_files = inputtxtFile.read().splitlines()
            #print(input_files)
            if len(input_files) == 0:
                raise RuntimeError(f"AnaTupleTask: no input files found for {sample_name}")
            for input_file in input_files:
                branches[n] = (sample_id, sample_name, self.samples[sample_name]['sampleType'], input_file)
                n += 1
        return branches

    def output(self, force_pre_output=False):
        sample_id, sample_name, sample_type, input_file = self.branch_data
        outFileName = os.path.basename(input_file)
        #outDir_tmp = os.path.join('/eos/user/a/androsov/valeria/','anaTuples', self.period, self.version)
        outDir_tmp = self.central_anaTuples_path()
        if not os.path.exists(outDir_tmp):
            os.makedirs(outDir_tmp)
        out = os.path.join(outDir_tmp, sample_name,outFileName)
        return law.LocalFileTarget(out)

    def run(self):
        thread = threading.Thread(target=update_kinit_thread)
        thread.start()
        try:
            job_home, remove_job_home = self.law_job_home()
            sample_id, sample_name, sample_type, input_file = self.branch_data
            if self.test: print(f"sample_id= {sample_id}\nsample_name = {sample_name}\nsample_type = {sample_type}\ninput_file = {input_file}")
            producer_anatuples = os.path.join(self.ana_path(), 'AnaProd', 'anaTupleProducer.py')
            anaCache = os.path.join(self.central_anaCache_path(), sample_name, 'anaCache.yaml')
            outdir_anatuples = os.path.join(job_home, 'anaTuples', sample_name)
            anatuple_cmd = [ 'python3', producer_anatuples, '--config', self.sample_config, '--inFile', input_file,
                        '--outDir', outdir_anatuples, '--sample', sample_name, '--anaCache', anaCache, '--customisations',
                        self.customisations, '--compute_unc_variations', 'True', '--store-noncentral']
            if self.test:
                anatuple_cmd.extend(['--nEvents', '100'])
            sh_call(anatuple_cmd, env=self.cmssw_env(),verbose=1)
            producer_skimtuples = os.path.join(self.ana_path(), 'AnaProd', 'SkimProducer.py')
            outdir_skimtuples = os.path.join(job_home, 'skim', sample_name)
            outFileName = os.path.basename(input_file)
            if self.test: print(f"outFileName is {outFileName}")
            if sample_type!='data':
                skimtuple_cmd = ['python3', producer_skimtuples, '--inputDir',outdir_anatuples, '--centralFile',outFileName, '--workingDir', outdir_skimtuples,
                        '--outputFile', outFileName]
                if self.test:
                    skimtuple_cmd.extend(['--test' , 'True'])
                sh_call(skimtuple_cmd,verbose=1)
            tmpFile = os.path.join(outdir_skimtuples, outFileName)
            if sample_type=='data':
                tmpFile = os.path.join(outdir_anatuples, outFileName)

            outDir_tmp = self.central_anaTuples_path()
            outdir_final = os.path.join(outDir_tmp, sample_name)
            os.makedirs(outdir_final, exist_ok=True)
            finalFile = self.output().path
            if self.test: print(f"finalFile is {finalFile}")
            shutil.copy(tmpFile, finalFile)
            if os.path.exists(finalFile):
                os.remove(tmpFile)
            else:
                raise RuntimeError(f"file {tmpFile} was not copied in {finalFile} ")
            if remove_job_home:
                shutil.rmtree(job_home)
            print(f'anaTuple for sample {sample_name} is created in {outdir_final}')
        finally:
            kInit_cond.acquire()
            kInit_cond.notify_all()
            kInit_cond.release()
            thread.join()


    @staticmethod
    def getOutputDir(central_anaTuples_path, sample_name):
        return os.path.join(central_anaTuples_path, sample_name)

class DataMergeTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 5.0)

    def workflow_requires(self):
        prod_branches = self.create_branch_map()
        workflow_dict = {}
        workflow_dict["anaTuple"] = {
            idx: AnaTupleTask.req(self, branches=tuple((br,) for br in branches))
            for idx, branches in prod_branches.items()
        }
        return workflow_dict

    def requires(self):
        prod_branches = self.branch_data
        deps = [AnaTupleTask.req(self, max_runtime=AnaCacheTask.max_runtime._default, branch=prod_br) for prod_br in prod_branches ]
        return deps


    def create_branch_map(self):
        anaProd_branch_map = AnaTupleTask.req(self, branch=-1, branches=()).create_branch_map()
        prod_branches = []
        for prod_br, (sample_id, sample_name, sample_type, input_file) in anaProd_branch_map.items():
            if sample_type != "data": continue
            prod_branches.append(prod_br)
        return { 0: prod_branches }

    def output(self, force_pre_output=False):
        outDir_tmp = self.central_anaTuples_path()
        out = os.path.join(outDir_tmp, 'data','nano.root')
        return law.LocalFileTarget(out)

    def run(self):
        producer_dataMerge = os.path.join(self.ana_path(), 'AnaProd', 'MergeNtuples.py')
        outDir_tmp = self.central_anaTuples_path()
        tmpFile = os.path.join(outDir_tmp, 'data', 'data_tmp.root')
        dataMerge_cmd = [ 'python3', producer_dataMerge, '--outFile', tmpFile, '--useUproot', 'True']
        dataMerge_cmd.extend([f.path for f in self.input()])
        sh_call(dataMerge_cmd,verbose=1)
        finalFile = self.output().path
        if self.test: print(f"finalFile is {finalFile}")
        print(f"finalFile is {finalFile}")
        shutil.copy(tmpFile, finalFile)
        if os.path.exists(finalFile):
            os.remove(tmpFile)



class AnaCacheTupleTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 30.0)
    n_cpus = copy_param(HTCondorWorkflow.n_cpus, 1)

    def workflow_requires(self):
        branch_map = self.create_branch_map()
        workflow_dict = {}
        workflow_dict["anaTuple"] = {
            idx: AnaTupleTask.req(self, branch=br, branches=())
            for idx, (sample, sample_type,br) in branch_map.items()
        }
        return workflow_dict

    def requires(self):
        sample_name,sample_type, prod_br = self.branch_data
        return [ AnaTupleTask.req(self, branch=prod_br, max_runtime=AnaTupleTask.max_runtime._default, branches=())]

    def create_branch_map(self):
        n = 0
        branches = {}
        anaProd_branch_map = AnaTupleTask.req(self, branch=-1, branches=()).create_branch_map()
        sample_id_data = 0
        for prod_br, (sample_id, sample_name, sample_type, input_file) in anaProd_branch_map.items():
            if sample_type =='QCD':
                continue
            branches[n] = (sample_name, sample_type,prod_br)
            n+=1
        #branches[n+1] = ('data', 0)
        return branches

    def output(self):
        sample_name, sample_type,prod_br = self.branch_data
        input_file = self.input()[0].path
        #print(f"inputFile is {input_file}")
        fileName  = os.path.basename(input_file)
        #print(f"filename is {fileName}")
        outDir = os.path.join(self.central_anaCache_path(), sample_name, self.version)
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        outFileName = os.path.join(outDir, fileName)
        return law.LocalFileTarget(outFileName)

    def run(self):
        thread = threading.Thread(target=update_kinit_thread)
        thread.start()
        finalFile = self.output().path
        try:
            job_home, remove_job_home = self.law_job_home()
            sample_name, sample_type,prod_br = self.branch_data
            if len(self.input()) > 1:
                raise RuntimeError(f"multple input files!! {' '.join(f.path for f in self.input())}")
            input_file = self.input()[0].path
            #print(f"inputFile is {input_file}")
            fileName_list = os.path.basename(input_file).split('.')
            #print(fileName_list)
            fileName = fileName_list[0]
            #print(fileName)
            #print(f"filename is {fileName}")
            sample_config = self.sample_config
            unc_config = os.path.join(os.getenv("ANALYSIS_PATH"), 'config', 'weight_definition.yaml')
            unc_cfg_dict = load_unc_config(unc_config)
            #print(f"\nsample_name = {sample_name}\ninput_file = {input_file}")
            outdir_anacacheTuples = os.path.join(job_home, 'anaCacheTuples', sample_name)
            #print(outdir_anacacheTuples)
            if os.path.isdir(outdir_anacacheTuples):
                shutil.rmtree(outdir_anacacheTuples)
            os.makedirs(outdir_anacacheTuples, exist_ok=True)
            outFileName_anacacheTuples = os.path.join(outdir_anacacheTuples, fileName)
            print(outFileName_anacacheTuples)
            anaCacheTupleProducer = os.path.join(self.ana_path(), 'AnaProd', 'anaCacheTupleProducer.py')
            anaCacheTupleProducer_cmd = ['python3', anaCacheTupleProducer,'--inFileName', input_file, '--outFileName',outFileName_anacacheTuples,  '--uncConfig', unc_config]
            if sample_name !='data':
                anaCacheTupleProducer_cmd.extend(['--compute_unc_variations', 'True'])
            if self.version == 'v9_deepTau2p5':
                anaCacheTupleProducer_cmd.extend([ '--deepTauVersion', 'v2p5'])
            sh_call(anaCacheTupleProducer_cmd, env=self.cmssw_env(),verbose=1)
            print(f"finished to produce anacachetuple")
            outdir_final = os.path.join(self.central_anaCache_path(), sample_name, self.version)
            if not os.path.exists(outdir_final):
                os.makedirs(outdir_final)
            os.makedirs(outdir_final, exist_ok=True)
            finalFile = self.output().path
            print(f"finalFile is {finalFile}")
            outFileName_anacacheTuples += '.root'
            shutil.copy(outFileName_anacacheTuples, finalFile)
            if os.path.exists(finalFile):
                os.remove(outFileName_anacacheTuples)
            else:
                raise RuntimeError(f"file {tmpFile} was not copied in {finalFile} ")
            if remove_job_home:
                shutil.rmtree(job_home)
            print(f'anaCache for sample {sample_name} is created in {outdir_final} in {finalFile}')

            '''
            producer_skimtuples = os.path.join(self.ana_path(), 'AnaProd', 'SkimProducer.py')
            outdir_skimcachetuples = os.path.join(job_home, 'skimcache', sample_name)

            if os.path.isdir(outdir_skimcachetuples):
                shutil.rmtree(outdir_skimcachetuples)
            os.makedirs(outdir_skimcachetuples, exist_ok=True)
            outFileName = os.path.basename(input_file)
            if self.test: print(f"outFileName is {outFileName}")
            if sample_name!='data':
                skimtuple_cmd = ['python3', producer_skimtuples, '--inputDir',outdir_anacacheTuples, '--centralFile',outFileName, '--workingDir', outdir_skimcachetuples,
                        '--outputFile', outFileName]
                if self.test:
                    skimtuple_cmd.extend(['--test' , 'True'])
                #sh_call(skimtuple_cmd,verbose=1)
                print(skimtuple_cmd)
            tmpFile = os.path.join(outdir_skimcachetuples, outFileName)
            if sample_name=='data':
                tmpFile = os.path.join(outdir_anacacheTuples, outFileName)
            outdir_final = os.path.join(self.central_anaCache_path(), sample_name, self.version)
            if not os.path.exists(outdir_final):
                os.makedirs(outdir_final)
            os.makedirs(outdir_final, exist_ok=True)
            finalFile = self.output().path
            if self.test: print(f"finalFile is {finalFile}")
            shutil.copy(tmpFile, finalFile)
            if os.path.exists(finalFile):
                os.remove(tmpFile)
            else:
                raise RuntimeError(f"file {tmpFile} was not copied in {finalFile} ")
            if remove_job_home:
                shutil.rmtree(job_home)
            print(f'anaTuple for sample {sample_name} is created in {outdir_final}')
            '''
        finally:
            kInit_cond.acquire()
            kInit_cond.notify_all()
            kInit_cond.release()
            thread.join()


class DataCacheMergeTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 5.0)

    def workflow_requires(self):
        prod_branches = self.create_branch_map()
        workflow_dict = {}
        workflow_dict["anaCacheTuple"] = {
            idx: AnaCacheTupleTask.req(self, branches=tuple((br,) for br in branches))
            for idx, branches in prod_branches.items()
        }
        return workflow_dict

    def requires(self):
        prod_branches = self.branch_data
        deps = [AnaCacheTupleTask.req(self, max_runtime=AnaCacheTask.max_runtime._default, branch=prod_br) for prod_br in prod_branches ]
        return deps


    def create_branch_map(self):
        anaProd_branch_map = AnaCacheTupleTask.req(self, branch=-1, branches=()).create_branch_map()
        prod_branches = []
        for prod_br, (sample_name, sample_type,branch) in anaProd_branch_map.items():
            if sample_type != "data": continue
            prod_branches.append(prod_br)
        return { 0: prod_branches }

    def output(self, force_pre_output=False):
        outDir = os.path.join(self.central_anaCache_path(), 'data', self.version)
        out = os.path.join(outDir, 'nano.root')
        return law.LocalFileTarget(out)

    def run(self):
        producer_dataMerge = os.path.join(self.ana_path(), 'AnaProd', 'MergeNtuples.py')
        outDir = os.path.join(self.central_anaCache_path(), 'data', self.version)
        tmpFile = os.path.join(outDir, 'data_tmp.root')
        dataMerge_cmd = [ 'python3', producer_dataMerge, '--outFile', tmpFile ]
        dataMerge_cmd.extend([f.path for f in self.input()])
        sh_call(dataMerge_cmd,verbose=1)
        finalFile = self.output().path
        if self.test: print(f"finalFile is {finalFile}")
        print(f"finalFile is {finalFile}")
        shutil.copy(tmpFile, finalFile)
        if os.path.exists(finalFile):
            os.remove(tmpFile)

