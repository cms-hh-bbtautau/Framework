import copy
import contextlib
import json
import law
import os
import shutil
import threading
import yaml


from RunKit.run_tools import ps_call, natural_sort
from RunKit.crabLaw import cond as kInit_cond, update_kinit_thread
from run_tools.law_customizations import Task, HTCondorWorkflow, copy_param
from Common.Utilities import SerializeObjectToString
from AnaProd.anaCacheProducer import addAnaCaches

class InputFileTask(Task, law.LocalWorkflow):
    def __init__(self, *args, **kwargs):
        kwargs['workflow'] = 'local'
        super(InputFileTask, self).__init__(*args, **kwargs)

    def create_branch_map(self):
        branches = {}
        for sample_id, sample_name in self.iter_samples():
            branches[sample_id] = sample_name
        return branches

    def output(self):
        sample_name = self.branch_data
        return self.local_target('input_files', f'{sample_name}.txt')

    def run(self):
        sample_name = self.branch_data
        print(f'Creating inputFile for sample {sample_name} into {self.output().path}')
        with self.output().localize("w") as out_local_file:
            input_files = []
            for file in natural_sort(self.fs_nanoAOD.listdir(sample_name)):
                if file.endswith(".root"):
                    input_files.append(file)
            with open(out_local_file.path, 'w') as inputFileTxt:
                for input_line in input_files:
                    inputFileTxt.write(input_line+'\n')
        print(f'inputFile for sample {sample_name} is created in {self.output().path}')

    @staticmethod
    def load_input_files(input_file_list, sample_name, fs=None, return_uri=False):
        input_files = []
        with open(input_file_list, 'r') as txt_file:
            for file in txt_file.readlines():
                file_path = os.path.join(sample_name, file.strip())
                file_full_path = fs.uri(file_path) if return_uri else file_path
                input_files.append(file_full_path)
        if len(input_files) == 0:
            raise RuntimeError(f"No input files found for {sample_name}")
        return input_files

class AnaCacheTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 10.0)

    def create_branch_map(self):
        branches = {}
        for sample_id, sample_name in self.iter_samples():
            isData = self.samples[sample_name]['sampleType'] == 'data'
            branches[sample_id] = (sample_name, isData)
        return branches

    def requires(self):
        return [ InputFileTask.req(self) ]

    def workflow_requires(self):
        return { "inputFile": InputFileTask.req(self) }

    def output(self):
        sample_name, isData = self.branch_data
        return self.remote_target('anaCache', self.period, f'{sample_name}.yaml', fs=self.fs_anaCache)

    def run(self):
        sample_name, isData = self.branch_data
        if isData:
            self.output().touch()
            return
        print(f'Creating anaCache for sample {sample_name} into {self.output().uri()}')
        producer = os.path.join(self.ana_path(), 'AnaProd', 'anaCacheProducer.py')
        input_files = InputFileTask.load_input_files(self.input()[0].path, sample_name)
        ana_caches = []
        global_params_str = SerializeObjectToString(self.global_params)
        n_inputs = len(input_files)
        for input_idx, input_file in enumerate(input_files):
            input_target = self.remote_target(input_file, fs=self.fs_nanoAOD)
            print(f'[{input_idx+1}/{n_inputs}] {input_target.uri()}')
            with input_target.localize("r") as input_local:
                returncode, output, err = ps_call([ 'python3', producer, '--input-files', input_local.path,
                                                    '--global-params', global_params_str, '--verbose', '1' ],
                                                  env=self.cmssw_env, catch_stdout=True)
            ana_cache = json.loads(output)
            print(json.dumps(ana_cache))
            ana_caches.append(ana_cache)
        total_ana_cache = addAnaCaches(*ana_caches)
        print(f'total anaCache: {json.dumps(total_ana_cache)}')
        with self.output().localize("w") as output_local:
            with open(output_local.path, 'w') as file:
                yaml.dump(total_ana_cache, file)
        print(f'anaCache for sample {sample_name} is created in {self.output().uri()}')

class AnaTupleTask(Task, HTCondorWorkflow, law.LocalWorkflow):
    max_runtime = copy_param(HTCondorWorkflow.max_runtime, 30.0)

    def create_branch_map(self):
        input_file_task_complete = InputFileTask.req(self).complete()
        if not input_file_task_complete:
            self.cache_branch_map = False
            if not hasattr(self, '_branches_backup'):
                self._branches_backup = copy.deepcopy(self.branches)
            return { 0: () }
        self.cache_branch_map = True
        if hasattr(self, '_branches_backup'):
            self.branches = self._branches_backup
        branch_idx = 0
        branches = {}
        for sample_id, sample_name in self.iter_samples():
            input_file_list = InputFileTask.req(self, branch=sample_id, branches=(sample_id,)).output().path
            input_files = InputFileTask.load_input_files(input_file_list, sample_name)
            for input_file in input_files:
                branches[branch_idx] = (sample_id, sample_name, self.samples[sample_name]['sampleType'],
                                        self.remote_target(input_file, fs=self.fs_nanoAOD))
                branch_idx += 1
        return branches

    def workflow_requires(self):
        return { "anaCache" : AnaCacheTask.req(self),
                 "inputFile": InputFileTask.req(self) }

    def requires(self):
        sample_id, sample_name, sample_type, input_file = self.branch_data
        return [ AnaCacheTask.req(self, branch=sample_id, max_runtime=AnaCacheTask.max_runtime._default, branches=()) ]

    def output(self):
        if len(self.branch_data) == 0:
            return self.local_target('dummy.txt')
        sample_id, sample_name, sample_type, input_file = self.branch_data
        output_name = os.path.basename(input_file.path)
        output_path = os.path.join('anaTuples', self.version, self.period, sample_name, output_name)
        return self.remote_target(output_path, fs=self.fs_anaTuple)

    def run(self):
        sample_id, sample_name, sample_type, input_file = self.branch_data
        producer_anatuples = os.path.join(self.ana_path(), 'AnaProd', 'anaTupleProducer.py')
        producer_skimtuples = os.path.join(self.ana_path(), 'AnaProd', 'SkimProducer.py')
        thread = threading.Thread(target=update_kinit_thread)
        thread.start()
        anaCache_remote = self.input()[0]
        try:
            job_home, remove_job_home = self.law_job_home()
            print(f"sample_id = {sample_id}\nsample_name = {sample_name}\nsample_type = {sample_type}\n"
                  f"input_file = {input_file.uri()}")

            print("step 1: nanoAOD -> anaTupleS")
            outdir_anatuples = os.path.join(job_home, 'anaTuples', sample_name)
            anaTupleDef = os.path.join(self.ana_path(), self.global_params['anaTupleDef'])
            with input_file.localize("r") as local_input, anaCache_remote.localize("r") as anaCache_input:
                anatuple_cmd = [ 'python3', producer_anatuples, '--period', self.period,
                                 '--inFile', local_input.path, '--outDir', outdir_anatuples, '--sample', sample_name,
                                 '--anaTupleDef', anaTupleDef, '--anaCache', anaCache_input.path ]
                if len(self.customisations) > 0:
                    anatuple_cmd.extend([ '--customisations', self.customisations ])
                if self.global_params.get('compute_unc_variations', False):
                    anatuple_cmd.append('--compute-unc-variations')
                if self.global_params.get('store_noncentral', False):
                    anatuple_cmd.append('--store-noncentral')
                centralFileName = os.path.basename(local_input.path)
                if self.test:
                    anatuple_cmd.extend(['--nEvents', '100'])
                ps_call(anatuple_cmd, env=self.cmssw_env, verbose=1)

            print("step 2: anaTupleS -> skimTuplE")
            outdir_skimtuples = os.path.join(job_home, 'skim', sample_name)
            outFileName = os.path.basename(input_file.path)

            if self.test: print(f"outFileName is {outFileName}")
            tmpFile = os.path.join(outdir_skimtuples, outFileName)
            if sample_type!='data':
                skimtuple_cmd = [ 'python', producer_skimtuples, '--inputDir', outdir_anatuples,
                                  '--centralFile',centralFileName, '--workingDir', outdir_skimtuples,
                                  '--outputFile', outFileName ]
                if self.test:
                    skimtuple_cmd.extend(['--test' , 'True'])
                ps_call(skimtuple_cmd,verbose=1)
            else:
                tmpFile = os.path.join(outdir_anatuples, centralFileName)

            with self.output().localize("w") as tmp_local_file:
                out_local_path = tmp_local_file.path
                shutil.move(tmpFile, out_local_path)

            if remove_job_home:
                shutil.rmtree(job_home)
            print(f'anaTuple for sample {sample_name} is created in {self.output().uri()}')
        finally:
            kInit_cond.acquire()
            kInit_cond.notify_all()
            kInit_cond.release()
            thread.join()

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
        outFileName = 'nanoHTT_0.root'
        output_path = os.path.join('anaTuples', self.version, self.period, 'data', outFileName)
        return self.remote_target(output_path, fs=self.fs_anaTuple)

    def run(self):
        producer_dataMerge = os.path.join(self.ana_path(), 'AnaProd', 'MergeNtuples.py')
        with contextlib.ExitStack() as stack:
            local_inputs = [stack.enter_context(inp.localize('r')).path for inp in self.input()]
            with self.output().localize("w") as tmp_local_file:
                tmpFile = tmp_local_file.path
                dataMerge_cmd = [ 'python3', producer_dataMerge, '--outFile', tmpFile]#, '--useUproot', 'True']
                dataMerge_cmd.extend(local_inputs)
                ps_call(dataMerge_cmd,verbose=1)



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
        for prod_br,(sample_id, sample_name, sample_type, input_file) in anaProd_branch_map.items():
            if sample_type =='QCD':
                continue
            branches[n] = (sample_name, sample_type,prod_br)
            n+=1
        #branches[n+1] = ('data', 0)
        return branches


    def output(self):
        sample_name, sample_type,prod_br = self.branch_data
        outFileName = os.path.basename(self.input()[0].path)
        #print(outFileName)
        outDir = os.path.join('anaCache', self.period, sample_name, self.version)
        finalFile = os.path.join(outDir, outFileName)
        return self.remote_target(finalFile, fs=self.fs_anaCache)

    def run(self):
        sample_name, sample_type,prod_br = self.branch_data
        unc_config = os.path.join(self.ana_path(), 'config',self.period, f'weights.yaml')
        producer_anacachetuples = os.path.join(self.ana_path(), 'AnaProd', 'anaCacheTupleProducer.py')

        thread = threading.Thread(target=update_kinit_thread)
        thread.start()
        try:
            job_home, remove_job_home = self.law_job_home()
            input_file = self.input()[0]
            with input_file.localize("r") as local_input, self.output().localize("w") as outFile:
                anaCacheTupleProducer_cmd = ['python3', producer_anacachetuples,'--inFileName', local_input.path, '--outFileName', outFile.path,  '--uncConfig', unc_config]
                if sample_name !='data':
                    anaCacheTupleProducer_cmd.extend(['--compute_unc_variations', 'True'])
                if self.version.split('_')[1]=='deepTau2p5':
                    anaCacheTupleProducer_cmd.extend([ '--deepTauVersion', 'v2p5'])
                ps_call(anaCacheTupleProducer_cmd, env=self.cmssw_env, verbose=1)
            print(f"finished to produce anacachetuple")

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
        for prod_br, (sample_name, sample_type, branch) in anaProd_branch_map.items():
            if sample_type != "data": continue
            prod_branches.append(prod_br)
        return { 0: prod_branches }

    def output(self, force_pre_output=False):
        outFileName = 'nanoHTT_0.root'
        output_path = os.path.join('anaCache', self.period, 'data',self.version, outFileName)
        return self.remote_target(output_path, fs=self.fs_anaCache)

    def run(self):
        producer_dataMerge = os.path.join(self.ana_path(), 'AnaProd', 'MergeNtuples.py')
        with contextlib.ExitStack() as stack:
            local_inputs = [stack.enter_context(inp.localize('r')).path for inp in self.input()]
            with self.output().localize("w") as tmp_local_file:
                tmpFile = tmp_local_file.path
                dataMerge_cmd = [ 'python3', producer_dataMerge, '--outFile', tmpFile]
                dataMerge_cmd.extend(local_inputs)
                #print(dataMerge_cmd)
                ps_call(dataMerge_cmd,verbose=1)



