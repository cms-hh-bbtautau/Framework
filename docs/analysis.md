# Common analysis steps

Remarks:

- commands bellow assume that `ERA` variable is set. E.g.
    ```sh
    ERA=Run2_2016
    ```
    Alternatively you can add `ERA=Run2_2016; ...` in front of each command.

- `version` argument alows to produce different versions of the same task. In the command below `--version dev` is used for illustration purposes. You can replace it with your version naming.
- `--workflow` can be `htcondor` or `local`. It is recommended to develop and test locally and then switch to `htcondor` for production. In examples below `--workflow local` is used for illustration purposes.
- when running on `htcondor` it is recommended to add `--transfer-logs` to the command to transfer logs to local.
- `--customisations` argument is used to pass custom parameters to the task in form param1=value1,param2=value2,...
    IMPORTANT for HHbbTauTau analysis: if running using deepTau 2p5 add `--customisations deepTauVersion=2p5`
- if you want to run only on few files, you can specify list of branches to run using `--branches` argument. E.g. `--branches 2,7-10,17`.
- to get status, use `--print-stauts N,K` where N is depth for task dependencies, K is depths for file dependencies. E.g. `--print-status 3,1`.
- to remove task output use `--remove-output N,a`, where N is depth for task dependencies. E.g. `--remove-output 0,a`.
- it is highly recommended to limitate the maximum number of parallel jobs running adding `--parallel-jobs M` where M is the number of the parallel jobs (e.g. M=100)

## Create input file list

```sh
law run InputFileTask  --period ${ERA} --version dev
```

## Create anaCache

```sh
law run AnaCacheTask  --period ${ERA} --version dev
```

## Create anaTuple

```sh
law run AnaTupleTask --period ${ERA} --version dev
```

## Merge data

```sh
law run DataMergeTask --period ${ERA} --version dev
```

- note: It's very important to first run `InputFileTask` then the other tasks dependencies are automatically fixed (e.g. if running `AnaTupleTask` without `AnaCacheTask`, it will first run `AnaCacheTask` then `AnaTupleTask`). If you do not run `InputFileTask`, running other tasks will raise error.