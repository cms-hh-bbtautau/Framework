# HH->bb$\tau$$\tau$ analysis steps

**Commands below assume that AnaTuples have already been produced. If not, please produce them following the instruction in the analysis section.**

Remember that:

- `ERA` variable is set. E.g.
    ```sh
    ERA=Run2_2016
    ```
    Alternatively you can add `ERA=Run2_2016; ...` in front of each command.
    Run2 possible eras are: `Run2_2016`,`Run2_2016_HIPM`,`Run2_2017` and `Run2_2018`
    <br/>
- when expliciting `VERSION_NAME` variable, its name contains explicitly the deepTau version: `VERSION_NAME= vXX_deepTauYY_ZZZ`, where:
    - XX is the anaTuple version (if not the first production it can be useful to have `v1,v2,..`),
    - YY is the deepTau version (`2p1` or `2p5`)
    - ZZZ are other eventual addition (e.g. if only tauTau channel `_onlyTauTau` or if `Zmumu` ntuples `_Zmumu`..)
    <br/>
- `--workflow` can be `htcondor` or `local`. It is recommended to develop and test locally and then switch to `htcondor` for production. In examples below `--workflow local` is used for illustration purposes.<br/> <br/>
- when running on `htcondor` it is recommended to add `--transfer-logs` to the command to transfer logs to local.<br/> <br/>
- `--customisations` argument is used to pass custom parameters to the task in form param1=value1,param2=value2,...
    **IMPORTANT for HHbbTauTau analysis:** if running using deepTau 2p5 add `--customisations deepTauVersion=2p5`<br/> <br/>
- if you want to run only on few files, you can specify list of branches to run using `--branches` argument. E.g. `--branches 2,7-10,17`.<br/> <br/>
- to get status, use `--print-stauts N,K` where N is depth for task dependencies, K is depths for file dependencies. E.g. `--print-status 3,1`.<br/> <br/>
- to remove task output use `--remove-output N,a`, where N is depth for task dependencies. E.g. `--remove-output 0,a`.<br/> <br/>
- it is highly recommended to limitate the maximum number of parallel jobs running adding `--parallel-jobs M` where M is the number of the parallel jobs (e.g. M=100)

## Create anaCacheTuple

For each Anatuple, an anaCacheTuple (storing observables which are computationally heavier) will be created.

```sh
law run AnaCacheTupleTask --period ${ERA} --version ${VERSION_NAME}
```
**Note**: at the `AnaCacheTupleTask` stage, the addition of customisation for specifying the version is still needed. For the other tasks, it won't be needed anymore.


#### Merge data in anaCache tuples

```sh
law run DataCacheMergeTask --period ${ERA} --version ${VERSION_NAME}
```


### Histograms Production

This has to be run after AnaTupleTask but **not necessairly** after AnaCacheTupleTask, if the variable to plot is not stored inside AnaCacheTuples.

These task will produce histograms with observables that need to be specified inside the `Analysis/tasks.py` file, specifically inside the `vars_to_plot` list.

The tasks to run are the following:

1. `HistProducerFileTask`: for each AnaTuple an histogram of the corresponding variable will be created.
    ```sh
    law run HistProducerFileTask --period $ERA --version ${VERSION_NAME}
    ```
1. `HistProducerSampleTask`: all the histogram belonging to a specific sample will be merged in one histogram.
    ```sh
    law run HistProducerSampleTask --period $ERA --version ${VERSION_NAME}
    ```
1. `MergeTask`: all the histogram will be merged from samples to only one histograms under the folder `${HISTOGRAMS}/all_histograms/` to a specific sample will be merged in one histogram. At this stage, for each norm/shape uncertainty (+ central scenario) will be created one histogram.
    ```sh
    law run MergeTask --period $ERA --version ${VERSION_NAME}
    ```
    Each histograms will be named as: `all_histograms_UNCERTAINTY.root` where uncertainty can be [Central, TauES_DM0, ecc....]

1. `HaddMergedTask`: all the merged histograms (produced separately for each uncertainty) will be merged in only one file.
    ```sh
    law run HaddMergedTask --period $ERA --version ${VERSION_NAME}
    ```
    Tip: It's very fast so it can be convenient to run this task in local.
    The final histogram will be named as: `all_histograms_Hadded.root`

## How to run HHbtag training skim ntuple production
```sh
python Studies/HHBTag/CreateTrainingSkim.py --inFile $CENTRAL_STORAGE/prod_v1/nanoAOD/2018/GluGluToBulkGravitonToHHTo2B2Tau_M-350.root --outFile output/skim.root --mass 350 --sample GluGluToBulkGraviton --year 2018 >& EventInfo.txt
python Common/SaveHisto.txt --inFile $CENTRAL_STORAGE/prod_v1/nanoAOD/2018/GluGluToBulkGravitonToHHTo2B2Tau_M-350.root --outFile output/skim.root
```
