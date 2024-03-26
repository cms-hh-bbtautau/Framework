import os
import sys
import yaml
import ROOT
import datetime

if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])

import Common.BaselineSelection as Baseline
import Corrections.Corrections as Corrections
import Common.Utilities as Utilities
from Corrections.CorrectionsCore import *
from Corrections.pu import *

def createAnaCache(inDir, outFile, global_params, isData, range=None, verbose=1):
    start_time = datetime.datetime.now()
    Baseline.Initialize(False, False)
    Corrections.Initialize(config=global_params, isData=isData,  load_corr_lib=True, load_pu=True, load_tau=False, load_trg=False, load_btag=False, loadBTagEff=False, load_met=False, load_mu = False, load_ele=False, load_puJetID=False, load_jet=False,load_fatjet=False)

    if os.path.exists(args.outFile):
        with open(args.outFile, 'r') as file:
            dict = yaml.safe_load(file)
    else:
        dict = {}
    dict['denominator']={}
    sources = [ central ] + puWeightProducer.uncSource
    input_files = []
    for f in os.listdir(inDir):
        if f.split('.')[1] != 'root': continue
        input_files.append(os.path.join(inDir, f))
    for tree in [ 'Events', 'EventsNotSelected']:
        df = ROOT.RDataFrame(tree, input_files)
        if range is not None:
            df = df.Range(range)
        df,syst_names = Corrections.getDenominator(df, sources)
        for source in sources:
            if source not in dict['denominator']:
                dict['denominator'][source]={}
            for scale in getScales(source):
                syst_name = getSystName(source, scale)
                dict['denominator'][source][scale] = dict['denominator'][source].get(scale, 0.) + df.Sum(f'weight_denom_{syst_name}').GetValue()
    end_time = datetime.datetime.now()
    dict['runtime'] = (end_time - start_time).total_seconds()
    if verbose > 0:
        print(dict)
    with open(outFile, 'w') as file:
        yaml.dump(dict, file)

if __name__ == "__main__":
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True, type=str)
    parser.add_argument('--inDir', required=True, type=str)
    parser.add_argument('--sample', required=True, type=str)
    parser.add_argument('--outFile', required=True, type=str)
    parser.add_argument('--nEvents', type=int, default=None)
    parser.add_argument('--customisations', type=str, default="")
    args = parser.parse_args()

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    if len(args.customisations) > 0:
        Utilities.ApplyConfigCustomisations(config['GLOBAL'], args.customisations)
    print(config[args.sample]['sampleType'])
    isData = True if config[args.sample]['sampleType'] == 'data' else False
    print(isData)
    createAnaCache(args.inDir, args.outFile, config['GLOBAL'],isData, range=args.nEvents)
