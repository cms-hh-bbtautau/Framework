import json
import os
import sys
import math

if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])

json_dir = '/afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/Studies/Triggers/files/pt_reordered/'
#json_dir = '/afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/Studies/Triggers/files/'
def GetFinalDict(dict_to_consider, final_dict):
    for key in dict_to_consider.keys():
        final_dict[key] = {}
        sum_key = 0
        sum_error_key = 0
        for sig_key in dict_to_consider[key].keys():
            if sig_key == 'all': continue
            sum_key+=dict_to_consider[key][sig_key]['sum']
            sum_error_key+=math.pow(dict_to_consider[key][sig_key]['error'],2)

        final_dict[key]['sum'] = sum_key
        final_dict[key]['error'] = math.sqrt(sum_error_key)

def print_table(dict_to_consider):
    for key in dict_to_consider.keys():
        key_split = key.split('_')
        for sig_key in dict_to_consider[key].keys():
            if sig_key == 'all': continue
            print(f"""{key_split[1]}\t{key_split[2]}\t{sig_key}\t{str(dict_to_consider[key][sig_key]["sum"]).replace('.',',')}\t{str(dict_to_consider[key][sig_key]["error"]).replace('.',',')}""")

def print_final_table(dict_to_consider):
    for key in dict_to_consider.keys():
        for sig_key in dict_to_consider[key].keys():
            sig_key_split = sig_key.split('_')
            print(f"""{sig_key_split[1]}\t{sig_key_split[2]}\t{key}\t{str(dict_to_consider[key][sig_key]["sum"]).replace('.',',')}\t{str(dict_to_consider[key][sig_key]["error"]).replace('.',',')}""")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--want2b', required=False, type=bool, default=False)
    parser.add_argument('--bckg_sel', required=False, type=str, default="0")
    args = parser.parse_args()

    all_sums_dict = {
        'radion':{},
        'graviton':{},
        'TT':{},
        'DY':{},

        }


    suffix = '_2btag' if args.want2b == True else ''
    sig_files_bulk_graviton = f"signals_bulk_graviton{suffix}_ptReordered.json"
    sig_files_radion = f"signals_radion{suffix}_ptReordered.json"
    bckg_files_DY = f"backgrounds_DY_{args.bckg_sel}{suffix}_ptReordered.json"
    bckg_files_TT = f"backgrounds_TT_{args.bckg_sel}{suffix}_ptReordered.json"

    #sig_files_bulk_graviton = f"signals_bulk_graviton{suffix}.json"
    #sig_files_radion = f"signals_radion{suffix}.json"
    #bckg_files_DY = f"backgrounds_DY_{args.bckg_sel}{suffix}.json"
    #bckg_files_TT = f"backgrounds_TT_{args.bckg_sel}{suffix}.json"

    dict_signals_radion = {}
    dict_signals_graviton = {}
    dict_bckg_DY = {}
    dict_bckg_TT = {}

    with open(os.path.join(json_dir, sig_files_bulk_graviton), 'r') as sig_file_bulk_graviton:
        dict_signals_graviton=json.load(sig_file_bulk_graviton)
    with open(os.path.join(json_dir, sig_files_radion), 'r') as sig_file_radion:
        dict_signals_radion=json.load(sig_file_radion)
    with open(os.path.join(json_dir, bckg_files_DY), 'r') as bckg_file_DY:
        dict_bckg_DY=json.load(bckg_file_DY)
    with open(os.path.join(json_dir, bckg_files_TT), 'r') as bckg_file_TT:
        dict_bckg_TT=json.load(bckg_file_TT)


    GetFinalDict(dict_signals_radion, all_sums_dict['radion'])
    GetFinalDict(dict_signals_graviton, all_sums_dict['graviton'])
    GetFinalDict(dict_bckg_TT, all_sums_dict['TT'])
    GetFinalDict(dict_bckg_DY, all_sums_dict['DY'])
    print_table(dict_signals_radion)
    print_table(dict_signals_graviton)
    print_table(dict_bckg_DY)
    print_table(dict_bckg_TT)
    print()
    print_final_table(all_sums_dict)
    with open(os.path.join(json_dir, f"all_sums_final_{args.bckg_sel}{suffix}.json"),'w') as final_file:
        json.dump(all_sums_dict, final_file, indent=4)
