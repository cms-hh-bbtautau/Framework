import json
import os
import sys
import math

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

if __name__ == "__main__":
    sys.path.append(os.environ['ANALYSIS_PATH'])

colors = ['black', 'k', 'dimgray', 'dimgrey', 'gray', 'grey', 'darkgray', 'darkgrey', 'silver', 'lightgray', 'lightgrey', 'gainsboro', 'whitesmoke', 'w', 'white', 'snow', 'rosybrown', 'lightcoral', 'indianred', 'brown', 'firebrick', 'maroon', 'darkred', 'r', 'red', 'mistyrose', 'salmon', 'tomato', 'darksalmon', 'coral', 'orangered', 'lightsalmon', 'sienna', 'seashell', 'chocolate', 'saddlebrown', 'sandybrown', 'peachpuff', 'peru', 'linen', 'bisque', 'darkorange', 'burlywood', 'antiquewhite', 'tan', 'navajowhite', 'blanchedalmond', 'papayawhip', 'moccasin', 'orange', 'wheat', 'oldlace', 'floralwhite', 'darkgoldenrod', 'goldenrod', 'cornsilk', 'gold', 'lemonchiffon', 'khaki', 'palegoldenrod', 'darkkhaki', 'ivory', 'beige', 'lightyellow', 'lightgoldenrodyellow', 'olive', 'y', 'yellow', 'olivedrab', 'yellowgreen', 'darkolivegreen', 'greenyellow', 'chartreuse', 'lawngreen', 'honeydew', 'darkseagreen', 'palegreen', 'lightgreen', 'forestgreen', 'limegreen', 'darkgreen', 'g', 'green', 'lime', 'seagreen', 'mediumseagreen', 'springgreen', 'mintcream', 'mediumspringgreen', 'mediumaquamarine', 'aquamarine', 'turquoise', 'lightseagreen', 'mediumturquoise', 'azure', 'lightcyan', 'paleturquoise', 'darkslategray', 'darkslategrey', 'teal', 'darkcyan', 'c', 'aqua', 'cyan', 'darkturquoise', 'cadetblue', 'powderblue', 'lightblue', 'deepskyblue', 'skyblue', 'lightskyblue', 'steelblue', 'aliceblue', 'dodgerblue', 'lightslategray', 'lightslategrey', 'slategray', 'slategrey', 'lightsteelblue', 'cornflowerblue', 'royalblue', 'ghostwhite', 'lavender', 'midnightblue', 'navy', 'darkblue', 'mediumblue', 'b', 'blue', 'slateblue', 'darkslateblue', 'mediumslateblue', 'mediumpurple', 'rebeccapurple', 'blueviolet', 'indigo', 'darkorchid', 'darkviolet', 'mediumorchid', 'thistle', 'plum', 'violet', 'purple', 'darkmagenta', 'm', 'fuchsia', 'magenta', 'orchid', 'mediumvioletred', 'deeppink', 'hotpink', 'lavenderblush', 'palevioletred', 'crimson', 'pink', 'lightpink']



actual_colors = ['black',  'silver', 'lightgray', 'lightcoral', 'maroon', 'mistyrose', 'salmon', 'olivedrab',  'limegreen', 'darkgreen','aquamarine', 'azure', 'cornflowerblue',  'navy', 'mediumslateblue', 'indigo', 'darkorchid', 'darkviolet', 'mediumorchid', 'thistle', 'violet', 'purple', 'darkmagenta', 'm', 'fuchsia', 'magenta', 'orchid', 'mediumvioletred', 'deeppink', 'hotpink', 'lavenderblush', 'palevioletred', 'crimson', 'pink', 'lightpink']


json_dir = '/afs/cern.ch/work/v/vdamante/hhbbTauTauRes/prod/Framework/Studies/Triggers/files/'
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


def GetSqrtBTot(all_sums):
    dict_sqrtB = {}
    for key in all_sums_dict['TT'].keys():
        if key not in all_sums_dict['DY']:
            print(f"{key} not in DY contrib!")
        if key not in dict_sqrtB:
            dict_sqrtB[key] = {}
        dict_sqrtB[key]['sum'] = math.sqrt(all_sums_dict['TT'][key]['sum']+all_sums_dict['DY'][key]['sum'])
        #dict_sqrtB[key]['error'] = math.sqrt(all_sums_dict['TT'][key]+all_sums_dict['DY'][key])
    return dict_sqrtB

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
    sig_files_bulk_graviton = f"signals_bulk_graviton{suffix}_VSeVLoose_onlyElectrons.json"
    sig_files_radion = f"signals_radion{suffix}_VSeVLoose_onlyElectrons.json"
    bckg_files_DY = f"backgrounds_DY_{args.bckg_sel}{suffix}_VSeVLoose_onlyElectrons.json"
    bckg_files_TT = f"backgrounds_TT_{args.bckg_sel}{suffix}_VSeVLoose_onlyElectrons.json"

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

    #print(dict_signals_graviton)

    GetFinalDict(dict_signals_radion, all_sums_dict['radion'])
    GetFinalDict(dict_signals_graviton, all_sums_dict['graviton'])
    GetFinalDict(dict_bckg_TT, all_sums_dict['TT'])
    GetFinalDict(dict_bckg_DY, all_sums_dict['DY'])

    #print_table(dict_signals_radion)
    #print_table(dict_signals_graviton)
    #print_table(dict_bckg_DY)
    #print_table(dict_bckg_TT)
    print()
    print_final_table(all_sums_dict)
    with open(os.path.join(json_dir, f"all_sums_final_{args.bckg_sel}{suffix}_VSeVLoose_onlyElectrons.json"),'w') as final_file:
        json.dump(all_sums_dict, final_file, indent=4)
    '''
    dict_sqrtB = GetSqrtBTot(all_sums_dict)
    #print(dict_sqrtB)
    dict_S_grav = {}
    dict_S_rad = {}
    masses_points_grav = []
    masses_points_rad = []
    wp_key_ref = 'tautau_Medium_Medium'
    for wp_key in dict_signals_graviton.keys():
        if wp_key not in dict_S_grav.keys():
            dict_S_grav[wp_key] = {}
        for sample_key in dict_signals_graviton[wp_key].keys():
            if sample_key == 'all' : continue
            #print( dict_signals_graviton[wp_key][sample_key])
            mass_grav = sample_key.split('-')[1]
            if 'masses' not in dict_S_grav[wp_key].keys():
                dict_S_grav[wp_key]['masses'] = []
            if 's/sqrtB' not in dict_S_grav[wp_key].keys():
                dict_S_grav[wp_key]['s/sqrtB'] = []
            dict_S_grav[wp_key]['masses'].append(int(mass_grav))

            ssqrtb=dict_signals_graviton[wp_key][sample_key]['sum']/dict_sqrtB[wp_key]['sum']
            ssqrtb_ref=dict_signals_graviton[wp_key_ref][sample_key]['sum']/dict_sqrtB[wp_key_ref]['sum']
            dict_S_grav[wp_key]['s/sqrtB'].append(ssqrtb/ssqrtb_ref)
        if wp_key not in dict_S_rad.keys():
            dict_S_rad[wp_key] = {}
        for sample_key in dict_signals_radion[wp_key]:

            if sample_key == 'all' : continue
            #print(sample_key)
            mass_rad = sample_key.split('-')[1]
            if 'masses' not in dict_S_rad[wp_key].keys():
                dict_S_rad[wp_key]['masses'] = []
            if 's/sqrtB' not in dict_S_rad[wp_key].keys():
                dict_S_rad[wp_key]['s/sqrtB'] = []
            dict_S_rad[wp_key]['masses'].append(int(mass_rad))
            ssqrtb=dict_signals_radion[wp_key][sample_key]['sum']/dict_sqrtB[wp_key]['sum']
            ssqrtb_ref=dict_signals_radion[wp_key_ref][sample_key]['sum']/dict_sqrtB[wp_key_ref]['sum']
            dict_S_rad[wp_key]['s/sqrtB'].append(ssqrtb/ssqrtb_ref)
    #print(dict_S_grav)
    #print(dict_S_rad)
    #colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    #print(colors)
    #by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
    #            for name, color in colors.items())
    #sorted_names = [name for hsv, name in by_hsv]
    #print(sorted_names)

    color_key = 0
    for wp_key in dict_S_grav.keys():
        #actual_color_key = int((len(colors)-1)/color_key) if color_key != 0 else 0
        actual_color_key = int(color_key * 2 + 0.5 * color_key/2)
        label_split = wp_key.split('_')
        label = f"{label_split[1]}_{label_split[2]}"
        plt.plot(dict_S_grav[wp_key]['masses'], dict_S_grav[wp_key]['s/sqrtB'], color=actual_colors[color_key], label=label)
        color_key +=1
        #print(color_key)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xlabel("M_{X}")
    plt.ylabel("s/sqrt(b)")
    plt.title("s/sqrt(b)")
    plt.subplots_adjust(right=0.7)
    plt.savefig(f"ssqrtb_rel_grav_{args.bckg_sel}{suffix}.pdf")
    plt.clf()
    color_key = 0
    for wp_key in dict_S_rad.keys():
        #actual_color_key = int((len(colors)-1)/color_key) if color_key != 0 else 0
        actual_color_key = int(color_key * 2 + 0.5 * color_key/2)
        label_split = wp_key.split('_')
        label = f"{label_split[1]}_{label_split[2]}"
        plt.plot(dict_S_rad[wp_key]['masses'], dict_S_rad[wp_key]['s/sqrtB'], color=actual_colors[color_key], label=label)
        color_key +=1
        #print(color_key)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xlabel("M_{X}")
    plt.ylabel("s/sqrt(b)")
    plt.title("s/sqrt(b)")
    plt.subplots_adjust(right=0.7)
    plt.savefig(f"ssqrtb_rel_rad_{args.bckg_sel}{suffix}.pdf")
    '''