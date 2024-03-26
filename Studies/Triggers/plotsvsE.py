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

def make_plot(dict_to_plot, bckgs = '01', isGrav= True):
    color_key = 0
    colors_vse = ['red','blue','green']
    for final_key in dict_to_plot.keys():
        #actual_color_key = int((len(colors)-1)/color_key) if color_key != 0 else 0
        color = colors_vse[color_key]
        label = final_key.split('_')[-1]
        new_y_values = [a/b for a,b in zip(dict_to_plot[final_key][f"s/sqrtB_{bckgs}"],dict_to_plot['tautau_Medium_Medium_VVLoose'][f"s/sqrtB_{bckgs}"])]
        plt.plot(dict_to_plot[final_key]['masses'], new_y_values, color=color, label=label)
        color_key +=1
        #print(color_key)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    title = 'relative s/sqrt(b) - Radion'
    fileName = f"ssqrtb_rel_rad_all_bckgs{suffix}_{bckgs}_Medium.pdf"
    if isGrav :
        title = 'relative s/sqrt(b) - Graviton'
        fileName = f"ssqrtb_rel_grav_all_bckgs{suffix}_{bckgs}_Medium.pdf"
    plt.title(title)
    plt.ylabel("s/sqrt(b)")
    plt.xlabel(r'$m_X$ (GeV/$c^2$)')
    plt.subplots_adjust(right=0.7)
    plt.savefig(fileName)
    plt.clf()

def print_final_table(dict_to_consider):
    for key in dict_to_consider.keys():
        key_split = key.split('_')
        for idx in range(0, len(dict_to_consider[key]['masses'])):
            print(f"""{key_split[1]}\t{key_split[2]}\t{str(dict_to_consider[key]['masses'][idx]).replace('.',',')}\t{str(dict_to_consider[key]['s/sqrtB_01'][idx]).replace('.',',')}\t{str(dict_to_consider[key]['s/sqrtB_012'][idx]).replace('.',',')}""")


def GetSqrtBTot(all_sums, keys_to_consider =[]):
    dict_sqrtB = {}
    for key_wp in all_sums_dict[keys_to_consider[0]].keys():
        all_b_sums = 0
        for key_bckg in keys_to_consider:
            if key_wp not in all_sums_dict[key_bckg].keys():
                print(key_bckg, key_wp)
                continue
            all_b_sums+=all_sums_dict[key_bckg][key_wp]['sum']
        if key_wp not in dict_sqrtB:
            dict_sqrtB[key_wp] = {}
        dict_sqrtB[key_wp]['sum'] = math.sqrt(all_b_sums)
        #dict_sqrtB[key]['error'] = math.sqrt(all_sums_dict['TT'][key]+all_sums_dict['DY'][key])
    return dict_sqrtB

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--want2b', required=False, type=bool, default=False)
    args = parser.parse_args()
    colors_vse = ['red','blue','green']
    color_key = 0

    dict_S_grav = {}
    dict_S_rad = {}
    for vsesuffix in ['Loose', 'VLoose', 'VVLoose']:
        all_sums_dict = {
            'radion':{},
            'graviton':{},
            'TT_0':{},
            'TT_1':{},
            'TT_2':{},
            'DY_0':{},
            'DY_1':{},
            'DY_2':{},

            }

        wp_key = 'tautau_Medium_Medium'
        full_key = f'tautau_Medium_Medium_{vsesuffix}'
        print(full_key)
        if full_key not in dict_S_rad.keys():
            dict_S_rad[full_key] = {}

        suffix = '_2btag' if args.want2b == True else ''
        sig_files_bulk_graviton = f"signals_bulk_graviton{suffix}_VSe{vsesuffix}.json"
        sig_files_radion = f"signals_radion{suffix}_VSe{vsesuffix}.json"
        bckg_files_DY_0 = f"backgrounds_DY_0{suffix}_VSe{vsesuffix}.json"
        bckg_files_DY_1 = f"backgrounds_DY_1{suffix}_VSe{vsesuffix}.json"
        bckg_files_DY_2 = f"backgrounds_DY_2{suffix}_VSe{vsesuffix}.json"
        bckg_files_TT_0 = f"backgrounds_TT_0{suffix}_VSe{vsesuffix}.json"
        bckg_files_TT_1 = f"backgrounds_TT_1{suffix}_VSe{vsesuffix}.json"
        bckg_files_TT_2 = f"backgrounds_TT_2{suffix}_VSe{vsesuffix}.json"

        dict_signals_radion = {}
        dict_signals_graviton = {}
        dict_bckg_DY_0 = {}
        dict_bckg_TT_0 = {}

        dict_bckg_DY_1 = {}
        dict_bckg_TT_1 = {}

        dict_bckg_DY_2 = {}
        dict_bckg_TT_2 = {}

        with open(os.path.join(json_dir, sig_files_bulk_graviton), 'r') as sig_file_bulk_graviton:
            dict_signals_graviton=json.load(sig_file_bulk_graviton)
        with open(os.path.join(json_dir, sig_files_radion), 'r') as sig_file_radion:
            dict_signals_radion=json.load(sig_file_radion)

        with open(os.path.join(json_dir, bckg_files_DY_0), 'r') as bckg_file_DY_0:
            dict_bckg_DY_0=json.load(bckg_file_DY_0)
        with open(os.path.join(json_dir, bckg_files_TT_0), 'r') as bckg_file_TT_0:
            dict_bckg_TT_0=json.load(bckg_file_TT_0)

        with open(os.path.join(json_dir, bckg_files_DY_1), 'r') as bckg_file_DY_1:
            dict_bckg_DY_1=json.load(bckg_file_DY_1)
        with open(os.path.join(json_dir, bckg_files_TT_1), 'r') as bckg_file_TT_1:
            dict_bckg_TT_1=json.load(bckg_file_TT_1)

        with open(os.path.join(json_dir, bckg_files_DY_2), 'r') as bckg_file_DY_2:
            dict_bckg_DY_2=json.load(bckg_file_DY_2)
        with open(os.path.join(json_dir, bckg_files_TT_2), 'r') as bckg_file_TT_2:
            dict_bckg_TT_2=json.load(bckg_file_TT_2)
        #print(dict_signals_graviton)

        GetFinalDict(dict_signals_radion, all_sums_dict['radion'])
        GetFinalDict(dict_signals_graviton, all_sums_dict['graviton'])
        GetFinalDict(dict_bckg_TT_0, all_sums_dict['TT_0'])
        GetFinalDict(dict_bckg_DY_0, all_sums_dict['DY_0'])
        GetFinalDict(dict_bckg_TT_1, all_sums_dict['TT_1'])
        GetFinalDict(dict_bckg_DY_1, all_sums_dict['DY_1'])
        GetFinalDict(dict_bckg_TT_2, all_sums_dict['TT_2'])
        GetFinalDict(dict_bckg_DY_2, all_sums_dict['DY_2'])
        #print(all_sums_dict.keys())
        dict_sqrtB_012 = GetSqrtBTot(all_sums_dict, ['TT_0','TT_1','TT_2','DY_0','DY_1','DY_2'])
        dict_sqrtB_01 = GetSqrtBTot(all_sums_dict, ['TT_0','TT_1','DY_0','DY_1'])
        #print(dict_sqrtB)
        masses_points_grav = []
        masses_points_rad = []
        print(dict_S_rad.keys())
        #print(dict_signals_radion)
        for sample_key in dict_signals_radion[wp_key].keys():
            if sample_key == 'all' : continue
            #print( dict_signals_radion[wp_key][sample_key])
            mass_rad = sample_key.split('-')[1]
            if 'masses' not in dict_S_rad[full_key].keys():
                dict_S_rad[full_key]['masses'] = []
            dict_S_rad[full_key]['masses'].append(int(mass_rad))
            if 's/sqrtB_01' not in dict_S_rad[full_key].keys():
                dict_S_rad[full_key]['s/sqrtB_01'] = []
            if wp_key not in dict_signals_radion.keys():
                print(wp_key)
            ssqrtb_01=dict_signals_radion[wp_key][sample_key]['sum']/dict_sqrtB_01[wp_key]['sum']
            #dict_S_rad[full_key]['s/sqrtB_01'].append(ssqrtb_01/ssqrtb_01_ref)
            dict_S_rad[full_key]['s/sqrtB_01'].append(ssqrtb_01)

            if 's/sqrtB_012' not in dict_S_rad[full_key].keys():
                dict_S_rad[full_key]['s/sqrtB_012'] = []
            if wp_key not in dict_signals_radion.keys():
                print(wp_key)
            ssqrtb_012=dict_signals_radion[wp_key][sample_key]['sum']/dict_sqrtB_012[wp_key]['sum']
            #dict_S_rad[full_key]['s/sqrtB_012'].append(ssqrtb_012/ssqrtb_012_ref)
            dict_S_rad[full_key]['s/sqrtB_012'].append(ssqrtb_012)
    print(dict_S_rad.keys())



    #make_plot(dict_S_grav, bckgs = '01', isGrav= True)
    #make_plot(dict_S_grav, bckgs = '012', isGrav= True)
    make_plot(dict_S_rad, bckgs = '01', isGrav= False)
    make_plot(dict_S_rad, bckgs = '012', isGrav= False)
