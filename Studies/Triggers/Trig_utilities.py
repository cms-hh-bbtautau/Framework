
import matplotlib.pyplot as plt
import math

tau_pt_limits_dict = {
    "2016":
    {
        "eTau":[26,26,26],
        "muTau":[25, 20, 25],
        "tauTau":[40, 40, 40],
    },
    "2017":
    {
        "eTau":[33, 25, 35],
        "muTau":[28, 21, 32],
        "tauTau":[40, 40, 40],
    },
    "2018":
    {
        "eTau":[34, 26, 35],
        "muTau":[26, 22, 32],
        "tauTau":[40, 40, 40],
    }
}

ggR_samples = [
                "GluGluToRadionToHHTo2B2Tau_M-250", "GluGluToRadionToHHTo2B2Tau_M-260", "GluGluToRadionToHHTo2B2Tau_M-270", "GluGluToRadionToHHTo2B2Tau_M-280", "GluGluToRadionToHHTo2B2Tau_M-300", "GluGluToRadionToHHTo2B2Tau_M-320", "GluGluToRadionToHHTo2B2Tau_M-350", "GluGluToRadionToHHTo2B2Tau_M-400", "GluGluToRadionToHHTo2B2Tau_M-450", "GluGluToRadionToHHTo2B2Tau_M-500", "GluGluToRadionToHHTo2B2Tau_M-550", "GluGluToRadionToHHTo2B2Tau_M-600", "GluGluToRadionToHHTo2B2Tau_M-650", "GluGluToRadionToHHTo2B2Tau_M-700", "GluGluToRadionToHHTo2B2Tau_M-750", "GluGluToRadionToHHTo2B2Tau_M-800", "GluGluToRadionToHHTo2B2Tau_M-850", "GluGluToRadionToHHTo2B2Tau_M-900" , "GluGluToRadionToHHTo2B2Tau_M-1000", "GluGluToRadionToHHTo2B2Tau_M-1250", "GluGluToRadionToHHTo2B2Tau_M-1500", "GluGluToRadionToHHTo2B2Tau_M-1750", "GluGluToRadionToHHTo2B2Tau_M-2000", "GluGluToRadionToHHTo2B2Tau_M-2500","GluGluToRadionToHHTo2B2Tau_M-3000" ]

ggBG_samples = [ # missing mass 400
    'GluGluToBulkGravitonToHHTo2B2Tau_M-250', 'GluGluToBulkGravitonToHHTo2B2Tau_M-260', 'GluGluToBulkGravitonToHHTo2B2Tau_M-270', 'GluGluToBulkGravitonToHHTo2B2Tau_M-280', 'GluGluToBulkGravitonToHHTo2B2Tau_M-300', 'GluGluToBulkGravitonToHHTo2B2Tau_M-320', 'GluGluToBulkGravitonToHHTo2B2Tau_M-350','GluGluToBulkGravitonToHHTo2B2Tau_M-450', 'GluGluToBulkGravitonToHHTo2B2Tau_M-500', 'GluGluToBulkGravitonToHHTo2B2Tau_M-550', 'GluGluToBulkGravitonToHHTo2B2Tau_M-600', 'GluGluToBulkGravitonToHHTo2B2Tau_M-650', 'GluGluToBulkGravitonToHHTo2B2Tau_M-700', 'GluGluToBulkGravitonToHHTo2B2Tau_M-750', 'GluGluToBulkGravitonToHHTo2B2Tau_M-800', 'GluGluToBulkGravitonToHHTo2B2Tau_M-850', 'GluGluToBulkGravitonToHHTo2B2Tau_M-900', 'GluGluToBulkGravitonToHHTo2B2Tau_M-1000', 'GluGluToBulkGravitonToHHTo2B2Tau_M-1250', 'GluGluToBulkGravitonToHHTo2B2Tau_M-1500', 'GluGluToBulkGravitonToHHTo2B2Tau_M-1750', 'GluGluToBulkGravitonToHHTo2B2Tau_M-2000', 'GluGluToBulkGravitonToHHTo2B2Tau_M-2500','GluGluToBulkGravitonToHHTo2B2Tau_M-3000']


def get_regions_dict_ptreordered(channel, singleTau_pt_limits=[190,190], met_th=180, bigtau=False):
    regions = {}
    tau1new_eta_sel = {
        "eTau": "abs(tau1new_eta) <= 2.5",
        "muTau": "abs(tau1new_eta) <= 2.4",
        "tauTau": "abs(tau1new_eta)<=2.1"
    }
    tau1new_2p1 = "abs(tau1new_eta)<=2.1"
    tau2new_2p1 = "abs(tau2new_eta)<=2.1"
    tau_pt_limits = tau_pt_limits_dict["2018"][channel]
    first_tau_sel = "tau1new_pt >= {0} && {1}"
    second_tau_sel = "tau2new_pt >= {0} && {1}"

    first_tau_sel_singleTau = first_tau_sel.format(singleTau_pt_limits[0], tau1new_eta_sel[channel])
    second_tau_sel_singleTau = second_tau_sel.format(singleTau_pt_limits[1], tau2new_2p1)

    first_tau_sel_other = first_tau_sel.format(tau_pt_limits[0], tau1new_eta_sel[channel])
    second_tau_sel_other = second_tau_sel.format(tau_pt_limits[1], tau2new_2p1)
    single_lepton_validity = first_tau_sel.format(tau_pt_limits[0], tau1new_eta_sel[channel])
    cross_lepton_validity_first = first_tau_sel.format(tau_pt_limits[1], tau1new_2p1)
    cross_lepton_validity_second = second_tau_sel.format(tau_pt_limits[2], tau2new_2p1)
    cross_lepton_validity  = f"{cross_lepton_validity_first} && {cross_lepton_validity_second}"


    if bigtau:
        if channel == "tauTau":
            #second_tau_sel_singleTau = second_tau_sel.format(singleTau_pt_limits[1], tau2_2p1)
            #first_tau_sel_singleTau = first_tau_sel.format(singleTau_pt_limits[0], tau1_eta_sel[channel])
            regions["singleTau_region"] = f"""(({first_tau_sel_singleTau}) || ({second_tau_sel_singleTau}))"""

            #first_tau_sel_other = first_tau_sel.format(tau_pt_limits[0], tau1_eta_sel[channel])
            #second_tau_sel_other = second_tau_sel.format(tau_pt_limits[1], tau2_2p1)
            regions["other_trg_region"] = f"""(({first_tau_sel_other}) &&( {second_tau_sel_other})) && !({regions["singleTau_region"]})"""
        else:
            #second_tau_sel_singleTau = second_tau_sel.format(singleTau_pt_limits[1], tau2_2p1)
            regions["singleTau_region"] = second_tau_sel_singleTau

            #single_lepton_validity = first_tau_sel.format(tau_pt_limits[0], tau1_eta_sel[channel])
            #cross_lepton_validity_first = first_tau_sel.format(tau_pt_limits[1], tau1_2p1)
            #cross_lepton_validity_second = second_tau_sel.format(tau_pt_limits[2], tau2_2p1)
            #cross_lepton_validity  = f"{cross_lepton_validity_first} && {cross_lepton_validity_second}"
            regions['other_trg_region'] = f"""(({single_lepton_validity}) || ({cross_lepton_validity})) && !({regions["singleTau_region"]})"""

    else:
        if channel == "tauTau":
            #first_tau_sel_other = first_tau_sel.format(tau_pt_limits[0], tau1_eta_sel[channel])
            #second_tau_sel_other = second_tau_sel.format(tau_pt_limits[1], tau2_2p1)
            #second_tau_sel_singleTau = second_tau_sel.format(singleTau_pt_limits[1], tau2_2p1)
            #first_tau_sel_singleTau = first_tau_sel.format(singleTau_pt_limits[0], tau1_eta_sel[channel])
            regions["other_trg_region"] = f"""(({first_tau_sel_other}) && ({second_tau_sel_other}))"""
            regions["singleTau_region"] = f"""(({first_tau_sel_singleTau}) || ({second_tau_sel_singleTau}))&& !({regions["other_trg_region"]})"""

        else:
            #single_lepton_validity = first_tau_sel.format(tau_pt_limits[0], tau1_eta_sel[channel])
            #cross_lepton_validity_first = first_tau_sel.format(tau_pt_limits[1], tau1_2p1)
            #cross_lepton_validity_second = second_tau_sel.format(tau_pt_limits[2], tau2_2p1)
            #cross_lepton_validity  = f"{cross_lepton_validity_first} && {cross_lepton_validity_second}"
            #second_tau_sel_singleTau = second_tau_sel.format(singleTau_pt_limits[1], tau2_2p1)
            regions["other_trg_region"] = f"""(({single_lepton_validity}) || ({cross_lepton_validity}))"""
            regions["singleTau_region"] = f"""({second_tau_sel_singleTau}) &&!({regions["other_trg_region"]}) """

    regions["MET_region"] = f"""(met_pt > {met_th} && !( ({regions["other_trg_region"]}) || ( {regions["singleTau_region"]} ) ) )"""
    return regions




def get_regions_dict(channel, singleTau_pt_limits=[190,190], met_th=180, bigtau=False):
    regions = {}
    tau1_eta_sel = {
        "eTau": "abs(tau1_eta) <= 2.5",
        "muTau": "abs(tau1_eta) <= 2.4",
        "tauTau": "abs(tau1_eta)<=2.1"
    }
    tau1_2p1 = "abs(tau1_eta)<=2.1"
    tau2_2p1 = "abs(tau2_eta)<=2.1"
    tau_pt_limits = tau_pt_limits_dict["2018"][channel]
    first_tau_sel = "tau1_pt >= {0} && {1}"
    second_tau_sel = "tau2_pt >= {0} && {1}"

    first_tau_sel_singleTau = first_tau_sel.format(singleTau_pt_limits[0], tau1_eta_sel[channel])
    second_tau_sel_singleTau = second_tau_sel.format(singleTau_pt_limits[1], tau2_2p1)

    first_tau_sel_other = first_tau_sel.format(tau_pt_limits[0], tau1_eta_sel[channel])
    second_tau_sel_other = second_tau_sel.format(tau_pt_limits[1], tau2_2p1)
    single_lepton_validity = first_tau_sel.format(tau_pt_limits[0], tau1_eta_sel[channel])
    cross_lepton_validity_first = first_tau_sel.format(tau_pt_limits[1], tau1_2p1)
    cross_lepton_validity_second = second_tau_sel.format(tau_pt_limits[2], tau2_2p1)
    cross_lepton_validity  = f"{cross_lepton_validity_first} && {cross_lepton_validity_second}"


    if bigtau:
        if channel == "tauTau":
            #second_tau_sel_singleTau = second_tau_sel.format(singleTau_pt_limits[1], tau2_2p1)
            #first_tau_sel_singleTau = first_tau_sel.format(singleTau_pt_limits[0], tau1_eta_sel[channel])
            regions["singleTau_region"] = f"""(({first_tau_sel_singleTau}) || ({second_tau_sel_singleTau}))"""

            #first_tau_sel_other = first_tau_sel.format(tau_pt_limits[0], tau1_eta_sel[channel])
            #second_tau_sel_other = second_tau_sel.format(tau_pt_limits[1], tau2_2p1)
            regions["other_trg_region"] = f"""(({first_tau_sel_other}) &&( {second_tau_sel_other})) && !({regions["singleTau_region"]})"""
        else:
            #second_tau_sel_singleTau = second_tau_sel.format(singleTau_pt_limits[1], tau2_2p1)
            regions["singleTau_region"] = second_tau_sel_singleTau

            #single_lepton_validity = first_tau_sel.format(tau_pt_limits[0], tau1_eta_sel[channel])
            #cross_lepton_validity_first = first_tau_sel.format(tau_pt_limits[1], tau1_2p1)
            #cross_lepton_validity_second = second_tau_sel.format(tau_pt_limits[2], tau2_2p1)
            #cross_lepton_validity  = f"{cross_lepton_validity_first} && {cross_lepton_validity_second}"
            regions['other_trg_region'] = f"""(({single_lepton_validity}) || ({cross_lepton_validity})) && !({regions["singleTau_region"]})"""

    else:
        if channel == "tauTau":
            #first_tau_sel_other = first_tau_sel.format(tau_pt_limits[0], tau1_eta_sel[channel])
            #second_tau_sel_other = second_tau_sel.format(tau_pt_limits[1], tau2_2p1)
            #second_tau_sel_singleTau = second_tau_sel.format(singleTau_pt_limits[1], tau2_2p1)
            #first_tau_sel_singleTau = first_tau_sel.format(singleTau_pt_limits[0], tau1_eta_sel[channel])
            regions["other_trg_region"] = f"""(({first_tau_sel_other}) && ({second_tau_sel_other}))"""
            regions["singleTau_region"] = f"""(({first_tau_sel_singleTau}) || ({second_tau_sel_singleTau}))&& !({regions["other_trg_region"]})"""

        else:
            #single_lepton_validity = first_tau_sel.format(tau_pt_limits[0], tau1_eta_sel[channel])
            #cross_lepton_validity_first = first_tau_sel.format(tau_pt_limits[1], tau1_2p1)
            #cross_lepton_validity_second = second_tau_sel.format(tau_pt_limits[2], tau2_2p1)
            #cross_lepton_validity  = f"{cross_lepton_validity_first} && {cross_lepton_validity_second}"
            #second_tau_sel_singleTau = second_tau_sel.format(singleTau_pt_limits[1], tau2_2p1)
            regions["other_trg_region"] = f"""(({single_lepton_validity}) || ({cross_lepton_validity}))"""
            regions["singleTau_region"] = f"""({second_tau_sel_singleTau}) &&!({regions["other_trg_region"]}) """

    regions["MET_region"] = f"""(met_pt > {met_th} && !( ({regions["other_trg_region"]}) || ( {regions["singleTau_region"]} ) ) )"""
    return regions


def GetTrgValidityRegion(channel, trg, singleTau_pt_limits=[190,190], met_th=150, year="2018"):
    trg_validity_regions = {
        "HLT_ditau":f"""(tau1_pt>{tau_pt_limits_dict[year]["tauTau"][0]} && abs(tau1_eta) < 2.1) && (tau2_pt>{tau_pt_limits_dict[year]["tauTau"][1]} && abs(tau2_eta) < 2.1)""",
        "HLT_singleTau": f"""(tau2_pt>{singleTau_pt_limits[1]} && abs(tau2_eta) < 2.1 )""",
        "HLT_MET": f"metnomu_pt > {met_th}",
        "HLT_singleMu":f"""tau1_pt>{tau_pt_limits_dict[year]["muTau"][0]}""",
        "HLT_singleMu50":"tau1_pt>52",
        "HLT_mutau":f"""tau1_pt > {tau_pt_limits_dict[year]["muTau"][1]} && abs(tau1_eta) < 2.1 && tau2_pt > {tau_pt_limits_dict[year]["muTau"][2]} && abs(tau2_eta) < 2.1""",
        "HLT_singleEle":f"""tau1_pt>{tau_pt_limits_dict[year]["eTau"][0]} && abs(tau1_eta)<2.5""",
        "HLT_etau":f"""tau1_pt > {tau_pt_limits_dict[year]["eTau"][1]} && abs(tau1_eta) < 2.1 && tau2_pt > {tau_pt_limits_dict[year]["eTau"][2]} && abs(tau2_eta) < 2.1""",
    }
    if channel == "tauTau" and trg == "HLT_singleTau":
        trg_validity_regions[trg] = f"(tau1_pt > {singleTau_pt_limits[0]} && abs(tau1_eta) < 2.1) || (tau2_pt>{singleTau_pt_limits[1]} && abs(tau2_eta) < 2.1 )"
    return trg_validity_regions[trg]



def GetTrgValidityRegion_ptreordered(channel, trg, singleTau_pt_limits=[190,190], met_th=150, year="2018"):
    trg_validity_regions = {
        "HLT_ditau":f"""(tau1new_pt>{tau_pt_limits_dict[year]["tauTau"][0]} && abs(tau1new_eta) < 2.1) && (tau2new_pt>{tau_pt_limits_dict[year]["tauTau"][1]} && abs(tau2new_eta) < 2.1)""",
        "HLT_singleTau": f"""(tau2new_pt>{singleTau_pt_limits[1]} && abs(tau2new_eta) < 2.1 )""",
        "HLT_MET": f"metnomu_pt > {met_th}",
        "HLT_singleMu":f"""tau1new_pt>{tau_pt_limits_dict[year]["muTau"][0]}""",
        "HLT_singleMu50":"tau1new_pt>52",
        "HLT_mutau":f"""tau1new_pt > {tau_pt_limits_dict[year]["muTau"][1]} && abs(tau1new_eta) < 2.1 && tau2new_pt > {tau_pt_limits_dict[year]["muTau"][2]} && abs(tau2new_eta) < 2.1""",
        "HLT_singleEle":f"""tau1new_pt>{tau_pt_limits_dict[year]["eTau"][0]} && abs(tau1new_eta)<2.5""",
        "HLT_etau":f"""tau1new_pt > {tau_pt_limits_dict[year]["eTau"][1]} && abs(tau1new_eta) < 2.1 && tau2new_pt > {tau_pt_limits_dict[year]["eTau"][2]} && abs(tau2new_eta) < 2.1""",
    }
    if channel == "tauTau" and trg == "HLT_singleTau":
        trg_validity_regions[trg] = f"(tau1new_pt > {singleTau_pt_limits[0]} && abs(tau1new_eta) < 2.1) || (tau2new_pt>{singleTau_pt_limits[1]} && abs(tau2new_eta) < 2.1 )"
    return trg_validity_regions[trg]


def GetSum(df,channel, trg_reg_list, weight_expression, weight_dict, weight_key, verbose=False):
    filter_expr = f' ('
    filter_expr+= ' || '.join(trg for trg in trg_reg_list)
    filter_expr+= ')'
    sum_channel = df.Define("weight_for_sum",weight_expression).Filter(filter_expr).Sum("weight_for_sum").GetValue()
    return sum_channel
    '''
    if weight_key not in weight_dict.keys():
        weight_dict[weight_key] = 0.
    weight_dict[weight_key]+=round(sum_channel,5)
    '''

def AddEfficiencyToDict(df,trg_list, channel, n_initial_channel,eff_key, eff_dict, eff_dict_errors,wantNEvents=False, verbose=False):
    filter_expr = f' ('
    filter_expr+= ' || '.join(trg for trg in trg_list)
    filter_expr+= ')'
    #total_weight_string = 'weight_total*weight_TauID_Central*weight_tau1_EleidSF_Central*weight_tau1_MuidSF_Central*weight_tau2_EleidSF_Central*weight_tau2_MuidSF_Central*weight_L1PreFiring_Central*weight_L1PreFiring_ECAL_Central*weight_L1PreFiring_Muon_Central'
    #n_channel = df.Define('total_weight_string', total_weight_string).Filter(filter_expr).Sum('total_weight_string').GetValue()
    #from statsmodels.stats.proportion import proportion_confint
    n_channel = df.Filter(filter_expr).Sum('total_weight_string').GetValue()
    eff_channel = n_channel / n_initial_channel
    conf_interval = [0,0]
    #conf_interval = proportion_confint(n_channel, n_initial_channel, alpha=0.05, method='beta')
    error_low = eff_channel - conf_interval[0]
    error_up = conf_interval[1] - eff_channel
    if wantNEvents:
        eff_channel = n_channel
        error_low = math.sqrt(eff_channel)
        error_up =  math.sqrt(eff_channel)
    if verbose:
        print(eff_key)
        #print(f"with {filter_expr} : n_initial{channel} = {n_initial_channel}, n_{channel} = {n_channel}, eff_{channel} = {round(eff_channel,4)}")
        #print(f"with {filter_expr}")
        print(f"denum \t num  \t eff ")
        print(f"{n_initial_channel} \t  = {n_channel} \t {round(eff_channel,4)}")
    if eff_key not in eff_dict.keys():
        eff_dict[eff_key] = []
        eff_dict_errors[f'{eff_key}_error'] = []
    eff_dict[eff_key].append(round(eff_channel,5))
    eff_dict_errors[f'{eff_key}_error'].append((error_low + error_up) / 2)

def makeplot(eff_dict,eff_dict_errors, labels, linestyles, channel, x_values, sample, deepTauWP, deepTauVersion, ratio_ref,ref_errors, suffix='', colors=[], wantNEvents=False, wantLegend=False, verbose=False):
    if colors==[]:
        colors = ['blue', 'green', 'red', 'orangered', 'purple', 'pink', 'yellow', 'cyan','black','brown','lime','navy','crimson']
    markers = ['o', '^', 's', 'D', 'x', 'v', 'p', '*','o']
    plt.figure(figsize=(20,15))
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    for i, (key, values) in enumerate(eff_dict.items()):
        label = None
        if i < len(labels):
            label = labels[i]
        ratio_values = [(val/ref) for val,ref in zip(values,ratio_ref)]
        errors = eff_dict_errors[f'{key}_error']
        ratio_errors = []
        for val,error,ref,err_ref,ratio in zip(values,errors,ratio_ref,ref_errors,ratio_values):
            ratio_errors.append(ratio*(math.sqrt(math.pow(error/val, 2) + math.pow(err_ref/ref,2))))
        if verbose:
            print(key)
            for mass,ratio,value,err in zip(x_values,ratio_values,values, errors):
                print(f"mass = {mass} \t value = {value} \t error = {err} \t ratio = {ratio}")
        ax1.errorbar(x_values, values, yerr=eff_dict_errors[f'{key}_error'],color=colors[i % len(colors)],marker=markers[i % len(markers)],markersize=2,linestyle=linestyles[i % len(markers)], label=label)
        ax2.errorbar(x_values, ratio_values, yerr=ratio_errors,color=colors[i % len(colors)],marker=markers[i % len(markers)],markersize=2, linestyle=linestyles[i % len(markers)])#,label=label)
        ax1.yaxis.get_label().set_fontsize(17)
        ax2.yaxis.get_label().set_fontsize(17)
        #ax1.grid(linestyle='-', linewidth='0.5')
        ax2.xaxis.get_label().set_fontsize(17)

    #### Legend + titles + axis adjustment ####
    #plt.title(f'efficiencies for channel {channel} and {sample}')
    plt.xlabel(r'$m_X$ (GeV/$c^2$)')
    ax1.set_ylabel('efficiency')
    if wantNEvents:
        ax1.set_ylabel('nEvents')
    ax2.set_ylabel('path/legacy')
    if not wantNEvents:
        ax1.set_ylim(0., 1.05)
    #ax2.set_yscale('log')
    ax2.set_xscale('log')
    xticks = [250, 500, 800, 1000, 1500, 2000, 2500 ,3000]  # Esempio di numeri da visualizzare sull'asse logaritmico
    ax2.set_xticks(xticks)
    ax2.set_xticklabels([str(num) for num in xticks])

    figName = f'Studies/Triggers/eff_{channel}_{sample}_{deepTauWP}_{deepTauVersion}{suffix}'
    #fig.legend(bbox_to_anchor=(1.3, 0.6))
    #ax1.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), shadow=True, ncol=1)
    if wantLegend:
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), shadow=True, ncol=1)
    plt.tight_layout()
    if wantNEvents:
        figName+='_nEvents'
    if wantLegend:
        figName += '_WithLeg'
    plt.savefig(f'{figName}.pdf', format='pdf')
