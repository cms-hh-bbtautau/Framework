import yaml
class Triggers():
    dict_legtypes = {"Electron":"Leg::e", "Muon":"Leg::mu", "Tau":"Leg::tau"}

    def __init__(self, triggerFile, deltaR_matching=0.4):
        with open(triggerFile, "r") as stream:
            self.trigger_dict= yaml.safe_load(stream)
        self.deltaR_matching = deltaR_matching

    def ApplyTriggers(self, df, isData = False):
        hltBranches = []
        matchedObjectsBranches= []
        for path, path_dict in self.trigger_dict.items():
            path_key = 'path'
            if 'path' not in path_dict:
                path_key += '_data' if isData else '_MC'
            keys = [k for k in path_dict[path_key]]
            # check that HLT path exists:
            for key in keys:
                trigger_string = key.split(' ')
                trigName = trigger_string[0]
                if trigName not in df.GetColumnNames():
                    print(f"{trigName} does not exist!!")
                    path_dict[path_key].remove(key)
            or_paths = " || ".join(f'({p})' for p in path_dict[path_key])
            or_paths = f' ( { or_paths } ) '
            additional_conditions = [""]
            total_objects_matched = []
            for leg_id, leg_tuple in enumerate(path_dict['legs']):
                leg_dict_offline= leg_tuple["offline_obj"]
                type_name_offline = leg_dict_offline["type"]
                var_name_offline = f'{type_name_offline}_offlineCut_{leg_id+1}_{path}'
                df = df.Define(var_name_offline, leg_dict_offline["cut"])
                if not leg_tuple["doMatching"]:
                    if not leg_dict_offline["type"].startswith('MET'):
                        var_name_offline = f'{leg_dict_offline["type"]}_idx[{var_name_offline}].size()>0'
                    additional_conditions.append(var_name_offline)
                else:
                    df = df.Define(f'{type_name_offline}_{var_name_offline}_sel', f'httCand.isLeg({type_name_offline}_idx, {self.dict_legtypes[type_name_offline]}) && ({var_name_offline})')
                    leg_dict_online= leg_tuple["online_obj"]
                    var_name_online =  f'{leg_dict_offline["type"]}_onlineCut_{leg_id+1}_{path}'
                    df = df.Define(f'{var_name_online}',f'{leg_dict_online["cut"]}')
                    matching_var = f'{leg_dict_offline["type"]}_Matching_{leg_id+1}_{path}'

                    df = df.Define(matching_var, f"""FindMatchingSet( {type_name_offline}_{var_name_offline}_sel,{var_name_online},
                                                                            {leg_dict_offline["type"]}_p4, TrigObj_p4,{self.deltaR_matching} )""")
                    total_objects_matched.append(f'{{ {self.dict_legtypes[type_name_offline]}, {matching_var} }}')

            legVector = f'{{ { ", ".join(total_objects_matched)} }}'
            df = df.Define(f'hasOOMatching_{path}_details', f'HasOOMatching({legVector} )')
            df = df.Define(f'hasOOMatching_{path}', f'hasOOMatching_{path}_details.first')
            for offline_leg_id in range(2):
                matching_var_bool = f'tau{offline_leg_id+1}_HasMatching_{path}'
                df = df.Define(matching_var_bool,
                               f'hasOOMatching_{path}_details.second.count(LegIndexPair(httCand.leg_type.at({offline_leg_id}), httCand.leg_index.at({offline_leg_id}) ) ) > 0')
                matchedObjectsBranches.append(matching_var_bool)
            fullPathSelection = f'{or_paths} &&  hasOOMatching_{path}'
            fullPathSelection += ' && '.join(additional_conditions)
            hltBranch = f'HLT_{path}'
            hltBranches.append(hltBranch)
            df = df.Define(hltBranch, fullPathSelection)
        total_or_string = ' || '.join(hltBranches)
        df = df.Filter(total_or_string)
        hltBranches.extend(matchedObjectsBranches)
        return df,hltBranches
