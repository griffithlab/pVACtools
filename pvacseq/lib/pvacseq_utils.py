def prediction_method_to_iedb_lookup_dict():
    return{
        'NetMHCpan' : 'netmhcpan',
        'NetMHC'    : 'ann',
        'SMMPMBEC'  : 'smmpmbec',
        'SMM'       : 'smm',
        'NetMHCcons': 'netmhccons',
        'PickPocket': 'pickpocket',
    }

def prediction_methods():
    prediction_method_lookup_dict = prediction_method_to_iedb_lookup_dict()
    return sorted(prediction_method_lookup_dict.keys())

def iedb_to_prediction_method_lookup_dict():
    prediction_method_lookup_dict = prediction_method_to_iedb_lookup_dict()
    return {v: k for k, v in prediction_method_lookup_dict.items()}

def iedb_prediction_methods():
    prediction_method_lookup_dict = iedb_to_prediction_method_lookup_dict()
    return sorted(prediction_method_lookup_dict.keys())
