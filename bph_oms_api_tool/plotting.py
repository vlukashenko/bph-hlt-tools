def reformat_for_plotting(data, keys_to_remove = 'run_id', keys_to_merge = ['lumisection', 'prescale_indices', 'recorded_lumi_per_lumisection' ]):
    reformatted_dict = {}
    for k in keys_to_merge:
        reformatted_dict[k] = []
        
    assert isinstance(data, list), print("data should be a valid list and not ", type(data))
    for it in data:
        assert isinstance(it, dict),  "{} should be a valid python dictionary and not {}".format(it, type(it))
        assert keys_to_remove in it.keys(), "{} is not in the keys of the dictionary {}".format(keys_to_remove, it)
        assert 'prescale_idx' in it.keys(), "prescale_idx is not in the keys of the dictionary {}".format(it)
        assert isinstance(it['prescale_idx'], dict), "{} is not a valid python dict, but {}".format(it['prescale_idx'], type(it['prescale_idx']))
    
        for k in keys_to_merge:
             assert k in it['prescale_idx'].keys(), "{} is not in dictionary {} keys: {}".format(k, it['prescale_idx'], it['prescale_idx'].keys())
             reformatted_dict[k].extend(it['prescale_idx'][k])
    return reformatted_dict



