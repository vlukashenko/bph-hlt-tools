import json

ALL_BPH_Datasets = ['ParkingDoubleMuonLowMass0',
                    "ParkingDoubleMuonLowMass1",
                    "ParkingDoubleMuonLowMass2",
                    "ParkingDoubleMuonLowMass3",
                    "ParkingDoubleMuonLowMass4",
                    "ParkingDoubleMuonLowMass5",
                    "ParkingDoubleMuonLowMass6",
                    "ParkingDoubleMuonLowMass7",
                    "ParkingSingleMuon0",
                    "ParkingSingleMuon1",
                    "ParkingSingleMuon2",
                    "ParkingSingleMuon3",
                    "ParkingSingleMuon4",
                    "ParkingSingleMuon5",
                    "ParkingSingleMuon6",
                    "ParkingSingleMuon7",]
NOT_BPH_Datasets = ['AlCaLowPtJet',
                    'AlCaLumiPixelsCountsExpress',
                    'AlCaLumiPixelsCountsPrompt',
                    'AlCaP0', 'AlCaPPSExpress',
                    'AlCaPPSPrompt', 'AlCaPhiSym',
                    'BTagMu', 'Commissioning',
                    'Cosmics', 'DQMGPUvsCPU',
                    'DQMOnlineBeamspot', 'DQMPPSRandom',
                    'DisplacedJet', 'EGamma0', 'EGamma1','EcalLaser', 
                    'EphemeralHLTPhysics0',"EphemeralHLTPhysics7",
                    "EphemeralHLTPhysics1","EphemeralHLTPhysics2",
                    "EphemeralHLTPhysics3","EphemeralHLTPhysics4",
                    "EphemeralHLTPhysics5","EphemeralHLTPhysics6",                    
                    "EphemeralZeroBias0","EphemeralZeroBias1",
                    "EphemeralZeroBias2","EphemeralZeroBias3",
                    "EphemeralZeroBias4","EphemeralZeroBias5",
                    "EphemeralZeroBias6","EphemeralZeroBias7",
                    "EventDisplay","ExpressAlignment","ExpressPhysics",
                    "HLTMonitor","HLTPhysics","HcalNZS","JetMET0","JetMET1",
                    "L1Accept","Muon0","Muon1","MuonEG","NoBPTX","OnlineMonitor"
                    ]


def create_json(dictionary, fname):
    with open(fname, 'w+') as file:
        json.dump(dictionary, file, indent=2)

def get_all_datasets(Menu):
    datasets_obj = Menu.process.datasets
    datasets_dict= datasets_obj.parameters_()
    return datasets_dict

def get_BPH_Datasets(Menu):
    global ALL_BPH_Datasets, NOT_BPH_Datasets
    BPH_datasets = dict()
    all_paths = list()
    datasets_dict = get_all_datasets(Menu)
    for dataset, paths in datasets_dict.items():
        print(dataset)
        if dataset in ALL_BPH_Datasets: 
            open_='y'
        elif dataset in NOT_BPH_Datasets: 
            open_='n'
        else: 
            open_ = input('is BPH (y/n): ')

        if open_.lower() == 'y':
            BPH_datasets[dataset] = paths.value()
            all_paths += paths.value()
            ALL_BPH_Datasets.append(dataset)
        else:
            NOT_BPH_Datasets.append(dataset)
    
    NOT_BPH_Datasets = list(set(NOT_BPH_Datasets))
    all_paths = list(set(all_paths))
    return BPH_datasets, all_paths

def get_hlt_obj(Menu, HLTpath):
    process = Menu.process
    path = getattr(process, HLTpath)
    return path

def get_dependencies(path):
    return path.directDependencies()

def get_path_components(Menu, HLTpath):
    obj =  get_hlt_obj(Menu, HLTpath)
    return get_dependencies(obj)

def unpack_value(value):
    
    try:
        unpack1 = value.value()
    except AttributeError:
        return value

    if 'PSet' in str(type(unpack1)):
        unpack1 =  unpack1.parameters_()
        if type(unpack1) == dict:
            for key, val in unpack1.items():
                unpack1[key]  = unpack_value(val)

    elif type(unpack1)==list:

        for i,entry in enumerate(unpack1):

            if 'PSet' in str(type(entry)):
                unpack1[i]=entry.parameters_()

            if type(unpack1[i]) == dict:
                for key_e, val_e in unpack1[i].items():
                    unpack1[i][key_e]  = unpack_value(val_e)

    return unpack1

def get_all_filters(Menu, path):
    process = Menu.process
    paths_filters = dict()
    # Get the object
    path_obj =  get_hlt_obj(Menu, path)
    # Iterate over all modules
    for module in path_obj.moduleNames():
        module_obj = getattr(process, module)
        #Only save filters
        if 'filter' in str(type(module_obj)).lower(): 
            paths_filters[module] = dict()
            dict_ = module_obj.parameters_()
            for key, val in dict_.items():
                paths_filters[module][key]  = unpack_value(val)
            paths_filters[module]['type_'] = module_obj.type_()\
    
    return paths_filters




if __name__=="__main__":

    import importlib

    import argparse
    parser = argparse.ArgumentParser(description="Get all paths of interest for selected datasets")
    parser.add_argument('--menu', type=str, help='Menu input')
    args = parser.parse_args()

    menu_str = args.menu.replace('.py','').replace('/','.')
    print(menu_str)
    
    menu     = importlib.import_module(menu_str)
    datasets = get_all_datasets(menu)
    
    BPH_datasets  = dict()
    paths_defs    = dict()
    paths_filters = dict()
    all_paths     = list()


    for dataset, paths in datasets.items():
        print(dataset, end=' ')
        if dataset in ALL_BPH_Datasets: 
            print(' -> Included')
            open_='y'
        elif dataset in NOT_BPH_Datasets: 
            print(' -> XXX Ignored')
            open_='n'
        else: open_ = input('is BPH (y/n): ')
        if open_.lower() == 'y':
            BPH_datasets[dataset] = paths.value()
            all_paths += paths.value()
        else:
            NOT_BPH_Datasets.append(dataset)
            NOT_BPH_Datasets = list(set(NOT_BPH_Datasets))

    # Remove duplicates
    all_paths = list(set(all_paths))

    # Get all components for each path
    for path in all_paths:
        paths_defs[path] = get_path_components(menu, path)

    # Get all filters and unpack the variables
    process = menu.process
    for path in all_paths:
        paths_filters[path] = dict()
        # Get the object
        path_obj =  get_hlt_obj(menu, path)
        # Iterate over all modules
        for module in path_obj.moduleNames():
            module_obj = getattr(process, module)
            #Only save filters
            if 'filter' in str(type(module_obj)).lower(): 
                paths_filters[path][module] = dict()
                dict_ = module_obj.parameters_()
                for key, val in dict_.items():
                    paths_filters[path][module][key]  = unpack_value(val)
                paths_filters[path][module]['type_'] = module_obj.type_()

    os.makedirs('Paths', exist_ok=True)
    create_json(paths_filters, f'Paths/Path_filters_{menu_str.split(".")[1]}.json')



# for name, menu in zip(Names, Menus):
#     print('\n\n\n\n\n')
#     print(name, menu)
#     datasets, paths = get_BPH_Datasets(menu)
#     create_json(datasets, f'Datasets_{name}.json')
#     filters = dict()
#     for path in paths:
#         filters[path] = get_all_filters(menu, path)
#     create_json(filters, f'Path_filters_{name}.json')











# datasets_menu = get_all_datasets(menu)

# BPH_datasets  = dict()
# paths_defs    = dict()
# paths_filters = dict()
# all_paths     = list()

# # Get all BPH datasets
# for dataset, paths in datasets_menu.items():
#     print(dataset)
#     if dataset in ALL_BPH_Datasets: open_='y'
#     elif dataset in NOT_BPH_Datasets: open_='n'
#     else: open_ = input('is BPH (y/n): ')
#     if open_.lower() == 'y':
#         BPH_datasets[dataset] = paths.value()
#         all_paths += paths.value()
#     else:
#         NOT_BPH_Datasets.append(dataset)
#         NOT_BPH_Datasets = list(set(NOT_BPH_Datasets))

# # Remove duplicates
# all_paths = list(set(all_paths))

# # Get all components for each path
# for path in all_paths:
#     paths_defs[path] = get_path_components(menu, path)

# # Get all filters and unpack the variables
# process = menu.process
# for path in all_paths:
#     paths_filters[path] = dict()
#     # Get the object
#     path_obj =  get_hlt_obj(menu, path)
#     # Iterate over all modules
#     for module in path_obj.moduleNames():
#         module_obj = getattr(process, module)
#         #Only save filters
#         if 'filter' in str(type(module_obj)).lower(): 
#             paths_filters[path][module] = dict()
#             dict_ = module_obj.parameters_()
#             for key, val in dict_.items():
#                 paths_filters[path][module][key]  = unpack_value(val)
#             paths_filters[path][module]['type_'] = module_obj.type_()


#create_json(paths_filters, 'Path_filters.json')
