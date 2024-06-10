import omsapi
from messages import error, warning, info
from generic_helpers import removeDuplicates

def getPrescaleIdx_per_run(omsapi, run_num, attributes = ['recorded_lumi_per_lumisection', 'prescale_name', 'prescale_index']):
    query = omsapi.query("lumisections")
    query.filter("run_number", run_num)
    query.set_verbose(True)
    query.per_page = 10000  # to get all names in one go
    query.sort("lumisection_number", asc=True)
    query.attrs(attributes)
    internal_data = query.data().json()
    current_prescale_index = None
    scan_per_run = {'lumisection':[], 'recorded_lumi_per_lumisection': [],  'prescale_indices': []}
    if internal_data['data'] == []:
        return None
    for d in internal_data['data']:
        #print(run_num, d['id'], d['attributes']['recorded_lumi_per_lumisection'], d['attributes']['prescale_name'], d['attributes']['prescale_index'])
        if current_prescale_index != d['attributes']['prescale_index'] and d['attributes']['prescale_index'] != None:
            current_prescale_index = d['attributes']['prescale_index'] 
            scan_per_run['lumisection'].append(d['id'])
            scan_per_run['recorded_lumi_per_lumisection'].append(d['attributes']['recorded_lumi_per_lumisection'])
            scan_per_run['prescale_indices'].append(current_prescale_index)
    return scan_per_run


def getRunInfo(omsapi, run_num, run_type = "collisions2024", remove_if = "None", attributes = ['l1_hlt_mode', 'components', 'recorded_lumi', 'initial_prescale_index']):
    query = omsapi.query("runs")
    query.filter("run_number", run_num)
    query.set_verbose(True)
    query.per_page = 10000  # to get all names in one go
    query.attrs(attributes)
    internal_data = query.data().json()
    if internal_data['data'] == []:
        return None, None, None
    for d in internal_data['data']:
        if d['attributes']['l1_hlt_mode'] != run_type:
            print('Run {}  is not of type {}, but of {}'.format(run_num, run_type, d['attributes']['l1_hlt_mode']))
            return None, None, None
        else:
            if remove_if == "missing_components"  and len(d['attributes']['components']) != len(baseline_components):
                    missing_components = [v for v in baseline_components if v not in d['attributes']['components']]
                    print('Following components are missing from run {} : {}'.format(run_num, missing_components))
                    return None, None, None
            return d['attributes']['recorded_lumi'], d['attributes']['initial_prescale_index'], d['attributes']['components'] 


def getLastLumiSec(omsapi, run_num):
    query = omsapi.query("lumisections")
    query.filter("run_number", run_num)
    #query.set_verbose(True)
    query.per_page = 10000  # to get all names in one go
    last_lumisection = query.data().json()['data'][-1]
    return last_lumisection


def getL1RateInfo_helper(harvested_data, base_path, triggers):
    for entry in base_path:
        n =  entry['attributes']["name"]
        if n not in triggers:  
            info("Tigger {} is not found.".format(n))
            continue
        if entry['attributes']["last_lumisection_number"] != entry['attributes']["first_lumisection_number"] : print("ATTENTION")
        
        if n not in harvested_data.keys():
            harvested_data[n] = {
                    "lumisec" : [entry['attributes']["first_lumisection_number"]], 
                    "pre_rate" :  [entry["attributes"]["pre_dt_before_prescale_rate"]], 
                    "post_rate": [entry["attributes"]["post_dt_rate"]], 
                    "bit": [entry["attributes"]["bit"]], 
                    "initial_prescale": [entry["attributes"]["initial_prescale"]], 
                    "final_prescale": [entry["attributes"]["final_prescale"]]
                    }
        else:
            harvested_data[n]["lumisec"].append(entry['attributes']["first_lumisection_number"])
            harvested_data[n]["pre_rate"].append(entry["attributes"]["pre_dt_before_prescale_rate"])
            harvested_data[n]["post_rate"].append(entry["attributes"]["post_dt_rate"])
            harvested_data[n]["bit"].append(entry["attributes"]["bit"])
            harvested_data[n]["initial_prescale"].append(entry["attributes"]["initial_prescale"])
            harvested_data[n]["final_prescale"].append(entry["attributes"]["final_prescale"])
  
    last_lumisec = entry['attributes']["last_lumisection_number"]
    return harvested_data, last_lumisec


def collectTriggerBits(omsapi, run_num, triggers):
    query = omsapi.query("l1algorithmtriggers")
    query.filter("run_number", run_num)
    #query.set_verbose(True)
    query.per_page = 10000  # to get all names in one go
    base_path = query.data().json()['data']
    if len(base_path) == 0:
        warning("The l1algorithmtriggers endpoint returns empty fileds for Run {}. Therefore I skip.".format(run_num))
        return []
 
    bits = []
    for entry in base_path:
        n =  entry['attributes']["name"]
        if n not in triggers:  continue
        if entry['attributes']["last_lumisection_number"] != entry['attributes']["first_lumisection_number"] : print("ATTENTION")
        bits.append(entry['attributes']['bit'])
 
    bits = removeDuplicates(bits)
    return bits

def getL1RateInfo(omsapi, run_num, triggers, attributes = ["initial_prescale", "name", "bit", "pre_dt_before_prescale_rate", "final_prescale", "post_dt_rate", "last_lumisection_number", "first_lumisection_number"]):
    last_lumisection = getLastLumiSec(omsapi, run_num)
    trigger_bits = collectTriggerBits(omsapi, run_num, triggers)
    query = omsapi.query("l1algorithmtriggers")
    if len(trigger_bits) == 0:
        warning("No requested triggers are found in Run {}. Therefore I skip.".format(run_num))
        return {}
    harvested_data = {}
    for b in trigger_bits:
        query.clear_filter().filter("run_number", run_num).filter("bit", b)
        query.set_verbose(True)
        query.per_page = 10000  # to get all names in one go
        query.attrs(attributes)
        base_path = query.data().json()['data']
 
        if len(base_path) == 0:
             warning("The l1algorithmtriggers endpoint returns empty fileds for Run {}. Therefore I skip.".format(run_num))
             return {}
    
        single_harvested_data, single_lumisec = getL1RateInfo_helper(harvested_data, base_path, triggers)
        harvested_data.update(single_harvested_data)
    return harvested_data


def getL1TotalRateInfo(omsapi, run_num, attributes = ['l1a_physics', 'l1a_random', 'l1a_calibration', 'l1a_total']):
    query = omsapi.query("l1triggerrates")
    query.filter("run_number", run_num)
    query.set_verbose(True)
    query.per_page = 10000  # to get all names in one go
    #print("INSIDE: ", query.data().json()['data'])
    query.attrs(attributes)
    base_path = query.data().json()['data'][0]['attributes']
    return base_path['l1a_total']['rate'],  base_path['l1a_physics']['rate'],  base_path['l1a_calibration']['rate'], base_path['l1a_random']['rate']

def getLumiSecLumi(omsapi, lumisec_num, attributes = ["recorded_lumi_per_lumisection"]):
    query = omsapi.query("lumisections")
    query.filter("run_number", lumisec_num.split("_")[0])
    #query.set_verbose(True)
    query.per_page = 10000  # to get all names in one go
 
    query.filter("lumisection_number", lumisec_num.split("_")[1])
    query.attrs(attributes)
    return query.data().json()['data'][0]['attributes']["recorded_lumi_per_lumisection"]


def getLatestRun(omsapi):
    query = omsapi.query("runs")
    query.set_verbose(True)
    query.sort("run_number", asc=False)
    return query.data().json()['data'][0]['id']

def getPreviousRun(omsapi, r):
    query = omsapi.query("runs")
    query.set_verbose(True)
    query.sort("run_number", asc=False)
    query.filter("run_number", r, "LE")
    info = query.data().json()['data']
    assert r == info[0]['id'], error("This is a sanity check. Run {} is not {}".format(r, info['id']))
    return info[1]['id'] 

def getNextRun(omsapi, r):
    query = omsapi.query("runs")
    query.set_verbose(True)
    query.sort("run_number")
    query.filter("run_number", r, "GE")
    info = query.data().json()['data']
    assert r == info[0]['id'], error("This is a sanity check. Run {} is not {}".format(r, info['id']))
    if len(info) == 1:
        error("Run {} is the latest".format(r))
    return info[1]['id'] 


def getHLTRateInfo(omsapi, trig, run_num, attributes = []):
    query = omsapi.query("hltpathrates")
    query.filter("run_number", run_num)
    print(trig)
    query.filter("path_name", trig)
    query.set_verbose(True)
    query.per_page = 10000  # to get all names in one go
    print("HLT rates" , query.data().json())
    return 0,0,0,0
    #query.attrs(attributes)
    #base_path = query.data().json()['data'][0]['attributes']



def getAllFills(omsapi, era):
    query = omsapi.query("eras")


    query.set_verbose(True)
    query.per_page = 10000
    for d in query.data().json()['data']:
        if d["id"] != era:
            continue
        
        first_fill = d['attributes']['start_fill']
        last_fill = d['attributes']['end_fill'] 
    
    query = omsapi.query('fills')
    query.filter("fill_number", first_fill, operator="GT")
    query.filter("fill_number", last_fill, operator="LT")
    
    query.set_verbose(True)
    query.per_page = 10000
    
    good_fills = {}
    for f in query.data().json()['data']:
        if f['attributes']['stable_beams'] and f['attributes']['bunches_colliding'] != None:
           good_fills[f['id']] = [f['attributes']['first_run_number'], f['attributes']['last_run_number']]

    return good_fills

 
def isGoodRun(omsapi, r):
    query = omsapi.query("runs")
    query.filter('run_number', r)
    query.set_verbose(True)
    query.per_page = 10000
    query.attrs(['components', 'l1_hlt_mode',  'recorded_lumi']) 
    harvest_info = query.data().json()['data']
    flag_info_exist = len(harvest_info) > 0 and harvest_info[0]['attributes']['l1_hlt_mode'] != None and harvest_info[0]['attributes']['recorded_lumi'] != None
    flag = True if flag_info_exist and "collisions" in harvest_info[0]['attributes']['l1_hlt_mode'] and harvest_info[0]['attributes']['recorded_lumi'] > 0 else False
    if flag:
       info("Run {} is good. See more in helpers::isGoodRun. Continue.".format(r))
    else:
       info("Run {} is bad. See more in helpers::isGoodRun. Skipping.".format(r))

    return flag
