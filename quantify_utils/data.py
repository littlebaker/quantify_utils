from quantify_core.data.experiment import QuantifyExperiment
from quantify_core.data.types import TUID
import os
import re



def get_experiments(exp_dir: str):
    exp_path = os.path.abspath(exp_dir)
    dirs = os.listdir(exp_path)
    
    tuid_str = []
    for dirname in dirs:
        m = re.match(r"([\d]{8}-[\d]{6}-[\d]{3}-.{6}).*", dirname)
        if m is not None:
            tuid_str.append(m.group(1))
        
    
    tuid_list = map(
        lambda x: TUID(x),
        tuid_str
    )
    
    epmt_list = []
    for tuid in tuid_list:
        try:
            epmt = QuantifyExperiment(tuid)
        except:
            continue
        epmt.load_dataset()
        epmt_list.append(epmt)
        
    return epmt_list


