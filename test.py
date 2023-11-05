from qcodes.instrument import find_or_create_instrument
from quantify_utils.GS200 import YokogawaGS
from quantify_utils.KeysightVNA import KeysightVNA
   
   
# vna_name = "vna2"
# vna_trace = "all"
# vna = find_or_create_instrument(
#     KeysightVNA, name=vna_name, device_name=vna_name, trace=vna_trace
# )


gs_name = "DC6"
gs = find_or_create_instrument(YokogawaGS, gs_name, gs_name)