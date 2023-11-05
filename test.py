from qcodes.instrument import find_or_create_instrument
from quantify_utils.GS200 import YokogawaGS, gs_param_set
from quantify_utils.KeysightVNA import KeysightVNA
from quantify_utils.web_monitor2 import main
import panel as pn
from quantify_utils.visualization import *

# vna_name = "vna2"
# vna_trace = "all"
# vna = find_or_create_instrument(
#     KeysightVNA, name=vna_name, device_name=vna_name, trace=vna_trace
# )


# gs_name = "DC6"
# gs = find_or_create_instrument(YokogawaGS, gs_name, gs_name)
# gs_param_set(gs, 1e-4, "CURR", "v2")

# gs.level(-1)
if __name__ == "__main__":
    elec_delay = 0
    def _draw_undone(x):
        return  pn.Column(spec_draw_amp(x), spec_draw_phase(x, elec_delay))
    def _draw_done(x):
        return  pn.Column(spec_draw_widget_amp(x), spec_draw_widget_phase(x, elec_delay))
    
    main(
        "quantify-data/Gmon04",
        "20231105-190747-903-065e8e-qubit spectrum gmon from waveguide",
        _draw_undone,
        _draw_done,
        port=5010,
    )
