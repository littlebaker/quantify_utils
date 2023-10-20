from matplotlib import colorbar
from matplotlib.pyplot import set_cmap
import numpy as np

import quantify_core.data.handling as dh
from quantify_core.measurement import MeasurementControl
from quantify_core.data.experiment import QuantifyExperiment

from qcodes.instrument_drivers.Keysight import KeysightN5232B
from qcodes.parameters import ManualParameter, Parameter


dh.set_datadir("./quantify-data/")

meas_ctrl = MeasurementControl("Scan_Cavity_Frequency")


vna = KeysightN5232B("VNA2", "TCPIP0::172.25.146.89")
points = 201

vna.power(-15)
vna.points(points)

vna.if_bandwidth(100)

vna.electrical_delay(0)
vna.sweep_mode("SING")

vna.trigger_source("IMM")

vna.trace("S21")

class batched_param:
    def __init__(self, param, batch_points) -> None:
        self.name = param.name
        self.label = param.label
        self.unit = param.unit
        self.param = param
        self.batched = True
        self.batch_size = batch_points
        
    def get(self):
        return self.param.get()
    
    def set(self, x):
        self.param.set(x)
        
def set_freq_range(fr):
    start = fr[0]
    stop = fr[-1]
    points = len(fr)
    
    vna.start(start)
    vna.stop(stop)
    vna.points(points)
        
freq_range = Parameter(
    "freq_range", 
    label="Frequency", 
    unit='V',
    set_cmd=set_freq_range
)

gs_cur = ManualParameter(
    name = "current_value",
    label = "Current",
    unit = "mA"
)

meas_ctrl.settables([batched_param(freq_range, points),gs_cur])
# meas_ctrl.setpoints(np.linspace(6e9, 7e9, points))
meas_ctrl.setpoints_grid([np.linspace(6e9, 7e9, points),np.linspace(1,5,5)])
meas_ctrl.gettables([batched_param(vna.magnitude, points), batched_param(vna.frequency_axis, points)])
dset = meas_ctrl.run("freq_scan_test")

import hvplot.pandas

df = dset.to_dataframe()
df.hvplot(
    colorbar = True,
    kind = "heatmap",
    x = "x0",
    y = "y0",
)




# %%
