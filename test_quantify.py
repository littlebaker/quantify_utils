#%%
from time import sleep
import quantify_core
import tempfile
from pathlib import Path
import os

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from qcodes import Instrument, ManualParameter, Parameter, validators
from qcodes.instrument import find_or_create_instrument
from scipy.optimize import minimize_scalar

import quantify_core.data.handling as dh
from quantify_core.analysis import base_analysis as ba
from quantify_core.analysis import cosine_analysis as ca
from quantify_core.measurement import Gettable, Settable, MeasurementControl
from quantify_core.utilities.dataset_examples import mk_2d_dataset_v1
from quantify_core.utilities.examples_support import mk_cosine_instrument
from quantify_core.utilities.inspect_utils import display_source_code
from quantify_core.visualization.pyqt_plotmon import PlotMonitor_pyqt
from quantify_core.visualization.pyqt_plotmon_remote import RemotePlotmon

from quantify_utils.measurement_utils import BatchedGettable, BatchedParam, MP
from quantify_utils.web_plotmon2 import PlotMonitor_web


dh.set_datadir(os.path.join("./", "quantify-data"))
meas_ctrl = MeasurementControl("meas_ctrl")
# meas_ctrl.set_experiment_data(MP("RT_ATT", 60), False)
plotmon = PlotMonitor_web("plotmon")
meas_ctrl.instr_plotmon("plotmon")




#%%
# initialize source 
mw_source1 = find_or_create_instrument(Instrument, "ms_source1")

mw_source1.freq = ManualParameter(
    name="freq",
    label="Frequency",
    unit="Hz",
    # vals=validators.Numbers(),
    initial_value=1e9,
)

mw_source1.power = ManualParameter(
    name="power",
    label="Power",
    unit="V",
    # vals=validators.Numbers(-5, 5),
    initial_value=1
)
p = ManualParameter("p", initial_value=0)
def _get_sig():
    sleep(0.3)
    x = mw_source1.power() * mw_source1.freq() * 1e-8 + p.get()
    return x

pulsar_QRM = find_or_create_instrument(Instrument, "pulsar_QRM")
# NB: for brevity only, this not the proper way of adding parameters to QCoDeS instruments
pulsar_QRM.signal = Parameter(
    name="sig_a", label="Signal", unit="V", get_cmd=_get_sig
)


meas_ctrl.settables(p)
meas_ctrl.setpoints(np.linspace(0, 10, 100)[::-1])
meas_ctrl.gettables(pulsar_QRM.signal)
dset = meas_ctrl.run(name="Frequency sweep")

dset_grid = dh.to_gridded_dataset(dset)

import hvplot.xarray





