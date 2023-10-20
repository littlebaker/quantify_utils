#%%
import quantify_core
import tempfile
from pathlib import Path
import os

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from qcodes import Instrument, ManualParameter, Parameter, validators
from scipy.optimize import minimize_scalar

import quantify_core.data.handling as dh
from quantify_core.analysis import base_analysis as ba
from quantify_core.analysis import cosine_analysis as ca
from quantify_core.measurement import Gettable, Settable, MeasurementControl
from quantify_core.utilities.dataset_examples import mk_2d_dataset_v1
from quantify_core.utilities.examples_support import mk_cosine_instrument
from quantify_core.utilities.inspect_utils import display_source_code

dh.set_datadir(os.path.join("./", "quantify-data"))
meas_ctrl = MeasurementControl("test_meas_ctrl")
# from hunanu.Instrument.VNA_keysightN5231b_4Portsand2Ports import keysight_vna
# a = keysight_vna("xxx")
import hunanu.Instrument.M4iCard 
from hunanu.Instrument.VNA_keysightN5231b_4Portsand2Ports import keysight_vna
a = keysight_vna("dfadfa")

#%%
# initialize source 
mw_source1 = Instrument("ms_source1")

mw_source1.freq = ManualParameter(
    name="freq",
    label="Frequency",
    unit="Hz",
    vals=validators.Numbers(),
    initial_value=1.0,
)

mw_source1.power = ManualParameter(
    name="power",
    label="Power",
    unit="V",
    vals=validators.Numbers(-5, 5),
    initial_value=0
)

pulsar_QRM = Instrument("pulsar_QRM")
# NB: for brevity only, this not the proper way of adding parameters to QCoDeS instruments
pulsar_QRM.signal = Parameter(
    name="sig_a", label="Signal", unit="V", get_cmd=lambda: mw_source1.power() * mw_source1.freq() * 1e-8 + np.random.random()*5
)


meas_ctrl.settables(mw_source1.freq)
meas_ctrl.setpoints_grid([np.linspace(5e9, 5.2e9, 100), np.linspace(-3, 3, 100)])
meas_ctrl.gettables(pulsar_QRM.signal)
dset = meas_ctrl.run(name="Frequency sweep")
dset_grid = dh.to_gridded_dataset(dset)





# with tempfile.TemporaryDirectory() as tempdir:
#     # old_dir = dh.get_datadir()
#     # dh.set_datadir(os.path.join(tempdir, "quantify-data"))
#     print(dh.get_datadir())
#     (Path(dh.get_datadir()) / "20210301").mkdir()
#     (Path(dh.get_datadir()) / "20210428").mkdir()
#     qset = mk_2d_dataset_v1()
#     ba.BaseAnalysis(dataset=qset).run()
#     # dh.set_datadir(old_dir)