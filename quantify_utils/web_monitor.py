
import time
import panel as pn
import quantify_core.data.handling as dh
from quantify_core.data.handling import load_dataset, load_snapshot, get_tuids_containing, to_gridded_dataset
import hvplot.xarray  # noqa: F401
import hvplot.pandas
import numpy as np
import xarray
import holoviews as hv
import sys
import quantify_core.data.handling as dh
from quantify_core.data.handling import load_dataset
from quantify_core.measurement.control import _DATASET_LOCKS_DIR
from quantify_core.visualization.pyqt_plotmon_remote import _safe_load_dataset

hvplot.extension("matplotlib")

pn.extension(template="fast")

if len(sys.argv) > 1:
    print(sys.argv)
    tuid = sys.argv[1]
    dh.set_datadir(sys.argv[2])
else:    
    tuid = "20231028-225013-687-88b211"
    dh.set_datadir("/Users/lbaker/Projects/PycharmProjects/qcodes_test/quantify-data")



def create_json_pane():
    if tuid == "":
        return None
    try:
        snapshot = load_snapshot(tuid)
    except FileNotFoundError:
        return None
    json_pane = pn.pane.JSON(snapshot)
    return json_pane

def create_data_plot(interval, run):
    if not run:
        yield "Calculation did not run yet"
        return
    if tuid == "":
        yield None
    while True:
        time.sleep(interval)
        print("update")
    
        dset = _safe_load_dataset(tuid, _DATASET_LOCKS_DIR)
        print(np.any(dset.y0.isnull().to_numpy()))
        if not np.any(dset.y0.isnull().to_numpy()):
            exit(1)

        if dset.attrs["grid_2d"]:
            dset_grid = to_gridded_dataset(dset)
            
            yield dset_grid.hvplot.image(x="x0", y="x1", z="y0", dynamic=True)
            # return
        else:
            yield pn.pane.HoloViews(dset.hvplot.line(x="x0", y="y0"))



# def main():
# template = pn.template.BootstrapTemplate(title='Bootstrap Template')
run = pn.widgets.Button(name="Press to run calculation", align='center').servable(target="sidebar")
update_interval = pn.widgets.FloatInput(
    name="Update Interval",
    value = 1,
    start = 0.05
).servable(target="sidebar")
# template.sidebar.append(run)
# template.sidebar.append(update_interval)

plot = pn.Column(pn.bind(create_data_plot, update_interval, run)).servable(target="main")
# template.main.append(plot)

# template.servable()

# main()



