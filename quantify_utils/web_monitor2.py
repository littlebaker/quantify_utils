import sys
import panel as pn
import quantify_core.data.handling as dh
from quantify_core.data.handling import load_snapshot, to_gridded_dataset
import hvplot.xarray  # noqa: F401
import numpy as np
import xarray
import holoviews as hv  # noqa: F401


from quantify_core.measurement.control import _DATASET_LOCKS_DIR
from quantify_core.visualization.pyqt_plotmon_remote import _safe_load_dataset, _last_modified
from .visualization import heatmap_crange_cmap


def create_json_pane(tuid, depth=1, preprocess=None):
    if tuid == "":
        return None
    try:
        snapshot = load_snapshot(tuid)
        if preprocess is not None:
            snapshot = preprocess(snapshot)
    except FileNotFoundError:
        return None
    json_pane = pn.pane.JSON(snapshot, depth=depth)

    return json_pane


class DataPlot:
    def __init__(self, tuid, data_var="y0", draw_undone=None, draw_done=None) -> None:
        self.tuid = tuid
        self.data_var = data_var
        self.draw_undone = draw_undone
        self.draw_done = draw_done
        
        self.is_finished = False
        
        self.dset_cache = _safe_load_dataset(self.tuid, _DATASET_LOCKS_DIR)
        self.dset_last_modified = 0
        
        self.periodic_callback = None
        
        self.plot = pn.pane.panel(self.create_data_plot())
        
        self.interval = 1.0
        self.run = True

    def update_dset(self):
        last_modified = _last_modified(self.tuid)
        if last_modified > self.dset_last_modified:
            self.dset_last_modified = last_modified
            self.dset_cache = _safe_load_dataset(self.tuid, _DATASET_LOCKS_DIR)
            return True
        else:
            return False
        
    def add_periodic_callback(self, *events):
        # if not run:
        #     self.periodic_callback.stop()
        #     self.periodic_callback = None
        for event in events:
            if event.obj.name == "Update Interval":
                self.interval = event.new
            elif event.obj.name == "Press to run calculation":
                self.run = event.new
            else:
                return
        
        if self.periodic_callback is not None:
            self.periodic_callback.stop()
        if self.run:
            self.periodic_callback = pn.state.add_periodic_callback(self._callback, int(self.interval*1000))
        else:
            self.periodic_callback = None
        
    def _callback(self):
        # update dataset
        _modified = self.update_dset()
        print("modified: ", _modified)
        if not _modified:
            return
        else:
            _fig = self.create_data_plot() 
            print(hash(_fig))
            self.plot = _fig
            self.plot.param.trigger("object")
            
        

    def create_data_plot(self):
        self.is_finished = not (list(self.dset_cache.data_vars.values())[0].isnull().any().values)
        is_finished = self.is_finished

        dset = self.dset_cache
        

        if dset.attrs["grid_2d"]:  # 2d-data
            dset_grid = to_gridded_dataset(dset)

            if not is_finished:
                fig = None
                if self.draw_undone is not None:
                    fig = self.draw_undone(dset_grid)
                else:
                    fig = self._default_draw_undone(dset_grid)
                return pn.Column(fig)
            else:
                fig = None
                if self.draw_done is not None:
                    fig = self.draw_done(dset_grid)
                else:
                    fig = self._default_draw_done(dset_grid)
                return pn.Column(fig, pn.pane.Str("Finished"))
            # return
        else:  # 1d-data
            if not is_finished:
                if self.draw_undone is not None:
                    fig = self.draw_undone(dset)
                else:
                    fig = pn.pane.HoloViews(
                        dset.hvplot.line(
                            x="x0",
                            y=self.data_var,
                        )
                        * dset.hvplot.scatter(
                            x="x0",
                            y=self.data_var,
                        )
                    )

                return pn.Column(fig)
            else:
                if self.draw_done is not None:
                    fig = self.draw_done(dset)
                else:
                    fig = pn.Column(
                        pn.pane.HoloViews(
                            dset.hvplot.line(
                                x="x0",
                                y=self.data_var,
                            )
                            * dset.hvplot.scatter(
                                x="x0",
                                y=self.data_var,
                            )
                        ),
                        pn.pane.Str("Finished"),
                    )

                    return fig
        return fig

    def _default_draw_undone(self, dset):
        return dset.hvplot.image(
            x="x0", y="x1", z=self.data_var, dynamic=False, cmap="bwr"
        )

    def _default_draw_done(self, dset):
        return heatmap_crange_cmap(dset, self.data_var, cmap="bwr")


def snapshot_prune(snapshot):
    sm = snapshot["instruments"]["meas_ctrl"]["submodules"]
    if "experiment_data" in sm and len(sm["experiment_data"]["parameters"].keys()) != 0:
        return sm["experiment_data"]["parameters"]
    else:
        sm = snapshot["instruments"]["meas_ctrl"]
        if "metadata" in sm:
            return sm["metadata"]
        else:
            return {}


def main(datadir, tuid, draw_undone=None, draw_done=None, snapshot_preprocess=None, port=5006):
    dh.set_datadir(datadir)
    template = pn.template.BootstrapTemplate(title="Web Monitor")
    run = pn.widgets.Checkbox(name="Press to run calculation", align="center", value=False)
    update_interval = pn.widgets.FloatInput(name="Update Interval", value=1, start=0.05)
    shutdown_button = pn.widgets.Button(
        name="Shutdown the server", button_type="danger"
    )

    def shutdown(event):
        print("shutdown called")
        sys.exit(0)

    shutdown_button.on_click(shutdown)

    template.sidebar.append(run)
    template.sidebar.append(update_interval)
    template.sidebar.append(shutdown_button)

    # json snapshot `about`
    template.sidebar.append(create_json_pane(tuid, -1, snapshot_preprocess))

    # create data panel
    data_plot = DataPlot(tuid, draw_undone=draw_undone, draw_done=draw_done)
    
    update_interval.param.watch(data_plot.add_periodic_callback, 'value')
    run.param.watch(data_plot.add_periodic_callback, "value")
    
    
    pn.Column
    template.main.append(data_plot.plot)
    update_interval.param.trigger("value")
    run.param.trigger("value")

    server = pn.serve(template, show=False, port=port, threaded=True, admin=True)

    server.join()


# main()
if __name__ == "__main__":

    def data_preprocess(x: xarray.Dataset):
        x["y0"] = 20 * np.log10(np.abs(x.y0))
        return x

    main("quantify-data", "20231103-114104-952-084e5b", data_preprocess, snapshot_prune)
