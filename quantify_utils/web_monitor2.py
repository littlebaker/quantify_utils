import sys
import time
import panel as pn
import quantify_core.data.handling as dh
from quantify_core.data.handling import load_snapshot, to_gridded_dataset
import hvplot.xarray  # noqa: F401
import numpy as np
import xarray
import holoviews as hv  # noqa: F401


from quantify_core.measurement.control import _DATASET_LOCKS_DIR
from quantify_core.visualization.pyqt_plotmon_remote import _safe_load_dataset
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

    def create_data_plot(self, interval, run, cmap):
        is_finished = False
        if self.tuid == "":
            yield None
        while True:
            if not run:
                yield "Calculation did not run yet"
                return
            time.sleep(interval)

            dset = _safe_load_dataset(self.tuid, _DATASET_LOCKS_DIR)
            is_finished = not np.any(dset.y0.isnull().to_numpy())

            if dset.attrs["grid_2d"]:  # 2d-data
                dset_grid = to_gridded_dataset(dset)

                if not is_finished:
                    fig = None
                    if self.draw_undone is not None:
                        fig = self.draw_undone(dset_grid)
                    else:
                        fig = self._default_draw_undone(dset_grid)
                    yield pn.Column(fig)
                else:
                    fig = None
                    if self.draw_done is not None:
                        fig = self.draw_done(dset_grid)
                    else:
                        fig = self._default_draw_done(dset_grid)
                    yield pn.Column(fig, pn.pane.Str("Finished"))
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

                    yield pn.Column(fig)
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

                        yield fig

            if is_finished:
                return

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


def main(datadir, tuid, draw_undone=None, draw_done=None, snapshot_preprocess=None):
    dh.set_datadir(datadir)
    template = pn.template.BootstrapTemplate(title="Web Monitor")
    run = pn.widgets.Checkbox(name="Press to run calculation", align="center")
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

    data_plot = DataPlot(tuid, draw_undone=draw_undone, draw_done=draw_done)
    plot = pn.bind(data_plot.create_data_plot, update_interval, run, "bwr")
    template.main.append(plot)

    server = pn.serve(template, show=False, port=5006, threaded=True)
    server.join()


# main()
if __name__ == "__main__":

    def data_preprocess(x: xarray.Dataset):
        x["y0"] = 20 * np.log10(np.abs(x.y0))
        return x

    main("quantify-data", "20231103-114104-952-084e5b", data_preprocess, snapshot_prune)
