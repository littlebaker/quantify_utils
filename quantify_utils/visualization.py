import operator
from functools import reduce
import panel as pn
import hvplot.xarray  # noqa: F401
import holoviews as hv
import numpy as np
import xarray as xr
from .vna_scan import phase_handle

from .analysis import _cal_spec



def spec_draw_amp(dset: xr.Dataset):
    dset_amp = _cal_spec(dset, "amp")

    # add all plots up
    return reduce(
        operator.add,
        (
            dset_amp.hvplot.image(x="x1", y="x0", z=str(key), dynamic=False, cmap="bwr")
            for key in dset_amp.data_vars.keys()
        ),
    )


def spec_draw_phase(dset: xr.Dataset, elec_delay=0):
    dset_amp = _cal_spec(dset, "phase", elec_delay)

    # add all plots up
    return reduce(
        operator.add,
        (
            dset_amp.hvplot.image(x="x1", y="x0", z=str(key), dynamic=False, cmap="bwr")
            for key in dset_amp.data_vars.keys()
        ),
    )


def spec_draw_widget_amp(dset: xr.Dataset):
    amp = xr.DataArray(
        20 * np.log10(np.abs(dset.y0.values + dset.y1.values * 1j)), coords=dset.coords
    )
    dset_copy = dset.copy()
    dset_copy["amp"] = amp

    return heatmap_crange_cmap(dset_copy, "amp")


def spec_draw_widget_phase(dset: xr.Dataset, elec_delay=0):
    phase = xr.DataArray(
        phase_handle(
            np.angle(dset.y0.values + dset.y1.values * 1j), dset.x0.values, elec_delay
        ),
        coords=dset.coords,
    )
    dset_copy = dset.copy()
    dset_copy["phase"] = phase

    return heatmap_crange_cmap(dset_copy, "phase")


def heatmap_with_colorbar_range(dset_grid, data_var="y0", backend=None, cmap="bwr"):
    low, high = np.min(dset_grid[data_var].data), np.max(dset_grid[data_var].data)
    print(low, high)
    color_bar_range = pn.widgets.RangeSlider(
        name="color bar range", start=low, end=high
    )

    def create_heatmap(bounds):
        return hv.render(
            dset_grid[data_var].hvplot.image(dynamic=True, clim=bounds, cmap=cmap),
            backend=backend,
        )

    main_pane = pn.bind(create_heatmap, color_bar_range)

    return pn.Column(main_pane, color_bar_range, sizing_mode="stretch_both")


def heatmap_crange_cmap(dset_grid, data_var="y0", backend=None, cmap="bwr"):
    cmap_selector = pn.widgets.Select(
        name="color bar selector", value=cmap, options=["bwr", "hot", "hot_r", "jet"]
    )
    low, high = np.min(dset_grid[data_var].data), np.max(dset_grid[data_var].data)
    print(low, high)
    color_bar_range = pn.widgets.RangeSlider(
        name="color bar range", start=low, end=high
    )

    def create_heatmap(bounds, cmap):
        return hv.render(
            dset_grid[data_var].hvplot.image(dynamic=True, clim=bounds, cmap=cmap),
            backend=backend,
        )

    main_pane = pn.bind(create_heatmap, color_bar_range, cmap_selector)

    return pn.Column(main_pane, color_bar_range, cmap_selector)
