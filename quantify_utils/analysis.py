from __future__ import annotations

import matplotlib.pyplot as plt
import os
import lmfit
from quantify_core.analysis.fitting_models import LorentzianModel, lorentzian_func
import quantify_core.data.handling as dh
import quantify_core.analysis.base_analysis as ba
from xarray import Dataset
import xarray as xr
from quantify_core.data.types import TUID
import numpy as np

from quantify_utils.vna_scan import phase_handle


def _cal_spec(dset: xr.Dataset, cal_type: str, elec_delay=0):
    nof_channel = int(len(dset.data_vars) / 2)

    tp = cal_type

    if tp == "amp":
        dvars = {
            tp + str(i): xr.DataArray(
                20
                * np.log10(
                    np.abs(
                        dset["y" + str(2 * i)].values
                        + dset["y" + str(2 * i + 1)].values * 1j
                    )
                ),
                coords=dset.coords,
            )
            for i in range(nof_channel)
        }
    else:
        dvars = {
            tp + str(i): xr.DataArray(
                phase_handle(
                    np.angle(
                        dset["y" + str(2 * i)].values
                        + dset["y" + str(2 * i + 1)].values * 1j
                    ),
                    dset.x0.values,
                    elec_delay,
                ),
                coords=dset.coords,
            )
            for i in range(nof_channel)
        }
    dset_new = xr.Dataset(dvars)

    return dset_new


class LorentzianAnalysis(ba.BasicAnalysis):
    def __init__(
        self,
        dataset: Dataset = None,
        tuid: TUID | str = None,
        label: str = "",
        settings_overwrite: dict = None,
        plot_figures: bool = True,
        positive=False,
    ):
        self.positive = positive
        super().__init__(dataset, tuid, label, settings_overwrite, plot_figures)

    def run_fitting(self):
        model = LorentzianModel()
        _y = (1 if self.positive else -1) * self.dataset.amp.values
        _x = self.dataset.x0.values
        guess = model.guess(_y, x=_x)
        result = model.fit(_y, x=_x, params=guess)
        result.params["a"].value = (1 if self.positive else -1) * result.params[
            "a"
        ].value
        result.params["c"].value = (1 if self.positive else -1) * result.params[
            "c"
        ].value
        self.fit_results.update({"lorentzian": result})

    def create_figures(self):
        _x = self.dataset.x0.values
        _amp = self.dataset.amp.values
        _phase = self.dataset.phase.values

        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        fig_id = "Cavity scan"
        self.figs_mpl.update({fig_id: fig})
        self.axs_mpl.update({fig_id + "_" + str(i): axes[i] for i in range(len(axes))})

        fitres = self.fit_results["lorentzian"]

        def pv(x):
            return fitres.params.get(x)

        x0, width, a, c = pv("x0"), pv("width"), pv("a"), pv("c")

        _amp_predict = lorentzian_func(_x, x0, width, a, c)

        axes[0].plot(_x, _amp)
        axes[0].plot(_x, _amp_predict, "r--")
        axes[0].set_title(rf"Amp $\kappa$={round(width*1000,3)}MHz")
        axes[0].set_xlabel("Freqency (GHz)")
        axes[0].set_ylabel(r"Amp (dB)")
        axes[0].legend(["original", "fit"])
        axes[1].plot(_x, _phase)
        axes[1].set_title("Phase")
        axes[1].set_xlabel("Freqency (GHz)")
        axes[1].set_ylabel("Degree")

    def save_fit_results(self):
        for fr_name, fit_result in self.fit_results.items():
            path = os.path.join(self.results_dir, f"{fr_name}-report.txt")
            with open(path, "a") as f:
                f.write(lmfit.fit_report(fit_result))

        return super().save_fit_results()


class SpectrumAnalysis(ba.Basic2DAnalysis):
    def __init__(
        self,
        dataset: Dataset = None,
        tuid: TUID | str = None,
        label: str = "",
        settings_overwrite: dict = None,
        plot_figures: bool = True,
        elec_delay: float = 0,
    ):
        self.elec_delay = elec_delay
        super().__init__(dataset, tuid, label, settings_overwrite, plot_figures)

    def process_data(self):
        dset_grid = dh.to_gridded_dataset(self.dataset)
        _amps = _cal_spec(dset_grid, "amp")
        _phases = _cal_spec(dset_grid, "phase", elec_delay=self.elec_delay)

        data_vars = {}
        data_vars.update(_amps)
        data_vars.update(_phases)
        dset_processed = Dataset(data_vars)

        self.dataset_processed = dset_processed

    def create_figures(self):
        nof_channel = int(len(self.dataset_processed) / 2)
        for i in range(nof_channel):
            fig, axes = plt.subplots(1, 2, figsize=(10, 4))
            fig_id = "Spectrum_scan_" + str(i)
            self.figs_mpl.update({fig_id: fig})
            self.axs_mpl.update(
                {fig_id + "_" + str(i): axes[i] for i in range(len(axes))}
            )

            _freq = self.dataset_processed.x0.values
            _current = self.dataset_processed.x1.values

            im1 = axes[0].imshow(
                self.dataset_processed["amp" + str(i)].values.T,
                cmap="bwr",
                extent=[_current[0], _current[-1], _freq[0], _freq[-1]],
            )
            im2 = axes[1].imshow(
                self.dataset_processed["phase" + str(i)].values.T,
                cmap="bwr",
                extent=[_current[0], _current[-1], _freq[0], _freq[-1]],
            )
            axes[0].set_xlabel("Current(mA)")
            axes[0].set_ylabel("Frequency(GHz)")
            axes[0].set_title("Amplitude")
            axes[1].set_xlabel("Current(mA)")
            axes[1].set_ylabel("Frequency(GHz)")
            axes[1].set_title("Phase")
            
            fig.suptitle(f"Channel {i}")

            cbar1 = fig.colorbar(im1, ax=axes[0])
            cbar2 = fig.colorbar(im2, ax=axes[1])
            cbar1.set_label("dB")
            cbar2.set_label("rad")
