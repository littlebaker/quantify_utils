from qcodes.instrument import find_or_create_instrument
from quantify_utils.KeysightVNA import KeysightVNA
from quantify_utils.GS200 import YokogawaGS
from typing import Iterable
from quantify_utils.vna_scan import vna_scan_1D
from quantify_utils.analysis import LorentzianAnalysis, SpectrumAnalysis
from quantify_utils.energy_spectrum_scan import Energy_spectrum_scan_2D
import quantify_core.data.handling as dh
from quantify_utils.web_plotmon2 import PlotMonitor_web
from quantify_utils.visualization import (
    spec_draw_amp,
    spec_draw_phase,
    spec_draw_widget_amp,
    spec_draw_widget_phase,
)
import panel as pn


def spectrum_scan(
    RTatt: int,
    vna: KeysightVNA,
    gs: YokogawaGS,
    gs_step: float,
    elec_delay: float,
    volt: Iterable[float],
    fstart: float,
    fstop: float,
    points: int,
    ifband: int,
    vnapower: int,
    measurement_name="Energy spectrum 2D",
    analysis_model=SpectrumAnalysis,
):
    def _draw_undone(x):
        return  pn.Column(spec_draw_amp(x), spec_draw_phase(x, elec_delay))
    def _draw_done(x):
        return  pn.Column(spec_draw_widget_amp(x), spec_draw_widget_phase(x, elec_delay))
    
    plotmon = find_or_create_instrument(
        PlotMonitor_web,
        "plotmon",
        draw_undone=_draw_undone,
        draw_done=_draw_done,
    )
    dset = Energy_spectrum_scan_2D(
        RTatt,
        vna,
        gs,
        gs_step,
        elec_delay,
        volt,
        fstart,
        fstop,
        points,
        ifband,
        vnapower,
        measurement_name,
        plotmon=plotmon.name,
    )
    
    ana_model = analysis_model(dset, elec_delay=elec_delay)
    ana_model.run()
    
    
    return dset


def scan_cavity_freq(
    vna_name: str,
    vna_trace: str,
    RTatt: int,
    elec_delay: float,
    fstart: float,
    fstop: float,
    points: int,
    ifband: int,
    vnapower: float,
    measurement_name: str,
    analysis_model=LorentzianAnalysis,
):
    vna = find_or_create_instrument(
        KeysightVNA, name=vna_name, device_name=vna_name, trace=vna_trace
    )
    dset = vna_scan_1D(
        RTatt=RTatt,
        vna=vna,
        elec_delay=elec_delay,
        fstart=fstart,
        fstop=fstop,
        points=points,
        ifband=ifband,
        vnapower=vnapower,
        measurement_name=measurement_name,
    )
    if analysis_model is not None:
        analysis = analysis_model(dset)
        a_obj = analysis.run()
        a_obj.display_figs_mpl()


if __name__ == "__main__":
    import numpy as np

    dh.set_datadir("./quantify-data/CELv1")
    vna_name = "vna2"
    vna_trace = "all"
    gs_name = "DC7"
    vna = find_or_create_instrument(KeysightVNA, vna_name, vna_name, vna_trace)
    gs = find_or_create_instrument(YokogawaGS, gs_name, gs_name)
    # scan_cavity_freq(
    #     RTatt=60,
    #     vna_name = "vna2",
    #     vna_trace = "S21",
    #     elec_delay=60.05,
    #     fstart=6.833-0.01,
    #     fstop=6.833+0.01,
    #     points=101,
    #     ifband=20,
    #     vnapower=-15,
    #     measurement_name="Left Cavity freq scan",
    # )

    spectrum_scan(
        RTatt=60,
        vna=vna,
        gs=gs,
        gs_step=0.0001,
        elec_delay=0,
        volt=np.linspace(0, -1, 1001),
        fstart=6.45 - 0.025,
        fstop=6.45 + 0.025,
        points=51,
        ifband=40,
        vnapower=-15,
        measurement_name="Left cavity Global with reflection",
    )
