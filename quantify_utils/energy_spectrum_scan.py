from __future__ import annotations

from types import SimpleNamespace
from typing import Iterable
import numpy as np
from .KeysightVNA import KeysightVNA, vna_param_set
from .GS200 import YokogawaGS, gs_param_set

from qcodes.instrument import find_or_create_instrument
from qcodes.parameters import ManualParameter
from quantify_core.measurement import MeasurementControl
from .measurement_utils import BatchedParam, BatchedGettable



def checkbox(timeused):
    print(f"大概需要{timeused}秒({round(timeused/60)}分钟)")
    print("请输入Y/N确认程序执行:")
    in_content = input("请输入：").upper()
    if in_content != "Y":
        print("User abort.")
        return False
    elif in_content == "Y":
        print("Continuing executing...")
        return True
    



def Energy_spectrum_scan_2D(
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
    vnapower: float,
    measurement_name="Energy spectrum 2D",
    plotmon: str | None = None
):
    timeused = round(
        (points / ifband + (volt[-1] - volt[0])/(len(volt)) / gs_step) * (len(volt) + 1), 1
    )
    # if not checkbox(timeused):
        # return None

    freq_list = np.linspace(fstart, fstop, points)

    vna_param_set(vna, vnapower, fstart, fstop, ifband, points, elec_delay)
    gs_param_set(gs, gs_step, "CURR", "v2")
    gs.level(volt[0])

    about = {
        "Attenuation": {
            "vna": f"{RTatt}+{-vnapower}dB",
            "inline": "60dB=42+8+10dB",
            "outline": "76dB(RT)+36dB(4K)",
        },
        "VNA": {
            "fstart": fstart,
            "fstop": fstop,
            "points": points,
            "ifband": f"{ifband}Hz",
            "vnapower": f"{vnapower}dB",
            "elec_delay": elec_delay,
        },
        "DC": {"current": f"from {volt[0]}mA to {volt[-1]}mA"},
    }

    meas_ctrl = find_or_create_instrument(
        MeasurementControl,
        "meas_ctrl",
    )
    meas_ctrl.metadata["about"] = about
    
    if plotmon is not None:
        meas_ctrl.instr_plotmon(plotmon)

    temp_freq_param = ManualParameter("Frequency", unit="GHz")
    meas_ctrl.settables([BatchedParam(temp_freq_param, points), gs.level])
    meas_ctrl.setpoints_grid([freq_list, volt])
    
    def _get_extract_data():
        x = vna.data()
        if vna.trace == "all":
            sdata21,sdata43,_,_ = x
            return np.real(sdata21), np.imag(sdata21), np.real(sdata43), np.imag(sdata43),
                
        else:
            return np.real(x), np.imag(x)


    if vna.trace == "all":
        gettable = SimpleNamespace(  # A Gettable-compatible object
            name=["s21_real", "s21_imag", "s43_real", "s43_imag",],
            get=_get_extract_data,
            label=["re(s21)", "im(s21)", "re(s43)", "im(s43)",],
            unit=["V", "V", "V", "V",],
            batched=True,
        )
    else:
        gettable = SimpleNamespace(  # A Gettable-compatible object
            name=["z_real", "z_imag"],
            get=_get_extract_data,
            label=["re(z)", "im(z)"],
            unit=["V", "V"],
            batched=True,
        )


    meas_ctrl.gettables(
        BatchedGettable(gettable, points),
    )

    dset = meas_ctrl.run(measurement_name)
    return dset


if __name__ == "__main__":
    vna_name = "vna2"
    vna_trace = "S21"
    gs_name = "DC2"
    vna = find_or_create_instrument(
        KeysightVNA, name=vna_name, device_name=vna_name, trace=vna_trace
    )
    gs = find_or_create_instrument(YokogawaGS, name=gs_name, device_name=gs_name)
    dset, da_amp, da_phase = Energy_spectrum_scan_2D(
        RTatt=60,
        vna=vna,
        gs=gs,
        gs_step=1e-3,
        elec_delay=0,
        volt=[0],
        fstart=6.4,
        fstop=6.6,
        points=101,
        ifband=20,
        vnapower=-15,
        measurement_name="Left Cavity freq scan",
    )

    # from visualization import heatmap_with_colorbar_range
    # heatmap_with_colorbar_range(dset_abs)
