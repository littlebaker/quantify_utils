from __future__ import annotations

from types import SimpleNamespace
import numpy as np

import xarray as xr

from .KeysightVNA import KeysightVNA, vna_param_set


from qcodes.instrument import find_or_create_instrument
from qcodes.parameters import ManualParameter
from quantify_core.measurement import MeasurementControl

from .measurement_utils import BatchedParam, BatchedGettable
from quantify_core.data.experiment import QuantifyExperiment




def phase_handle(Phase, f1, phase_delay):
    # phase_delay=67.86#62.03
    phase_offset = 0
    if len(Phase.shape) == 1:
        _phase = np.array(Phase[:, None])
    else:
        _phase = np.array(Phase)
    phase_correct = np.zeros((_phase.shape[0], _phase.shape[1]))

    for i in range(len(f1)):
        phase_temp = _phase[i] + 2 * np.pi * f1[i] * phase_delay + phase_offset
        phase_correct[i] = np.angle(np.cos(phase_temp) + 1j * np.sin(phase_temp))
    return phase_correct


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





def concat_data(dset):
    _y = dset.y0.data + 1j * dset.y1.data
    _x = dset.x0.data
    return _x, _y

def _amp_process(_y):
    S_dataplot_amp = 20 * np.log10(np.abs(_y))
    return S_dataplot_amp

def _phase_process(_y, _x, elec_delay):
    Phase = np.angle(_y)
    S_dataplot_phase = phase_handle(Phase, _x, elec_delay)
    return S_dataplot_phase


def vna_scan_1D(
    RTatt: int,
    vna: KeysightVNA,
    elec_delay: float,
    fstart: float,
    fstop: float,
    points: int,
    ifband: int,
    vnapower: int,
    measurement_name="Energy spectrum 2D",
):
    timeused = round((points / ifband) * 1, 2)
    if not checkbox(timeused):
        return None

    freq_list = np.linspace(fstart, fstop, points)

    vna_param_set(vna, fstart, fstop, ifband, points)

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
    }

    meas_ctrl = find_or_create_instrument(
        MeasurementControl,
        "meas_ctrl",
    )
    meas_ctrl.metadata["about"] = about

    temp_freq_param = ManualParameter("Frequency", unit="GHz")
    meas_ctrl.settables([BatchedParam(temp_freq_param, points)])
    meas_ctrl.setpoints(freq_list)

    def _get_extract_data():
        x = vna.data()
        return np.real(x), np.imag(x)

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

    # post process
    _x, _y = concat_data(dset)
    _amp = _amp_process(_y)
    _phase = _phase_process(_y, _x, elec_delay=elec_delay).squeeze(0)
    dset["amp"] = xr.DataArray(_amp, coords=dict(x0=_x))
    dset["phase"] = xr.DataArray(_phase, coords=dict(x0=_x))

    exp = QuantifyExperiment(dset.tuid)
    exp.write_dataset(dset)

    return dset


if __name__ == "__main__":


    dset = vna_scan_1D(
        RTatt=60,
        vna_name = "vna2",
        vna_trace = "S21",
        elec_delay=50.14,
        fstart=6.4,
        fstop=6.5,
        points=101,
        ifband=20,
        vnapower=-15,
        measurement_name="Left Cavity freq scan",
    )

