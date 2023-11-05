from __future__ import annotations


import numpy as np

from collections.abc import Mapping, Sequence
from typing import Any
from qcodes.instrument import Instrument
from qcodes.parameters import ManualParameter
from deprecated import deprecated
from qcodes.validators import Numbers

from devices.VNA_keysightN5231b_4Portsand2Ports import keysight_vna


class KeysightVNA(Instrument):
    def __init__(
        self,
        name: str,
        device_name: str,
        trace: str,
        metadata: Mapping[Any, Any] | None = None,
        label: str | None = None,
    ) -> None:
        super().__init__(name, metadata, label)
        trace_lookup = {
            "S21": 21,
            "S43": 43,
            "S41": 41,
            "S23": 23,
            "all": 4,
        }
        trace_num = trace_lookup.get(trace, "S21")
        self.trace = trace
        self.vna = keysight_vna(device_name, trace_num, timeout=300)

        self.add_parameter(
            "device_name", parameter_class=ManualParameter, initial_value=device_name
        )
        self.add_parameter(
            "trace", parameter_class=ManualParameter, initial_value=trace
        )

        self.add_parameter(
            "power",
            label="Power",
            set_cmd=lambda x: self.vna.set_power(x),
            unit="dBm",
        )

        self.add_parameter(
            "ifband",
            label="IF Bandwidth",
            get_cmd="SENS:BAND?",
            get_parser=float,
            set_cmd="SENS:BAND:RES {:.2f}",
            unit="Hz",
            vals=Numbers(min_value=1, max_value=15e6),
        )

        # Number of points in a sweep
        self.add_parameter(
            "points",
            label="Points",
            get_cmd="SENS:SWE:POIN?",
            get_parser=int,
            set_cmd="SENS:SWE:POIN {}",
            unit="",
            vals=Numbers(min_value=1, max_value=100001),
        )

        # Electrical delay
        self.add_parameter(
            "electrical_delay",
            label="Electrical Delay",
            get_cmd="CALC:CORR:EDEL:TIME?",
            get_parser=float,
            set_cmd="CALC:CORR:EDEL:TIME {:.6e}ns",
            unit="ns",
            vals=Numbers(min_value=0, max_value=100000),
        )

        # Setting frequency range
        self.add_parameter(
            "fstart",
            label="Start Frequency",
            get_cmd="SENS:FREQ:STAR?",
            get_parser=lambda x: float(x) / 1e9,
            set_cmd="SENS:FREQ:STAR {}GHz",
            unit="GHz",
            vals=Numbers(min_value=300e-6, max_value=13.5),
        )
        self.add_parameter(
            "fstop",
            label="Stop Frequency",
            get_cmd="SENS:FREQ:STOP?",
            get_parser=lambda x: float(x) / 1e9,
            set_cmd="SENS:FREQ:STOP {}GHz",
            unit="GHz",
            vals=Numbers(min_value=300e-6, max_value=13.5),
        )

        self.add_parameter(
            "data",
            label="Data",
            get_cmd=lambda: self.vna.get_data(),
        )

        self.add_parameter(
            "data_S21",
            label="S21 Data",
            get_cmd=lambda: self.vna.get_data_S21(),
        )

        self.add_parameter(
            "data_S43",
            label="S43 Data",
            get_cmd=lambda: self.vna.get_data_S43(),
        )

        self.add_parameter(
            "data_S41",
            label="S41 Data",
            get_cmd=lambda: self.vna.get_data_S41(),
        )

        self.add_parameter(
            "data_S23",
            label="S23 Data",
            get_cmd=lambda: self.vna.get_data_S23(),
        )

        self.add_parameter(
            "data_all",
            label="All Data",
            get_cmd=lambda: self.vna.get_data_all(),
        )

    def write(self, cmd: str) -> None:
        return self.vna.Inst.write(cmd)

    def ask(self, cmd: str) -> str:
        return self.vna.Inst.query(cmd)

    def __getattr__(self, key: str) -> Any:
        try:
            return self.vna[key]
        except Exception:
            pass
        return super().__getattr__(key)

    @deprecated
    def _get_extract_data(self):
        x = self.vna.get_data()
        return np.real(x), np.imag(x)

    def snapshot_base(
        self,
        update: bool | None = False,
        params_to_skip_update: Sequence[str] | None = None,
    ) -> dict[Any, Any]:

        if params_to_skip_update is None:
            ptsu = []
        else:
            ptsu = (
                list(params_to_skip_update)
            )
        if "data" not in ptsu:
            ptsu.append("data")
        return super().snapshot_base(update, ptsu)


def vna_param_set(
    vna: KeysightVNA, vnapwr:float, fstart: float, fstop: float, ifband: int, points: int, edelay=0
):
    vna.power(vnapwr)
    vna.fstart(fstart)
    vna.fstop(fstop)
    vna.ifband(ifband)
    vna.points(points)
    vna.electrical_delay(edelay)


if __name__ == "__main__":
    vna = KeysightVNA("VNA1", "vna1", "S21")
    print(vna.address(), vna.trace())

    vna_param_set(vna, 6e9, 7e9, 30, 101)

    from pprint import pprint  # noqa: F401
    # pprint(vna.snapshot())

    print(vna.data_all())
    # import sys, traceback, threading
    # thread_names = {t.ident: t.name for t in threading.enumerate()}
    # print(thread_names)
    # for thread_id, frame in sys._current_frames().items():
    #     print("Thread %s:" % thread_names.get(thread_id, thread_id))
    #     traceback.print_stack(frame)
    #     print()
    # sys.exit(1)
