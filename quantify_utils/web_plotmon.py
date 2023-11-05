from __future__ import annotations

from collections.abc import Mapping
import time
from typing import Any
import warnings

import pyqtgraph.multiprocess as pgmp
from qcodes import Instrument, Parameter
from qcodes import validators as vals
from qcodes.utils.helpers import strip_attrs

from quantify_core.data.handling import get_datadir
from quantify_core.measurement.control import _DATASET_LOCKS_DIR
import subprocess as sp

class PlotMonitor_web(Instrument):
    
    def __init__(self, name: str, metadata: Mapping[Any, Any] | None = None, label: str | None = None) -> None:
        super().__init__(name, metadata, label)
        
        # self.proc = sp.run(["panel", "serve", "/Users/lbaker/Projects/PycharmProjects/qcodes_test/web_monitor.py", "--port", "5006"], capture_output=True)
        self.datadir = get_datadir()
        self.tuid = None
        
        
    def update(self, tuid):
        self.tuid = tuid 
        self.proc = sp.Popen(["panel", "serve", "/Users/lbaker/Projects/PycharmProjects/qcodes_test/web_monitor.py", "--port", "5006", "--args",
            f"{self.tuid}", f"{self.datadir}"])
        
        
if __name__ == "__main__":
    import quantify_core.data.handling as dh
    dh.set_datadir("/Users/lbaker/Projects/PycharmProjects/qcodes_test/quantify-data")
    pm = PlotMonitor_web("pm")
    pm.update("20231028-225013-687-88b211")
    