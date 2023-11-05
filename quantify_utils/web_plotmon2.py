from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Callable

from quantify_core.data.handling import get_datadir
from qcodes.instrument import Instrument
from multiprocessing import Process
from .web_monitor2 import main, snapshot_prune


class PlotMonitor_web(Instrument):
    def __init__(
        self,
        name: str,
        metadata: Mapping[Any, Any] | None = None,
        label: str | None = None,
        draw_undone: Callable | None = None,
        draw_done: Callable | None = None,
    ) -> None:
        super().__init__(name, metadata, label)
        self.datadir = get_datadir()
        self.tuid = None
        self.proc = None
        self.draw_undone = draw_undone
        self.draw_done = draw_done
        self.port = 5006

    def update(self, tuid):
        self.tuid = tuid
        if self.proc is not None:
            self.proc.terminate()
            self.proc = None
        self.proc = Process(
            target=main,
            args=(self.datadir, tuid, self.draw_undone, self.draw_done, snapshot_prune, self.port),
            daemon=True,
        )
        self.proc.start()

        return self.proc


if __name__ == "__main__":
    import quantify_core.data.handling as dh

    dh.set_datadir("/Users/lbaker/Projects/PycharmProjects/qcodes_test/quantify-data")
    pm = PlotMonitor_web("pm")
    pm.update("20231028-225013-687-88b211")
