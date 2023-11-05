from __future__ import annotations

from collections.abc import Mapping
from logging import getLogger
from typing import Any
from qcodes.instrument import Instrument
from qcodes.parameters import ManualParameter

from qcodes.validators import Enum, Numbers

from devices.GS200_SPT import GS_200


class YokogawaGS(Instrument):
    def __init__(self, name: str, device_name: str,  
        metadata: Mapping[Any, Any] | None = None, label: str | None = None) -> None:
        super().__init__(name, metadata, label)
        self.logger = getLogger()
        self.gs = GS_200(device_name)
        
        self.add_parameter(
            "device_name",
            parameter_class=ManualParameter,
            initial_value=device_name
        )
        
        self.add_parameter(
            "mode",
            label='Mode',
            set_cmd=':SOUR:FUNC {}',
            get_cmd=':SOUR:FUNC?',
            vals=Enum('CURR', "VOLT")
        )
        

        
        self.add_parameter(
            "level_set_mode",
            parameter_class=ManualParameter,
            label="Level change speed",
            initial_value="v2",
            vals=Enum("v1", "v2")
        )
        
        self.add_parameter(
            "step",
            parameter_class=ManualParameter,
            label="Level step",
            initial_value=0.0,
            vals=Numbers(0)
        )
        
        self.add_parameter(
            "level",
            label="level",
            set_cmd=lambda x: self._set_level(x, self.step()),
            get_cmd=':SOUR:LEV?',
            get_parser=lambda x: float(x)*10**3 ,
            unit="mA/mV"
        )
        
    def _set_level(self, v1, velocity):
        if self.level_set_mode() == 'v1':
            return self.gs.setlevel_slow(v1, velocity)
        elif  self.level_set_mode() == 'v2':
            return self.gs.setlevel_slowV2(v1, velocity)
        
    def start_output(self, start=True):
        assert start in [True, False]
        if start:
            self.gs.Start_OutPut()
        else:
            self.gs.Stop_OutPut()
        
    def reset(self):
        self.write('*RST')
        
    def write(self, cmd: str) -> None:
        return self.gs.Inst.write(cmd)
    
    def ask(self, cmd: str) -> str:
        return self.gs.Inst.query(cmd)
    
    def __getattr__(self, key: str) -> Any:
        try:
            return self.gs[key]
        except Exception as e: 
            self.logger.warning(e)
        return super().__getattr__(key)
    
    
def gs_param_set(gs: YokogawaGS, step,  mode="CURR", level_mode="v2"):
    gs.mode(mode)
    gs.level_set_mode(level_mode)
    gs.start_output(True)
    gs.step(step)