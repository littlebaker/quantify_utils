from typing import Any, Callable, Literal
from qcodes.instrument.base import InstrumentBase
from qcodes.parameters import Parameter, DelegateParameter
from qcodes.validators import Validator


class BatchedParam(DelegateParameter):
    def __init__(self, param, batch_points=None) -> None:
        super().__init__(param.name, param, metadata=param.metadata)

        self.batched = True
        if batch_points is not None:
            self.batch_size = batch_points


class BatchedGettable:
    def __init__(self, gettable, batch_points=None) -> None:
        self.gettable = gettable
        self.name = gettable.name
        self.label = gettable.label
        self.unit = gettable.unit

        self.batched = True
        if batch_points is not None:
            self.batch_size = batch_points

    def get(self):
        return self.gettable.get()


def MP(name, value, label="", unit=""):
    return {name: {"value": value, "label": label, "unit": unit}}


# def fitlorentz(xData,yData):
#     x00,kappa0=findFWHM(xData,yData)
#     p0=[np.max(yData),(np.min(yData)-np.max(yData))*kappa0,x00,kappa0]
#     A0,A,x0,kappa=optimize.curve_fit(lorentz,xData,yData,p0=p0)[0]
#     return A0,A,x0,kappa