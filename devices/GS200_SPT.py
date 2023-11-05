# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 19:32:49 2017

@author: SPT
"""
# Current units is mA, Voltage unites is mV
import pyvisa
import time
import numpy as np


class GS_200(object):
    def __init__(self, name):
        dcnum = {
            "DC1": "TCPIP::172.25.146.71::7655::SOCKET",
            "DC2": "TCPIP::172.25.146.78::7655::SOCKET",
            "DC4": "TCPIP::172.25.146.77::7655::SOCKET",
            "DC3": "TCPIP::172.25.146.73::7655::SOCKET",
            "DC5": "TCPIP::172.25.146.75::7655::SOCKET",
            "DC6": "TCPIP::172.25.146.70::7655::SOCKET",
            "DC7": "TCPIP::172.25.146.69::7655::SOCKET",
            "DC8": "TCPIP::172.25.146.72::7655::SOCKET",
            "DC9": "TCPIP::172.25.146.79::7655::SOCKET",
            "DC10": "TCPIP::172.25.146.74::7655::SOCKET",
        }
        Rm = pyvisa.ResourceManager()
        self.lastLevel = 0
        self.Inst = Rm.open_resource(
            dcnum[name], write_termination="\n", read_termination="\n"
        )
        print(self.Inst.query("*IDN?"))
        print(name + " Connect OK")
        #        self.Inst.write('*RST')
        self.Inst.write(":SOUR:FUNC CURR")

    def rst(self):
        self.Inst.write("*RST")

    def setCURRmode(self):
        self.Inst.write(":SOUR:FUNC CURR")

    def setslopetime(self, times):
        self.Inst.write(":PROG:SLOP {}".format(time))

    def setintervaltime(self, times):
        self.Inst.write(":PROG:INT {}".format(time))

    def setstep(self, step):
        self.Inst.write(":PROG:STEP {}".format(step))

    def setVOLTmode(self):
        self.Inst.write(":SOUR:FUNC VOLT")  # default

    def setLevel(self, level=0.0):
        self.Inst.write(":SOUR:LEV {}".format(level / 10**3))

    def getlevel(self):
        return float(self.Inst.query(":SOUR:LEV?"))

    def setlevel_slow(self, v1, velocity):  # 单位mA/s
        v1 = v1 * 10**-3
        velocity = velocity * 10**-3
        v0 = float(self.Inst.query(":SOUR:LEV?"))
        delta = v1 - v0
        steps = abs(int(delta / (velocity))) + 1
        if steps == 1:
            steps = 2
        levlist = np.linspace(v0, v1, steps)
        for i in range(len(levlist)):
            self.Inst.write(":SOUR:LEV {}".format(levlist[i]))
            time.sleep(1)

    def setlevel_slowV2(self, v1, velocity):
        v1 = v1 * 10**-3
        velocity = velocity * 10**-3
        v0 = float(self.Inst.query(":SOUR:LEV?"))
        t1 = int(np.abs(v1 - v0) / velocity)
        print(t1)
        if t1 >= 1 and t1 < 3600:
            aa = self.Inst.write(":PROG:EDIT:START")
            #        aa+=self.Inst.write('::SOUR:FUNC VOLT')
            # aa+=self.Inst.write(':SOUR:LEV {}'.format(v0))
            aa += self.Inst.write(":SOUR:LEV {}".format(v1))
            aa += self.Inst.write(":PROG:EDIT:END")
            aa += self.Inst.write(":PROG:INT {}".format(t1))
            aa += self.Inst.write(":PROG:SLOP {}".format(t1))
            aa += self.Inst.write(":PROG:REP OFF")
            aa += self.Inst.write(":PROG:RUN")
            # print(time.strftime("%H:%M:%S"))
            for i in range(int(t1 + 1)):
                time.sleep(1)
        elif t1 > 3600:
            print("time is too long")
            # vs=float(self.Inst.query(':SOUR:LEV?'))
            # print(vs)
            # print(time.strftime("%H:%M:%S"))

    def StopSetLevel(self):
        self.Inst.write(":PROG:PAUS")

    def setRange_Lev(self, rang=1.2):
        rang = rang * 10**-3
        self.Inst.write(":SOUR:RANG {}".format(rang))

    def getRange_Lev(self):
        return self.Inst.query(":SOUR:RANG?")

    def Start_OutPut(self):
        self.Inst.write(":OUTP 1")

    def Stop_OutPut(self):
        self.Inst.write(":OUTP 0")

    def Close(self):
        self.Inst.close()

    def program(self):
        aa = self.Inst.write(":PROG:EDIT:START")
        #        aa+=self.Inst.write('::SOUR:FUNC VOLT')
        aa += self.Inst.write(":SOUR:LEV {}".format(0.1))
        aa += self.Inst.write(":SOUR:LEV {}".format(0.2))
        aa += self.Inst.write(":SOUR:LEV {}".format(0.3))
        aa += self.Inst.write(":SOUR:LEV {}".format(0.1))
        aa += self.Inst.write(":PROG:EDIT:END")
        aa += self.Inst.write(":PROG:INT {}".format(1))
        aa += self.Inst.write(":PROG:SLOP {}".format(1))
        #        aa+=self.Inst.write(':PROG:REP OFF')
        aa += self.Inst.write(":PROG:RUN")
        # DC1.setLevel(0.1)
        print(aa)  # %%

    def clear(self):
        self.Inst.write("*CLS")
