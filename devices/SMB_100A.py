# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 09:57:21 2017

@author: andy
"""

# from PyQt5.QtWidgets import QMessageBox
import pyvisa
import numpy as np


class MW(object):
    _instance = None
    clsname = []

    def __init__(self, name):
        Rm = pyvisa.ResourceManager()
        mwnum = {
            "mw1": "TCPIP0::172.25.146.80::5025::SOCKET",
            "mw2": "TCPIP0::172.25.146.82::5025::SOCKET",
            "mw3": "TCPIP0::172.25.146.81::5025::SOCKET",
        }
        try:
            self.Inst = Rm.open_resource(mwnum[name])
            print(name + " Connect OK")
            self.Inst.write("*RST")
        except Exception:
            print("找不到微波源")
            pass
        
        # For Serial and TCP/IP socket connections enable the read Termination Character, or read's will timeout
        if self.Inst.resource_name.startswith('ASRL') or self.Inst.resource_name.endswith('SOCKET'):
            self.Inst.read_termination = '\n'
        

    def MW_set(self, fre=6e9, amp=-20, pha=0):
        #        print("nini")
        self.Inst.write("FREQ:CW {}GHz".format(fre))
        self.Inst.write(":POW {}".format(amp))
        self.Inst.write("PHAS {}DEG".format(pha))

    def MW_setFre(self, fre):
        #        print("nini")
        self.Inst.write("FREQ:CW {}GHz".format(fre))

    def MW_setpower(self, amp):
        #        print("nini")
        self.Inst.write(":POW {}".format(amp))

    def start_Output(self):
        self.Inst.write("OUTP ON")

    def MW_get(self):
        #        print("ffi")
        str_Fre = str(float(self.Inst.query("FREQ:CW?")) / 1e9)
        str_Amp = self.Inst.query(":POW?")
        str_Pha = self.Inst.query(":PHAS?")
        str_On_off = self.Inst.query("OUTP?")
        return str_Fre, str_Amp, str_Pha, str_On_off

    def Close(self):
        self.Inst.close()

    def Close_Output(self):
        self.Inst.write("OUTP OFF")

    def VnaTrig(self, f0, f1, steps):
        freqlist = np.linspace(f0, f1, steps)
        A = str()
        for i in range(len(freqlist)):
            if i != len(freqlist) - 1:
                A = A + str(freqlist[i]) + "GHZ,"
            else:
                A = A + str(freqlist[i]) + "GHZ"  # construct the freq list str
        self.Inst.write(":LIST:FREQ " + A)
        self.Inst.write("LIST:TRIG:SOUR EXT")  #
        self.Inst.write("LIST:TYPE LIST")
        self.Inst.write(":INIT:CONT ON")
        self.Inst.write(
            "SWE:CONT:STAT 1"
        )  # The preceding example sets the sweep control state to on
        self.Inst.write(
            "FREQ:MODE LIST"
        )  # The preceding example selects a list frequency sweep

    # def Trig_mode_init(self, f0, amp):

    #     # A=str(f0)+'GHZ'

    #     self.Inst.write('FREQ:CW {}GHz'.format(f0))
    #     self.Inst.write(':POW {}'.format(amp))
    #     self.Inst.write('TRIG:SOUR EXT') #
    #     # self.Inst.write('LIST:TYPE LIST')
    #     self.Inst.write(':INIT:CONT ON')
    #     #self.Inst.write('SWE:CONT:STAT 1') #The preceding example sets the sweep control state to on
    #     # self.Inst.write('FREQ:MODE LIST') #The pr

    # %%


if __name__ == "__main__":
    mw = MW("mw1")
    print(id(mw))
    mw11 = mw
    print(mw.MW_get())
    mw11.Close()

# %%
# Rm=pyvisa.ResourceManager()
# Inst=Rm.open_resource('TCPIP0::172.25.146.80::inst0::INSTR')
# Inst.query('FREQ:MODE?')
# Inst.query('LIST:DIR?')
# Inst.query('LIST:FREQ?')
# Inst.query('LIST:FREQ:POIN?')
# Inst.query('LIST:MODE?')
# Inst.query('LIST:TRIG:SOUR?')
# Inst.query('LIST:TYPE?')
# Inst.query('LIST:RETR?')
# Inst.query('SWE:CONT:STAT?')
# Inst.query('LIST:MAN?')
#
# Inst.write(':LIST:FREQ 4Ghz,5Ghz')
#
# Inst.write('*RST')
# Inst.write('LIST:TRIG:SOUR EXT')
# Inst.write('FREQ:MODE LIST')
# Inst.write('LIST:TYPE LIST')
# Inst.write('SWE:CONT:STAT 1')
#
#
# Inst.write('SWE:CONT:STAT 1')
#
#
#
#
#
#
# Inst.write('OUTP OFF')
