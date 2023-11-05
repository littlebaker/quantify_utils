import pyvisa
import numpy as np
from pyvisa import constants

"""
Created on Mon Oct 23 21:07:40 2017

@author: SPT
"""
# import visa


class keysight_vna(object):
    def __init__(self, name, trace=1, timeout=300):
        Rm = pyvisa.ResourceManager()
        vnanum = {
            # "vna1": "TCPIP::172.25.146.85::inst0::INSTR",
            # "vna2": "TCPIP::172.25.146.89::inst0::INSTR",
            "vna1": "TCPIP::172.25.146.85::5025::SOCKET",
            "vna2": "TCPIP::172.25.146.89::5025::SOCKET",
        }
        self.trace = trace
        self.Inst = Rm.open_resource(
            vnanum[name],
            open_timeout=timeout,
            access_mode=constants.AccessModes.shared_lock,
        )
        # For Serial and TCP/IP socket connections enable the read Termination Character, or read's will timeout
        if self.Inst.resource_name.startswith('ASRL') or self.Inst.resource_name.endswith('SOCKET'):
            self.Inst.read_termination = '\n'

        self.Inst.write(
            "*RST",
        )
        print(name + " Connect OK")
        self.Inst.write("SYSTem:PREset")
        self.Inst.write("SYSTem:FPReset")
        self.Inst.write("INIT:CONT OFF")
        self.Inst.query("*OPC?")
        if trace == 4:
            self.Inst.write("CALCulate1:PARameter:DEFine " "MySMeaS21" ",S21")
            self.Inst.write("CALCulate1:PARameter:DEFine " "MySMeaS43" ",S43")
            self.Inst.write("CALCulate1:PARameter:DEFine " "MySMeaS41" ",S41")
            self.Inst.write("CALCulate1:PARameter:DEFine " "MySMeaS23" ",S23")
            self.Inst.write("DISPlay:WINDow1:STATE ON")
            self.Inst.write("DISPlay:WINDow1:TRACe1:FEED " "MySMeaS21" "")
            self.Inst.write("DISPlay:WINDow1:TRACe2:FEED " "MySMeaS43" "")
            self.Inst.write("DISPlay:WINDow1:TRACe3:FEED " "MySMeaS41" "")
            self.Inst.write("DISPlay:WINDow1:TRACe4:FEED " "MySMeaS23" "")
            self.Inst.write("DISPlay:WINDow1:TITLe:STATe ON")
            self.Inst.write("DISPlay:ANNotation:FREQuency ON")
        elif trace == 41:
            self.Inst.write("CALCulate1:PARameter:DEFine " "MySMeaS41" ",S41")
            self.Inst.write("DISPlay:WINDow1:STATE ON")
            self.Inst.write("DISPlay:WINDow1:TRACe3:FEED " "MySMeaS41" "")
            self.Inst.write("DISPlay:WINDow1:TITLe:STATe ON")
            self.Inst.write("DISPlay:ANNotation:FREQuency ON")
        elif trace == 23:
            self.Inst.write("CALCulate1:PARameter:DEFine " "MySMeaS23" ",S23")
            self.Inst.write("DISPlay:WINDow1:STATE ON")
            self.Inst.write("DISPlay:WINDow1:TRACe4:FEED " "MySMeaS23" "")
            self.Inst.write("DISPlay:WINDow1:TITLe:STATe ON")
            self.Inst.write("DISPlay:ANNotation:FREQuency ON")
        elif trace == 43:
            self.Inst.write("CALCulate1:PARameter:DEFine " "MySMeaS43" ",S43")
            self.Inst.write("DISPlay:WINDow1:STATE ON")
            self.Inst.write("DISPlay:WINDow1:TRACe2:FEED " "MySMeaS43" "")
            self.Inst.write("DISPlay:WINDow1:TITLe:STATe ON")
            self.Inst.write("DISPlay:ANNotation:FREQuency ON")
        elif trace == 21:
            self.Inst.write("CALCulate1:PARameter:DEFine " "MySMeaS21" ",S21")
            self.Inst.write("DISPlay:WINDow1:STATE ON")
            self.Inst.write("DISPlay:WINDow1:TRACe1:FEED " "MySMeaS21" "")
            self.Inst.write("DISPlay:WINDow1:TITLe:STATe ON")
            self.Inst.write("DISPlay:ANNotation:FREQuency ON")

        elif trace == 24:
            self.Inst.write("CALCulate1:PARameter:DEFine " "MySMeaS21" ",S21")
            self.Inst.write("CALCulate1:PARameter:DEFine " "MySMeaS11" ",S11")
            self.Inst.write("CALCulate1:PARameter:DEFine " "MySMeaS12" ",S12")
            self.Inst.write("CALCulate1:PARameter:DEFine " "MySMeaS22" ",S22")
            self.Inst.write("DISPlay:WINDow1:STATE ON")
            self.Inst.write("DISPlay:WINDow1:TRACe1:FEED " "MySMeaS21" "")
            self.Inst.write("DISPlay:WINDow1:TRACe2:FEED " "MySMeaS11" "")
            self.Inst.write("DISPlay:WINDow1:TRACe3:FEED " "MySMeaS12" "")
            self.Inst.write("DISPlay:WINDow1:TRACe4:FEED " "MySMeaS22" "")
            self.Inst.write("DISPlay:WINDow1:TITLe:STATe ON")
            self.Inst.write("DISPlay:ANNotation:FREQuency ON")

        else:
            self.Inst.write("CALCulate1:PARameter:DEFine " "MySMeaS21" ",S21")
            self.Inst.write("DISPlay:WINDow1:STATE ON")
            self.Inst.write("DISPlay:WINDow1:TRACe1:FEED " "MySMeaS21" "")
            self.Inst.write("DISPlay:WINDow1:TITLe:STATe ON")
            self.Inst.write("DISPlay:ANNotation:FREQuency ON")
        self.Inst.write("TRIG:SCOPe ALL")
        self.Inst.write("SENSe:SWEep:MODE CONTinuous")

    def rst(self):
        self.Inst.write("*RST")
        pass  # zzz只要s21

    def get_allSet(self):
        strafre = self.Inst.query("SENSe:FREQuency:STARt?")
        Stopfre = self.Inst.query("SENSe:FREQuency:STOP?")
        points = self.Inst.query("SENSe:SWEep:POINts?")
        avg = self.Inst.query("SENSe:AVERage:COUNt?")
        Ifband = self.Inst.query("SENS:BAND:RES?")
        power = self.Inst.query("SOURce:POWer?")
        #        nomal_set={strafre,Stopfre,points,avg,Ifband,power}
        return strafre, Stopfre, points, avg, Ifband, power

    def _setAverage(self, N):
        N = int(N)
        if N > 1:
            self.instrhandle.write("SENS:AVER ON")
            self.instrhandle.write("SENS:AVER:CLE")
            self.instrhandle.write("SENS:AVER:MODE AUTO")
            self.instrhandle.write("SENSe:AVERage:COUNt {}".format(N))

    def set_startstopFre(self, strafre, stopfre):
        self.Inst.write("SENSe:FREQuency:STARt {}".format(float(strafre) * 1e9))
        self.Inst.write("SENSe:FREQuency:STOP {}".format(float(stopfre) * 1e9))  # Ghz

    #        self.Inst.write('INITiate:IMMediate;*wai')
    #        self.Inst.write('Display:WINDow1:TRACe1:Y:Scale:AUTO')
    def set_power(self, power):
        self.Inst.write("SOURce:POWer {}".format(float(power)))

    def set_points_band(self, points, ifband):
        self.Inst.write("SENSe:SWE:POIN {}".format(int(points)))
        self.Inst.write("SENSe:BAND:RES {}".format(float(ifband)))

    def set_allSetting(self, strafre, stopfre, points, avg, ifband, power):
        self.Inst.write("SENSe:FREQuency:STARt {}".format(float(strafre) * 1e9))
        self.Inst.write("SENSe:FREQuency:STOP {}".format(float(stopfre) * 1e9))  # Ghz
        self.Inst.write("SENSe:SWE:POIN {}".format(int(points)))
        self.Inst.write("SENSe:AVERage:COUNt {}".format(int(avg)))
        self.Inst.write("SENSe:BAND:RES {}".format(float(ifband)))
        #        self.Inst.write('SOURce:POWer{}'.format(power))
        self.Inst.write("SOURce:POWer {}".format(float(power)))
        print("write data")
        print(strafre, stopfre, points, avg, ifband, power)

    def testwrite_comand(self, string):
        return self.Inst.write(string)

    def testquery_comand(self, string):
        return self.Inst.query(string)

    def get_data(self):
        if self.trace == 21:
            return self.get_data_S21()
        elif self.trace == 23:
            return self.get_data_S23()
        elif self.trace == 41:
            return self.get_data_S41()
        elif self.trace == 43:
            return self.get_data_S43()
        elif self.trace == 4:
            return self.get_data_all()

        raise Exception("Unknown trace: ", self.trace)

        self.Inst.write("FORMat:DATA REAL,64")
        self.Inst.write("FORMat:DATA ASCII")
        self.Inst.write("CALCulate:PARameter:SELect " "MySMeaS21" "")
        self.Inst.write("INITiate:IMMediate;*wai")
        self.Inst.write("CALCulate:DATA? SDATA")
        sawdata = self.Inst.read()
        sawdata1 = sawdata.split(",")
        Sdata = []
        for i in range(len(sawdata1)):
            if i % 2 == 0:
                Sdata.append(complex(float(sawdata1[i]), float(sawdata1[i + 1])))
        Sdata = np.array(Sdata)
        return Sdata
        # print(data,f)

    def get_data_S21(self):
        self.Inst.write("FORMat:DATA REAL,64")
        self.Inst.write("FORMat:DATA ASCII")
        self.Inst.write("INITiate:IMMediate;*wai")
        self.Inst.write("CALCulate:PARameter:SELect " "MySMeaS21" "")
        self.Inst.write("CALCulate:DATA? SDATA")
        sawdata = self.Inst.read()
        sawdata1 = sawdata.split(",")
        Sdata = []
        for i in range(len(sawdata1)):
            if i % 2 == 0:
                Sdata.append(complex(float(sawdata1[i]), float(sawdata1[i + 1])))
        Sdata = np.array(Sdata)
        return Sdata
        # print(data,f)

    def get_data_S43(self):
        self.Inst.write("FORMat:DATA REAL,64")
        self.Inst.write("FORMat:DATA ASCII")
        self.Inst.write("INITiate:IMMediate;*wai")
        self.Inst.write("CALCulate:PARameter:SELect " "MySMeaS43" "")
        self.Inst.write("CALCulate:DATA? SDATA")
        sawdata = self.Inst.read()
        sawdata1 = sawdata.split(",")
        Sdata = []
        for i in range(len(sawdata1)):
            if i % 2 == 0:
                Sdata.append(complex(float(sawdata1[i]), float(sawdata1[i + 1])))
        Sdata = np.array(Sdata)
        return Sdata
        # print(data,f)

    def get_data_S41(self):
        self.Inst.write("FORMat:DATA REAL,64")
        self.Inst.write("FORMat:DATA ASCII")
        self.Inst.write("INITiate:IMMediate;*wai")
        self.Inst.write("CALCulate:PARameter:SELect " "MySMeaS41" "")
        self.Inst.write("CALCulate:DATA? SDATA")
        sawdata = self.Inst.read()
        sawdata1 = sawdata.split(",")
        Sdata = []
        for i in range(len(sawdata1)):
            if i % 2 == 0:
                Sdata.append(complex(float(sawdata1[i]), float(sawdata1[i + 1])))
        Sdata = np.array(Sdata)
        return Sdata
        # print(data,f)

    def get_data_S23(self):
        self.Inst.write("FORMat:DATA REAL,64")
        self.Inst.write("FORMat:DATA ASCII")
        self.Inst.write("INITiate:IMMediate;*wai")
        self.Inst.write("CALCulate:PARameter:SELect " "MySMeaS23" "")
        self.Inst.write("CALCulate:DATA? SDATA")
        sawdata = self.Inst.read()
        sawdata1 = sawdata.split(",")
        Sdata = []
        for i in range(len(sawdata1)):
            if i % 2 == 0:
                Sdata.append(complex(float(sawdata1[i]), float(sawdata1[i + 1])))
        Sdata = np.array(Sdata)
        return Sdata
        # print(data,f)

    def fixdata(self, sawdata):
        sawdata1 = sawdata.split(",")
        Sdata = []
        for i in range(len(sawdata1)):
            if i % 2 == 0:
                Sdata.append(complex(float(sawdata1[i]), float(sawdata1[i + 1])))
        Sdata = np.array(Sdata)
        return Sdata

    def get_data_all(self):
        self.Inst.write("FORMat:DATA REAL,64")
        self.Inst.write("FORMat:DATA ASCII")
        self.Inst.write("INITiate:IMMediate;*wai")
        self.Inst.write("CALCulate:PARameter:SELect " "MySMeaS21" "")
        sawdata21 = self.Inst.query("CALCulate:DATA? SDATA")
        self.Inst.write("CALCulate:PARameter:SELect " "MySMeaS43" "")
        sawdata43 = self.Inst.query("CALCulate:DATA? SDATA")
        self.Inst.write("CALCulate:PARameter:SELect " "MySMeaS41" "")
        sawdata41 = self.Inst.query("CALCulate:DATA? SDATA")
        self.Inst.write("CALCulate:PARameter:SELect " "MySMeaS23" "")
        sawdata23 = self.Inst.query("CALCulate:DATA? SDATA")
        Sdata21 = self.fixdata(sawdata21)
        Sdata43 = self.fixdata(sawdata43)
        Sdata41 = self.fixdata(sawdata41)
        Sdata23 = self.fixdata(sawdata23)
        return Sdata21, Sdata43, Sdata41, Sdata23

    def get_data_all_2port(self):
        self.Inst.write("FORMat:DATA REAL,64")
        self.Inst.write("FORMat:DATA ASCII")
        self.Inst.write("INITiate:IMMediate;*wai")
        self.Inst.write("CALCulate:PARameter:SELect " "MySMeaS21" "")
        sawdata21 = self.Inst.query("CALCulate:DATA? SDATA")
        self.Inst.write("CALCulate:PARameter:SELect " "MySMeaS11" "")
        sawdata11 = self.Inst.query("CALCulate:DATA? SDATA")
        self.Inst.write("CALCulate:PARameter:SELect " "MySMeaS12" "")
        sawdata12 = self.Inst.query("CALCulate:DATA? SDATA")
        self.Inst.write("CALCulate:PARameter:SELect " "MySMeaS22" "")
        sawdata22 = self.Inst.query("CALCulate:DATA? SDATA")
        Sdata21 = self.fixdata(sawdata21)
        Sdata11 = self.fixdata(sawdata11)
        Sdata12 = self.fixdata(sawdata12)
        Sdata22 = self.fixdata(sawdata22)
        return Sdata21, Sdata11, Sdata12, Sdata22

    def Scan_freSpectrum(self, powers):
        S_data = list()
        for i in range(len(powers)):
            self.set_power(powers[i])
            f, data = self.get_data()
            S_data.append(data)
            if i % (int(len(powers) / 10)) == 0:
                print("scan_Fr: {:.0f}% completed.".format((i + 1) / len(powers) * 100))
        return f, np.array(S_data)

    def close(self):
        self.Inst.close()
        pass

    def two_tone_vna(self):
        self.Inst.write(
            "TRIGger:SEQuence:SOURce MANual"
        )  #'TRIGger:SEQuence:SOURce ''External'''
        self.Inst.write(
            "TRIG:SCOP CURRENT"
        )  # TRIG:SCOP CURRENT   #TRIGger:SEQuence:SCOPe ''Chan''
        self.Inst.write("SENSe:SWEep:TRIGger:DELay {}".format(float(1) * 1e-3))
        self.Inst.write("TRIG:CHAN:AUX2 0")
        self.Inst.write("TRIG:CHAN:AUX1 1")
        self.Inst.write("TRIG:CHAN:AUX1:INTerval POINt")
        pass

    def two_tone_vna2(self):
        self.Inst.write(
            "TRIGger:SEQuence:SOURce MANual"
        )  #'TRIGger:SEQuence:SOURce ''External'''
        self.Inst.write(
            "TRIG:SCOP CURRENT"
        )  # TRIG:SCOP CURRENT   #TRIGger:SEQuence:SCOPe ''Chan''
        self.Inst.write("SENSe:SWEep:TRIGger:DELay {}".format(float(1) * 1e-3))
        self.Inst.write("TRIG:CHAN:AUX1 0")
        self.Inst.write("TRIG:CHAN:AUX2 1")
        self.Inst.write("TRIG:CHAN:AUX2:INTerval POINt")
        pass

    def poweroff(self, i):
        self.Inst.write("SOURce:POWer{}:MODE OFF".format(int(i)))

    def SingleSweep(self):
        self.Inst.write("SENSe:SWEep:MODE Single")
