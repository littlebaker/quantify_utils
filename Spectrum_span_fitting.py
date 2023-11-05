#%%
from hunanu.Instrument.VNA_keysightN5231b_4Portsand2Ports import keysight_vna
from hunanu.Instrument.GS200_SPT  import GS_200
from hunanu.Instrument.SA_keysightN9020B import keysightSA_N9020B as SA
from tqdm import tqdm
from pprint import pprint, pformat
# from tqmd import tqmd
from fit_show.T1_T2_fit_show import Image_plot
from fit_show.Vna_phase_unwarp import handlePhase_2D ,handlePhase
import time
import os
import numpy as np
import hvplot.xarray
import xarray as xr
import holoviews as hv
renderer = hv.renderer('matplotlib')
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator, MaxNLocator,LogLocator,FixedFormatter
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)

from hunanu.Instrument.SMB_100A import MW
from def_plotly import heatmap
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator, MaxNLocator,LogLocator,FixedFormatter
from hunanu.Functions.function0 import savetxt,plot_removebg,todB,imshow,phase_handle,Lorentz_fitting_VNA,half_value_width
from scipy.optimize import leastsq
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
from hunanu.Functions.function_ying import heatmap_crange_cmap,Plot1D,Plot2D,dir_name,CheckBox,Record_EXP_CONT,Energy_spectrum_scan_2D
from scipy.optimize import leastsq
def todB(S21):
    dBS21=10*np.log10(np.abs(S21))
    return(dBS21)
def NorAmp(S21):
    Amp21=(np.abs(S21))/np.max(np.abs(S21))
    return(Amp21**2)
def Phase_handle(Phase,f1,phase_delay,phase_offset):
    phase_correct=np.zeros((Phase.shape[0]))
    for i in range(len(f1)):
        phase_temp=Phase[i]+2*np.pi*f1[i]*phase_delay+phase_offset
        phase_correct[i]=np.angle(np.cos(phase_temp)+1j*np.sin(phase_temp))
    return phase_correct
def deldata(data,index):
    newdata=[data[i] for i in range(len(data)) if i!=index]
    return newdata
from scipy import optimize
def lorentz(x,A0,A,x0,kappa):
    y=A0-1*np.abs(A*kappa)/(4*(x-x0)**2+kappa**2)
    return y
def findFWHM(xData,yData):
    sort=np.argsort(yData)
    x0=sort[0]
    x0_data=xData[x0]
    ymin=np.min(yData)
    hight=np.max(yData)-np.min(yData)
    halfhight=ymin+hight/2
    diff=yData-halfhight
    Left=diff[0:x0]
    Right=diff[x0:-1]
    if x0==0:
        x2=np.argsort(Right)[0]+x0
        kappa=2*x2
    else:
        x1=np.argsort(Left)[0]
        x2=np.argsort(Right)[0]+x0
        kappa=xData[x2]-xData[x1]
    return x0_data,kappa
def fitlorentz(xData,yData):
    x00,kappa0=findFWHM(xData,yData)
    p0=[np.max(yData),(np.min(yData)-np.max(yData))*kappa0,x00,kappa0]
    A0,A,x0,kappa=optimize.curve_fit(lorentz,xData,yData,p0=p0)[0]
    return A0,A,x0,kappa

def Energy_spectrum_scan_fitting(RTatt,vna_name, vna_trace,gs, gs_step, elec_delay,
    volt,fstart,fstop,points,ifband,vnapower,dirname,name=None):

    about={
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
        "DC": {
            "current": f" {volt[0]}mA",
            "step": gs_step
        },
    }
    pprint(about)
    if gs!=None:
        current_time = np.abs(volt[0]-round(float(gs.Inst.query(':SOUR:LEV?'))*10e2,3))/gs_step
    else:
        current_time=0
    pprint(f"电流表偏置到初始值需要{current_time}秒({round(current_time/60)}分钟)")
    timenow = time.strftime("%y%m%d%H%M")
    timeused = np.abs(round((points/ifband+(volt[-1]-volt[0])/(len(volt))/gs_step)*(len(volt)+1),1)) + current_time
    if not CheckBox(timeused):
        return None, None, None
    if name==None:
        name = f"Energy Spectrum freq {fstart} to {fstop} GHz vnapower {vnapower} {timenow}"
    else:
        name=name
    print("\n")
    Record_EXP_CONT(dirname,name)
    vna = keysight_vna(vna_name,trace=vna_trace) 
    f = np.linspace(fstop,fstart,points)
    f1 = f[::-1]  
    vna.set_power(vnapower)
    vna.set_startstopFre(fstart,fstop)
    vna.set_points_band(points,ifband)
    
    if gs!=None:
        gs.setlevel_slowV2(volt[0],gs_step)
    
    
    Sdata43=vna.get_data()
    S_dataplot_amp=20*np.log10(np.abs(Sdata43))

    Phase=np.angle(np.array(Sdata43))
    S_dataplot_phase=phase_handle(Phase[None,:],f1,elec_delay)[0]

    Amp=(S_dataplot_amp)
    A0,A,x0,kappa=fitlorentz(f1[:],Amp[:])
    fig,axes=plt.subplots(1,2,figsize = (10,4))
    axes[0].plot(f1, Amp)
    axes[0].plot(f1[:],lorentz(f1[:],A0,A,x0,kappa),'r--')
    axes[0].set_title(rf"Amp $\kappa$={round(kappa*1000,3)}MHz")
    axes[0].set_xlabel("Freqency (GHz)")
    axes[0].set_ylabel(r"Amp (dB)")
    axes[1].plot(f1, S_dataplot_phase.T)
    axes[1].set_title("Phase")
    axes[1].set_xlabel("Freqency (GHz)")
    axes[1].set_ylabel("Degree")
    plt.savefig(dirname+name+'.png',format='png',dpi=600)
    plt.savefig(dirname+name+'.pdf',format='pdf',dpi=600)
    plt.show()

    about_txt=open(dirname+name+'.txt','a')
    about_txt.write(pformat(about))
    about_txt.close()
    np.savez(dirname+name,frequency=f1,volt=volt,Sdata=np.array(Sdata43),vnapower=vnapower,about=about)
    fig_amp = (Plot1D(f1,S_dataplot_amp,labels=["Frequency(GHz)","Amplitude(dB)"],save_data=dirname+name+"_amp",save_fig=dirname+name+"_amp")*Plot1D(f1,lorentz(f1[:],A0,A,x0,kappa),labels=["Frequency(GHz)","Amplitude(dB)"],save_data=dirname+name+"_amp",save_fig=dirname+name+"_amp"))
    fig_phase = Plot1D(f1,S_dataplot_phase,labels=["Frequency(GHz)","Phase(deg)"],save_data=dirname+name+"_phase",save_fig=dirname+name+"_phase")
    return fig_amp, fig_phase


#%%


dirname = dir_name(f'./data/CEL_V1/Cavity_spectrum1D/{time.strftime("%Y%m%d")}/')
fig1,fig2 = Energy_spectrum_scan_fitting(
    name="Right_cavity_spectrum_scan_and_fitting_{}".format(time.strftime("%m%d%H%M")),
    RTatt=60,
    vna_name = "vna2", vna_trace=43, elec_delay=49.99,
    gs=None, gs_step=0.001,
    volt=[0],
    fstart=6.43, fstop=6.47, points=201,
    ifband=20, vnapower=-15,
    dirname=dirname
)
# %%
f=np.linspace(1,10,10)
S_dataplot_amp=np.linspace(1,20,10)
fig1 = Plot1D(f,S_dataplot_amp,labels=["Frequency(GHz)","Amplitude(dB)"],save_data=dirname+"name"+"_amp",save_fig=dirname+"name"+"_amp")
# %%
