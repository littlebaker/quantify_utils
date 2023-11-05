# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 07:06:55 2022

@author: DELL
"""
#%%
from hunanu.Instrument.VNA_keysightN5231b_4Portsand2Ports import keysight_vna
from hunanu.Instrument.GS200_SPT  import GS_200
from hunanu.Instrument.SA_keysightN9020B import keysightSA_N9020B as SA
from tqdm import tqdm
# from tqmd import tqmd
from fit_show.T1_T2_fit_show import Image_plot
from fit_show.Vna_phase_unwarp import handlePhase_2D ,handlePhase
import time
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

filename = f'/home/machine1/E/Data/LSY_new/data/CELv1/'
if not os.path.exists(filename):
    os.makedirs(filename)
from hunanu.Instrument.SMB_100A import MW
from def_plotly import heatmap
matplotlib.rcParams['text.usetex'] = True
# matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['font.family']='sans-serif'
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator, MaxNLocator,LogLocator,FixedFormatter
from hunanu.Functions.function0 import savetxt,plot_removebg,todB,imshow,phase_handle,Lorentz_fitting_VNA,half_value_width
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
#%%

def GS(num,setvalue,velocity):
    if velocity==0:
        GS_200('DC{}'.format(num)).Start_OutPut()
        print("GS_200_DC{} connected, work at {}mA".format(num,round(float(GS_200('DC{}'.format(num)).Inst.query(':SOUR:LEV?'))*10e2,3)))
    else:
        GS_200('DC{}'.format(num)).setlevel_slow(setvalue,velocity) 
        print("GS_200_DC{} has been set to {}mA".format(num,round(float(GS_200('DC{}'.format(num)).Inst.query(':SOUR:LEV?'))*10e2,3)))
        
GS(3,-0.245,0.001)

#%% DC-control-1
 
# gs_qubit_local=GS_200('DC2') 
#gs3=GS_200('DC3') 
gs4=GS_200('DC4') 
gs4.setCURRmode()
gs4.Start_OutPut()
gs4.getlevel()
# gs_gmon_local=GS_200('DC1') 
# vna=keysight_vna('vna1',trace=21)
#%%
# gs_qubit_local.Start_OutPut()
gs_global.Start_OutPut()
# gs_gmon_local.Start_OutPut()
#%% DC-control 2
# gs_qubit_local.setlevel_slow(0.00,0.001)#QUBIT4
gs_global.setlevel_slowV2(-0.556,0.0001)#QUBIT4
# gs_gmon_local.setlevel_slow(-0.57,0.002)#QUBIT4

#%% VNA-sweep-setup-1port
which='S21'
# whichtrace=int(which[1]+which[2])
cavityfre=6.745
measurewidth=0.01

vna=keysight_vna('vna1',trace=21)
which='S21'
RTatt=70
# fstart=cavityfre - measurewidth/2
# fstop=cavityfre + measurewidth/2

fstart=6
fstop=8

points=2001
ifband=10
vnapwr=-20
f=np.linspace(fstart,fstop,points)

vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
exec('test=vna.get_data_'+which+'()')
plt.figure('test')
plt.plot(f,todB(test))    
plt.xlabel('Frequency(Ghz)')
plt.ylabel('Amplitude(dB)')
plt.title('qubit')
np.savez(filename+'qubit_freq'+str(time.strftime("%H%M%S"))+'.npz',fre=f,vnapwr=vnapwr,Amp=todB(test),RTatt=RTatt)
plt.savefig(filename+'qubit_freq'+str(time.strftime("%H%M%S"))+".eps",format="eps",dpi=300)#%% VNA-sweep-setup-1port
plt.savefig(filename+'qubit_freq'+str(time.strftime("%H%M%S"))+".jpg",format="jpg",dpi=300)#%% VNA-sweep-setup-1port
#%%
which='S21'
# whichtrace=int(which[1]+which[2])
cavityfre=6.772
measurewidth=0.06

vna=keysight_vna('vna2',trace=43)
which='S21'
RTatt=60
fstart=cavityfre - measurewidth/2
fstop=cavityfre + measurewidth/2

# fstart=6
# fstop=8
Data_ALL = []
vnapwr_list = np.linspace(-25,-5,21)
for i in range(len(vnapwr_list)):
    points=601
    ifband=30
    vnapwr=vnapwr_list[i]
    f=np.linspace(fstart,fstop,points)
    
    vna.set_power(vnapwr)
    vna.set_startstopFre(fstart,fstop)
    vna.set_points_band(points,ifband)
    
    Sdata43=vna.get_data() 
    Data_ALL.append(Sdata43)
    S_dataplot_amp=20*np.log10(np.abs(Sdata43))
    
    Phase=np.angle(np.array(Sdata43))
    S_dataplot_phase=phase_handle(Phase[None,:],f,67.12)
    
    Amp=(S_dataplot_amp)
    A0,A,x0,kappa=fitlorentz(f[:],Amp[:])
    fig,axes=plt.subplots(1,2,figsize = (10,4))
    axes[0].plot(f, Amp)
    axes[0].plot(f[:],lorentz(f[:],A0,A,x0,kappa),'r--')
    axes[0].set_title(rf"Amp $\kappa$={round(kappa*1000,3)}MHz")
    axes[0].set_xlabel("Freqency (GHz)")
    axes[0].set_ylabel(r"Amp (dB)")
    axes[1].plot(f, S_dataplot_phase.T)
    axes[1].set_title("Phase")
    axes[1].set_xlabel("Freqency (GHz)")
    axes[1].set_ylabel("Degree")
    
    import os
    filename = f'/home/machine1/E/Data/Hedong/tanmo_cavity/'
    if not os.path.exists(filename):
        os.makedirs(filename)
    plt.savefig(filename+f'cavity_spectrum_vnapwr_{vnapwr}dB'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    
    Sdata=np.array(Sdata43)
    about=f'''
    BOX?; setup: in:port 2,out: 1;
    Attenuation: Vna={RTatt-vnapwr} dB={RTatt}+{-vnapwr} dB, inline = 60 dB = 42+8+10 dB,
    Amp: outline=76 dB(RT)+36 dB(4K) 
    
    VNA Paras: fstart = {fstart} GHz,fstop = {fstop} GHz,points = {points},ifband = {ifband} Hz,vnapwr = {vnapwr} dB 
    '''
    
    
    # filename1=filename+filenameafter+'.npz'
    np.savez(filename+f'cavity_spectrum_vnapwr_{vnapwr}dB'+str(time.strftime("%H%M%S"))+'.png',cavityfre=cavityfre,Sdata=Sdata,vnapwr=vnapwr,about=about)
    about_txt=open(filename+f'cavity_spectrum_vnapwr_{vnapwr}dB'+str(time.strftime("%H%M%S"))+'.txt','a')
    about_txt.write(str(about))
    about_txt.close()
np.savez(filename+'cavity_spectrum_Data_all'+str(time.strftime("%H%M%S"))+'.png',cavityfre=cavityfre,Data_ALL=Data_ALL,vnapwr_list=vnapwr_list)

#%%sweep setting
qubitnum=4
# gs1.setlevel_slow(0,0.01)#QUBIT2
# gs7.setlevel_slow(-1.648,0.01)#QUBIT3
# gs6.setlevel_slow(-0.75,0.01)#QUBIT4
# gs9.setlevel_slow(-1,0.01)#QUBIT5
# gs3.setlevel_slow(-0.25,0.01)#QUBIT6

filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_230609A/Energyspectrum/'
if not os.path.exists(filename):
    os.makedirs(filename)
    
qubitfre=6                                    #
qubit_bandwith=4                              #
vna=keysight_vna('vna2',trace=21)
which='S21'

Currentvalue_list = np.linspace(0,0.56,1)     #
gs_global.setlevel_slow(Currentvalue_list[0],0.0002)
#VNA_setting
RTatt=70
fstart=qubitfre - qubit_bandwith/2
fstop=qubitfre + qubit_bandwith/2
points=2001
ifband=10
vnapwr=-5
f=np.linspace(fstart,fstop,points)
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)



Line_width = []
Fit_condition = False
for i in tqdm(range(len(Currentvalue_list))):

    Currentvalue=Currentvalue_list[i]
    gs_global.setlevel_slow(Currentvalue,0.0001)
    
    # whichtrace=int(which[1]+which[2])
    exec('test=vna.get_data_'+which+'()')
    if Fit_condition==True:
        fg1 =plt.figure('test')
        origindata = todB(test)
        fitdata = Lorentz_fitting_VNA(f,origindata)
        Width = round(half_value_width(f, fitdata)*1e3,3)
        plt.plot(f,origindata)
        plt.plot(f,fitdata)
        plt.xlabel('Frequency(Ghz)')
        plt.ylabel('Amplitude(dB)')
        plt.title('VNA Result_Width{}MHz'.format(Width))
    
        plt.savefig(filename+f'I_{Currentvalue}mA_Width_{Width}MHz_{time.strftime("%m%d%H%M%S")}'+'.png',format='png',dpi=600)
    else:
        fg1 =plt.figure('test')
        plt.plot(f,todB(test))    
        plt.xlabel('Frequency(Ghz)')
        plt.ylabel('Amplitude(dB)')
        plt.title('VNA Result')
        plt.savefig(filename+f'I_{Currentvalue}mA_{time.strftime("%m%d%H%M%S")}'+'.png',format='png',dpi=600)
        

#%%
filename = f'E:\\Data\\YingHu\\tanmo\\20230114\\'
if not os.path.exists(filename):
    os.makedirs(filename)
vna=keysight_vna('vna1',trace=21)
fstart = 10.8  # 4.066#4.251#6.476#6.4668
fstop = 11.2  # 4.066#4.251#6.476#6.4668
points = 200
ifband = 10  # 每秒扫描次数
vnapwr = -30
f = np.linspace(fstop, fstart, points)
f1 = f[::-1]
vna.set_power(vnapwr)
vna.set_startstopFre(fstart, fstop)
vna.set_points_band(points, ifband)
data = vna.get_data_S43()
plt.figure()
plt.plot(f1, 20 * np.log10(np.abs(data)))
about = ' BOX:A; Reflection:setup: in:port31,out:20; \n Attenuation: Vna=20+30dB=50dB,\
    inline=60dB=42+8(line)+20dB(add), outline=72dB(RT)+36dB(4K) \n\
        fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n \
        '.format(fstart, fstop, points, ifband, vnapwr)
#%%sweep setting
qubitnum=2

qubitfre=5.71
qubit_bandwith=0.5
Currentvalue=-0.56



for i in range(1):
    
    Currentvalue=Currentvalue+i*0.05
    print(Currentvalue)
    gs3.setlevel_slow(Currentvalue,0.005)
    
    vna=keysight_vna('vna1',trace=21)
    which='S21'
    # whichtrace=int(which[1]+which[2])
    RTatt=60
    fstart=5.71-0.03
    fstop=5.71+0.03
    points=401
    ifband=5
    vnapwr=-13
    f=np.linspace(fstart,fstop,points)
    
    vna.set_power(vnapwr)
    vna.set_startstopFre(fstart,fstop)
    vna.set_points_band(points,ifband)
    exec('test=vna.get_data_'+which+'()')
    plt.figure('test')
    plt.plot(f,todB(test))    
    plt.xlabel('Frequency(Ghz)')
    plt.ylabel('Amplitude(dB)')
    plt.title('VNA Result')
#%%sweep qubit energy specturm_global
qubitnum = 2
vna = keysight_vna('vna2',trace=21)
which = 'S21'
RTatt = 70

fstart, fstop, points = 6.7, 7.1, 201 

ifband = 30
vnapwr = -15

f = np.linspace(fstop,fstart,points)
f1 = f[::-1]  
  
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_2309D_XS02_231010/Energyspectrum/{time.strftime("%m%d")}_{time.strftime("%H%M%S")}_freq_{fstart}GHz_to_{fstop}GHz/'
if not os.path.exists(filename):
    os.makedirs(filename)

Globalnum = 4
filenameafter = 'Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs_global_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))
filename1=filename+filenameafter

volt = np.linspace(0.5, -0.5, 201)
gs_global.Start_OutPut()
gs_global.setlevel_slow(volt[0],0.005)
time1=time.time()
S_data=[]
k=0
for i in tqdm(range(len(volt))):
    gs_global.setlevel_slow(volt[i],0.005)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    
    Phase = np.angle(np.array(S_data))
    S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 48.13))
    
    k += 1
    if k==10:
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_amp, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Amp', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_phase, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Phase', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        k=0

# gs3.setlevel_slow(-0,0.01) 
Sdata=np.array(S_data)

# filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

about=' BOX?; setup: in:port 18,out: 20; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])

filename1=filename+filenameafter+'.npz'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)

heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_AMP",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(volt,f,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_PHASE",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
#%%sweep qubit energy specturm_coupler_local
qubitnum = 2
vna = keysight_vna('vna2',trace=21)
which = 'S21'
RTatt = 70

fstart, fstop, points = 6.5, 7.1, 301 

ifband = 30
vnapwr = -15

f = np.linspace(fstop,fstart,points)
f1 = f[::-1]  
  
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_2309D_XS02_231010/Energyspectrum_Coupler_local/{time.strftime("%m%d")}_{time.strftime("%H%M%S")}_freq_{fstart}GHz_to_{fstop}GHz/'
if not os.path.exists(filename):
    os.makedirs(filename)

Globalnum = 4
filenameafter = 'Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_GS_coupler_local_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))
filename1=filename+filenameafter

volt = np.linspace(-0.86, -0.36, 51)
gs_gmon_local.Start_OutPut()
gs_gmon_local.setlevel_slow(volt[0],0.005)
time1=time.time()
S_data=[]
k=0
for i in tqdm(range(len(volt))):
    gs_gmon_local.setlevel_slow(volt[i],0.005)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    
    Phase = np.angle(np.array(S_data))
    S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 62.06))
    
    k += 1
    if k==10:
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_amp, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Amp', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_phase, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Phase', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        k=0

# gs3.setlevel_slow(-0,0.01) 
Sdata=np.array(S_data)

# filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

about=' BOX?; setup: in:port 18,out: 20; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])

filename1=filename+filenameafter+'.npz'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)

heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_AMP",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(volt,f,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_PHASE",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
#%%sweep qubit energy specturm_qubit_local
qubitnum = 2
vna = keysight_vna('vna2',trace=21)
which = 'S21'
RTatt = 60

fstart, fstop, points =6.5, 7.1, 301 

ifband = 40
vnapwr = -20

f = np.linspace(fstop,fstart,points)
f1 = f[::-1]  
  
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_2309D_XS02_231010/Energyspectrum_qubit_local/{time.strftime("%m%d")}_{time.strftime("%H%M%S")}_freq_{fstart}GHz_to_{fstop}GHz/'
if not os.path.exists(filename):
    os.makedirs(filename)

Globalnum = 4
filenameafter = 'Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_GS_qubit_local_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))
filename1=filename+filenameafter

volt = np.linspace(0, 0.5, 201)
gs_qubit_local.Start_OutPut()
gs_qubit_local.setlevel_slow(volt[0],0.0005)
time1=time.time()
S_data=[]
k=0
for i in tqdm(range(len(volt))):
    gs_qubit_local.setlevel_slow(volt[i],0.0005)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    
    Phase = np.angle(np.array(S_data))
    S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 62.06))
    
    k += 1
    if k==10:
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_amp, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Amp', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_phase, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Phase', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        k=0

# gs3.setlevel_slow(-0,0.01) 
Sdata=np.array(S_data)

# filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

about=' BOX?; setup: in:port 18,out: 20; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])

filename1=filename+filenameafter+'.npz'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)

heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_AMP",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(volt,f,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_PHASE",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
#%%sweep qubit energy specturm_qubit_local
gmon_current_list = np.linspace(-1 ,1,101)
for gmon_current in gmon_current_list:
    gs_gmon_local.setlevel_slow(gmon_current,0.005)#QUBIT4
    qubitnum = 2
    vna = keysight_vna('vna2',trace=43)
    which = 'S43'
    RTatt = 50
    
    fstart, fstop, points = 5.3, 5.45, 151 
    
    ifband = 5
    vnapwr = -20
    
    f = np.linspace(fstop,fstart,points)
    f1 = f[::-1]  
      
    vna.set_power(vnapwr)
    vna.set_startstopFre(fstart,fstop)
    vna.set_points_band(points,ifband)
    # exec('test=vna.get_data_'+which+'()')
    # plt.figure('test')
    # plt.plot(f[::-1],todB(test))    
    # plt.xlabel('Frequency(Ghz)')
    # plt.ylabel('Amplitude(dB)')   
    
    filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_230609D/Energyspectrum/{time.strftime("%m%d")}_freq_{fstart}GHz_to_{fstop}GHz_gmon/'
    if not os.path.exists(filename):
        os.makedirs(filename)
    
    Globalnum = 4
    filenameafter = 'Energyspetrum_Tran{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs_qubit_local_{}mA'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum,gmon_current)+str(time.strftime("%m%d%H%M%S"))
    filename1=filename+filenameafter
    
    volt = np.linspace(-0.3,0.3,61)
    gs_qubit_local.Start_OutPut()
    gs_qubit_local.setlevel_slow(volt[0],0.005)
    time1=time.time()
    S_data=[]
    k = 0
    for i in tqdm(range(len(volt))):
        gs_qubit_local.setlevel_slow(volt[i],0.002)
        data=vna.get_data() 
        S_data.append(data)
        S_dataplot=np.rot90(np.array(S_data))
        S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
        
        Phase = np.angle(np.array(S_data))
        S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 68.04))
    
        
        k=k+1
        if k==100:
            heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_amp, corlor_scale='bwr', xlabel='Current[mA]',
                    ylabel='Frequency[Ghz]', title='sweep_Amp', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
            heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_phase, corlor_scale='bwr', xlabel='Current[mA]',
                    ylabel='Frequency[Ghz]', title='sweep_Phase', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
            k=0
    
    # gs3.setlevel_slow(-0,0.01) 
    Sdata=np.array(S_data)
    
    filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))
    
    about=' BOX?; setup: in:port 18,out: 20; \n\
        Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
        fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
        current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])
    
    filename1=filename+filenameafter+'.npz'
    np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
    about_txt=open(filename+filenameafter+'.txt','a')
    about_txt.write(str(about))
    about_txt.close()
    
    filename2=filename+filenameafter
    imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
    plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)
    
    imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
    plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)
    
    heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_AMP",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
    heatmap(volt,f,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_PHASE",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
    
#%%sweep qubit energy specturm_qubit_local2

qubitnum = 2
vna = keysight_vna('vna2',trace=43)
which = 'S43'
RTatt = 0

fstart, fstop, points = 5.28, 5.43, 151 

ifband = 15
vnapwr = -20

f = np.linspace(fstop,fstart,points)
f1 = f[::-1]  
  
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_2309A_XS03/Energyspectrum/{time.strftime("%m%d")}_freq_{fstart}GHz_to_{fstop}GHz_gmon/'
if not os.path.exists(filename):
    os.makedirs(filename)

Globalnum = 4
filenameafter = 'Energyspetrum_Tran{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs_gmon_local_{}mA'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum,gmon_current)+str(time.strftime("%m%d%H%M%S"))
filename1=filename+filenameafter

# volt = np.linspace(-3,3,601)
# gs_qubit_local.Start_OutPut()
# gs_qubit_local.setlevel_slow(volt[0],0.005)
time1=time.time()
S_data=[]
k = 0

gs_qubit_local.setlevel_slow(-0.09,0.002)
gmon_current_list = np.linspace(3, -3, 601)
volt = gmon_current_list
for gmon_current in tqdm(gmon_current_list):
    
# for i in tqdm(range(len(volt))):
    gs_gmon_local.setlevel_slow(gmon_current,0.005)#QUBIT4
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    
    Phase = np.angle(np.array(S_data))
    S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 68.04))

    
    k=k+1
    if k==20:
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_amp, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Amp', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_phase, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Phase', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        k=0

# gs3.setlevel_slow(-0,0.01) 
Sdata=np.array(S_data)

filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

about=' BOX?; setup: in:port 18,out: 20; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])

filename1=filename+filenameafter+'.npz'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
plt.clf()
imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[-1],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)
plt.clf()
imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[-1],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase",color='bwr',interpolations='None') 
plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)
plt.clf()
heatmap(volt,f,S_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_AMP",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(volt,f,S_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_PHASE",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
    
    #%%sweep qubit energy specturm_qubit_local
gmon_current_list = np.linspace(-1 ,1,101)

qubitnum = 2
vna = keysight_vna('vna2',trace=21)
which = 'S21'
RTatt = 50

fstart, fstop, points = 6, 7.5, 151 

ifband = 20
vnapwr = -15

f = np.linspace(fstop,fstart,points)
f1 = f[::-1]  
# for gmon_current in gmon_current_list:
    # gs_gmon_local.setlevel_slow(gmon_current,0.005)#QUBIT4      
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_230609D/Energyspectrum_qubit_local/{time.strftime("%m%d")}_freq_{fstart}GHz_to_{fstop}GHz_Vs_gmon/'
if not os.path.exists(filename):
    os.makedirs(filename)

Globalnum = 4
filenameafter = 'Energyspetrum_{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs_qubit_local_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))
filename1=filename+filenameafter

volt = np.linspace(-1 ,1,101)
# gs_qubit_local.Start_OutPut()
# gs_qubit_local.setlevel_slow(volt[0],0.005)
# time1=time.time()
S_data=[]
k = 0

gmon_current_list = np.linspace(-1 ,1,101)

# for i in tqdm(range(len(volt))):
for gmon_current in tqdm(volt):
    gs_gmon_local.setlevel_slow(gmon_current,0.005)#QUBIT4      
    # gs_qubit_local.setlevel_slow(volt[i],0.002)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    
    Phase = np.angle(np.array(S_data))
    S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 68.04))

    
    k=k+1
    if k==100:
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_amp, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Amp', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_phase, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Phase', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        k=0

# gs3.setlevel_slow(-0,0.01) 
Sdata=np.array(S_data)

filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

about=' BOX?; setup: in:port 18,out: 20; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])

filename1=filename+filenameafter+'.npz'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)

heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_AMP",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(volt,f,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_PHASE",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
#%%sweep cavity VS power
qubitnum = 2
vna = keysight_vna('vna2',trace=21)
which = 'S21'
RTatt = 60

fstart, fstop, points =6.45-0.004, 6.45+0.004, 101 

ifband = 20
# vnapwr = -20

f = np.linspace(fstop,fstart,points)
f1 = f[::-1]  
  

vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   
filename = f'/home/machine1/E/Data/LSY_new/data/CELv1/'
if not os.path.exists(filename):
    os.makedirs(filename)

Globalnum = 4
filenameafter = 'Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_GS_qubit_local_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))
filename1=filename+filenameafter

vnapwr_list = np.linspace(-20, 6, 27)
# gs_qubit_local.Start_OutPut()
# gs_qubit_local.setlevel_slow(volt[0],0.0005)
time1=time.time()
S_data=[]
k=0
for i in tqdm(range(len(vnapwr_list))):
    vna.set_power(vnapwr_list[i])
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    
    Phase = np.angle(np.array(S_data))
    S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 62.06))
    
    # k += 1
    # if k==10:
    #     heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_amp, corlor_scale='bwr', xlabel='Current[mA]',
    #             ylabel='Frequency[Ghz]', title='sweep_Amp', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
    #     heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_phase, corlor_scale='bwr', xlabel='Current[mA]',
    #             ylabel='Frequency[Ghz]', title='sweep_Phase', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
    #     k=0

# gs3.setlevel_slow(-0,0.01) 
Sdata=np.array(S_data)

# filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

# about=' BOX?; setup: in:port 18,out: 20; \n\
#     Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
#     fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
#     current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])

# filename1=filename+filenameafter+'.npz'
# np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
# about_txt=open(filename+filenameafter+'.txt','a')
# about_txt.write(str(about))
# about_txt.close()

# filename2=filename+filenameafter
imshow(S_dataplot_amp,figname=filename+"_AMP",extent=[vnapwr_list[0],vnapwr_list[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
# plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
# plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

# imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
# plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
# plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)

heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_AMP",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(volt,f,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_PHASE",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
#%%sweep cavity energy specturm_local
qubitnum = 2
vna = keysight_vna('vna2',trace=21)
which = 'S21'
RTatt = 70

fstart, fstop, points = 6.84, 6.86, 101 

ifband = 40
vnapwr = -15

f = np.linspace(fstop,fstart,points)
f1 = f[::-1]  
  
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_2309C_XS01/AMP_VS_local_Current{time.strftime("%m%d")}_{time.strftime("%H%M%S")}_freq_{fstart}GHz_to_{fstop}GHz{time.strftime("%m%d")}_{time.strftime("%H%M%S")}_freq_{fstart}GHz_to_{fstop}GHz/'
if not os.path.exists(filename):
    os.makedirs(filename)

Globalnum = 4
filenameafter = 'AMP_VS_local_Current{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs_global_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))
filename1=filename+filenameafter

volt = np.linspace(-1.15, 1.15, 231)
gs_global.Start_OutPut()
gs_global.setlevel_slow(volt[0],0.005)
time1=time.time()
S_data=[]
k=0
for i in tqdm(range(len(volt))):
    gs_global.setlevel_slow(volt[i],0.005)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    
    Phase = np.angle(np.array(S_data))
    S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 48.13))
    
    k += 1
    if k==10:
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_amp, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Amp', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_phase, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Phase', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        k=0

# gs3.setlevel_slow(-0,0.01) 
Sdata=np.array(S_data)

# filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

about=' BOX?; setup: in:port 18,out: 20; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])

filename1=filename+filenameafter+'.npz'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)

heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_AMP",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(volt,f,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_PHASE",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)

#%%sweep AMP Vs current global 
qubitnum = 2
vna = keysight_vna('vna2',trace=4)
which = 'S21'
RTatt = 60

fstart, fstop, points = 6.45,6.45,11 

ifband = 30
vnapwr = -15

f = np.linspace(fstop,fstart,points)
f1 = f[::-1]  
  
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

filename = f'/home/machine1/E/Data/LSY_new/data/CELv1/cavity_VS_Global current/'
if not os.path.exists(filename):
    os.makedirs(filename)

Globalnum = 4
filenameafter='AMPVSCurrent{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)
filename1=filename+filenameafter

volt = np.linspace(0,1,1001)
gs_global.Start_OutPut()
gs_global.setlevel_slowV2(volt[0],0.0001)
time1=time.time()
S21_dataplot_amp_data=[]
S21_dataplot_phase_data=[]
S43_dataplot_amp_data=[]
S43_dataplot_phase_data=[]
k = 0
for i in tqdm(range(len(volt))):

    
    S21_data=[]
    gs_global.setlevel_slowV2(volt[i],0.0001)
    Sdata21,Sdata43,_,_=vna.get_data_all()  
    S21_data.append(Sdata21)
    S21_dataplot=np.rot90(np.array(S21_data))
    S21_dataplot_amp=20*np.log10(np.abs(S21_dataplot))
    
    Phase = np.angle(np.array(S21_data))
    S21_dataplot_phase = np.rot90(phase_handle(Phase, f1, 50.11))

    S21_dataplot_amp_data.append( np.mean(S21_dataplot_amp))
    S21_dataplot_phase_data.append( np.mean(S21_dataplot_phase))
    
    S43_data=[]
    S43_data.append(Sdata43)
    S43_dataplot=np.rot90(np.array(S43_data))
    S43_dataplot_amp=20*np.log10(np.abs(S43_dataplot))
    
    Phase = np.angle(np.array(S43_data))
    S43_dataplot_phase = np.rot90(phase_handle(Phase, f1, 50.1))

    S43_dataplot_amp_data.append( np.mean(S43_dataplot_amp))
    S43_dataplot_phase_data.append( np.mean(S43_dataplot_phase))
    
    
    k=k+1
    if k==200:
        
        # Phase = np.angle(np.rot90(np.expand_dims(np.array(S_data),axis=0)))
        # cp.expand_dims(loadfiltertxt('/home/machine1/g2/python_spt/Qlab/200-20-60.txt'),axis=0)
        # S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 68.04))
        
        
        
        plt.figure(figsize=(12, 8), dpi=80)
        plt.subplot(2,2,1)
        plt.plot(volt[0:i],S21_dataplot_amp_data[0:i])
        plt.title("S21amp")
        plt.xlabel("Current(mA)")
        plt.ylabel("Amplitude(dB)")   
        plt.subplot(2,2,2)
        plt.plot(volt[0:i],S21_dataplot_phase_data[0:i])
        plt.title("S21phase")
        plt.xlabel("Current(mA)")
        plt.ylabel("Phase")  
        
        plt.subplot(2,2,3)
        plt.plot(volt[0:i],S43_dataplot_amp_data[0:i])
        plt.title("S43amp")
        plt.xlabel("Current(mA)")
        plt.ylabel("Amplitude(dB)")   
        plt.subplot(2,2,4)
        plt.plot(volt[0:i],S43_dataplot_phase_data[0:i])
        plt.title("S43phase")
        plt.xlabel("Current(mA)")
        plt.ylabel("Phase")  
        
        filenameafter='AMPVSCurrent{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)
        filename2=filename+filenameafter
        plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
        plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)
        plt.clf()
        k=0

# gs3.setlevel_slow(-0,0.01) 
# Sdata=np.array(S_data)

# filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

about=' BOX?; setup: in:port 7,out: 13; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])

# filename1=filename+filenameafter+'.npz'
np.savez(filename2,fre=f,volt=volt,S21_dataplot_amp_data=S21_dataplot_amp_data, S21_dataplot_phase_data=S21_dataplot_phase_data, S43_dataplot_amp_data=S43_dataplot_amp_data, S43_dataplot_phase_data=S43_dataplot_phase_data, vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

# filename2=filename+filenameafter
# imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
# plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
# plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

# imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
# plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
# plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)

# heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_AMP",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
# heatmap(volt,f,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_PHASE",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
plt.figure(figsize=(12, 8), dpi=80)
plt.subplot(2,2,1)
plt.plot(volt[0:i],S21_dataplot_amp_data[0:i])
plt.title("S21amp")
plt.xlabel("Current(mA)")
plt.ylabel("Amplitude(dB)")   
plt.subplot(2,2,2)
plt.plot(volt[0:i],S21_dataplot_phase_data[0:i])
plt.title("S21phase")
plt.xlabel("Current(mA)")
plt.ylabel("Phase")  

plt.subplot(2,2,3)
plt.plot(volt[0:i],S43_dataplot_amp_data[0:i])
plt.title("S43amp")
plt.xlabel("Current(mA)")
plt.ylabel("Amplitude(dB)")   
plt.subplot(2,2,4)
plt.plot(volt[0:i],S43_dataplot_phase_data[0:i])
plt.title("S43phase")
plt.xlabel("Current(mA)")
plt.ylabel("Phase")  

filenameafter='AMPVSCurrent{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_2'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)
filename2=filename+filenameafter
plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)
plt.clf()
#%%sweep cavity VS power
qubitnum = 2
vna = keysight_vna('vna2',trace=21)
which = 'S21'
RTatt = 70

fstart, fstop, points = 6.848,6.852,201

ifband = 10
# vnapwr = -15

f = np.linspace(fstop,fstart,points)
f1 = f[::-1]  
  
# vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_2309C_XS01/Cavity mode_VS_Power{time.strftime("%m%d")}_{time.strftime("%H%M%S")}_freq_{fstart}GHz_to_{fstop}GHz/'
if not os.path.exists(filename):
    os.makedirs(filename)

Globalnum = 4
filenameafter=f'cavity mode VS Power_qubit_with_qubit_f_{fstart}_to_{fstop}_vnapwr_-28to5dB_RTatt{RTatt}'
filename1=filename+filenameafter

# volt = np.linspace(-1,1,201)
# gs_global.Start_OutPut()
# gs_global.setlevel_slow(volt[0],0.005)
time1=time.time()
S_data=[]
# S_dataplot_phase_data=[]
k = 0

vnapwrstart = -28
vnapwrstop = 4
vnapwr_points = 17
vnapwr =np.linspace(vnapwrstart,vnapwrstop,vnapwr_points)

for i in tqdm(range(len(vnapwr))):

    vna.set_power(vnapwr[i])
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    
    Phase = np.angle(np.array(S_data))
    S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 48.13))
    
    k += 1
    if k==10:
        heatmap(np.linspace(vnapwr[0], vnapwr[i], i), f, S_dataplot_amp, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Amp', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        heatmap(np.linspace(vnapwr[0], vnapwr[i], i), f, S_dataplot_phase, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Phase', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        k=0

# gs3.setlevel_slow(-0,0.01) 
# Sdata=np.array(S_data)

# filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

about=' BOX?; setup: in:port 7,out: 13; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz \n\
    vnapwr from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr[0],vnapwr[-1])

# filename1=filename+filenameafter+'.npz'
np.savez(filename2,fre=f,volt=volt,S_dataplot_amp_data=S_dataplot_amp_data, S_dataplot_phase_data=S_dataplot_phase_data, vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

# filename2=filename+filenameafter
# imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
# plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
# plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

# imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
# plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
# plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)

# heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_AMP",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
# heatmap(volt,f,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_PHASE",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
#%%sweep AMP Vs current local 
qubitnum = 2
vna = keysight_vna('vna1',trace=21)
which = 'S21'
RTatt = 50

fstart, fstop, points = 6.5688,6.5688,11 

ifband = 30
vnapwr = -15

f = np.linspace(fstop,fstart,points)
f1 = f[::-1]  
  
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

volt = np.linspace(-6,6,6001)
gs_qubit_local.Start_OutPut()
gs_gmon_local.setlevel_slow(0,0.005)
gs_qubit_local.setlevel_slow(volt[0],0.005)
S_dataplot_amp_data=[]
S_dataplot_phase_data=[]
k = 0
for i in tqdm(range(len(volt))):
    if i == 0:
        gmon_CURRENT_Value = round(float(gs_gmon_local.Inst.query(':SOUR:LEV?'))*10e2,3)
        filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_230609C/AMP_VS_Current_qubit_local/{time.strftime("%m%d")}_freq_{fstart}GHz_to_{fstop}GHz_gmon_{gmon_CURRENT_Value}mA/'
        if not os.path.exists(filename):
            os.makedirs(filename)
        filenameafter=f'AMPVSCurrent_gmon_{gmon_CURRENT_Value}mA_with_qubit_f_{fstart}_to_{fstop}_vnapwr{vnapwr}_RTatt{RTatt}_gs_qubit_local'
        filename1=filename+filenameafter
    
    S_data=[]
    gs_qubit_local.setlevel_slow(volt[i],0.0005)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    
    Phase = np.angle(np.array(S_data))
    S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 50))

    S_dataplot_amp_data.append( np.mean(S_dataplot_amp))
    S_dataplot_phase_data.append( np.mean(S_dataplot_phase))
    k=k+1
    if k==1000:
        
        # Phase = np.angle(np.rot90(np.expand_dims(np.array(S_data),axis=0)))
        # cp.expand_dims(loadfiltertxt('/home/machine1/g2/python_spt/Qlab/200-20-60.txt'),axis=0)
        # S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 68.04))
        
        
        
        plt.figure(figsize=(12, 4), dpi=100)
        plt.subplot(1,2,1)
        plt.plot(volt[0:i],S_dataplot_amp_data[0:i])
        plt.title("amp")
        plt.xlabel("Current(mA)")
        plt.ylabel("Amplitude(dB)")   
        plt.subplot(1,2,2)
        plt.plot(volt[0:i],S_dataplot_phase_data[0:i])
        plt.title("phase")
        plt.xlabel("Current(mA)")
        plt.ylabel("Phase")  
        
        filenameafter='AMPVSCurrent_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)
        filename2=filename+filenameafter
        plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
        plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)
        plt.clf()
        k=0

# gs3.setlevel_slow(-0,0.01) 
# Sdata=np.array(S_data)

# filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

about=' BOX?; setup: in:port 7,out: 13; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])

# filename1=filename+filenameafter+'.npz'
np.savez(filename2,fre=f,volt=volt,S_dataplot_amp_data=S_dataplot_amp_data, S_dataplot_phase_data=S_dataplot_phase_data, vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

# filename2=filename+filenameafter
# imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
# plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
# plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

# imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
# plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
# plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)

# heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_AMP",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
# heatmap(volt,f,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_PHASE",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)

#%%sweep qubit energy specturm_global
qubitnum = 2
vna = keysight_vna('vna1',trace=21)
which = 'S21'
RTatt = 50

fstart, fstop, points = 6.5688-0.005,6.5688+0.005,101 

ifband = 30
vnapwr = -20

f = np.linspace(fstop,fstart,points)
f1 = f[::-1]  
  
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_230609C/AMP_VS_Current_global/{time.strftime("%m%d")}_{time.strftime("%H%M%S")}_freq_{fstart}GHz_to_{fstop}GHz/'
if not os.path.exists(filename):
    os.makedirs(filename)

Globalnum = 4
filenameafter = 'AMP_VS_Current{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs_global_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))
filename1=filename+filenameafter

volt = np.linspace(-1, 1, 501)
gs_global.Start_OutPut()
gs_global.setlevel_slow(volt[0],0.005)
time1=time.time()
S_data=[]
k=0
for i in tqdm(range(len(volt))):
    gs_global.setlevel_slow(volt[i],0.001)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    
    Phase = np.angle(np.array(S_data))
    S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 60.12))
    
    k += 1
    if k==10:
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_amp, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Amp', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_phase, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Phase', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        k=0

# gs3.setlevel_slow(-0,0.01) 
Sdata=np.array(S_data)

# filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

about=' BOX?; setup: in:port 18,out: 20; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])

filename1=filename+filenameafter+'.npz'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)

heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_AMP",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(volt,f,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_PHASE",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)

#%%sweep qubit energy specturm_local
qubitnum = 2
vna = keysight_vna('vna1',trace=21)
which = 'S21'
RTatt = 50

fstart, fstop, points = 6.5691-0.001,6.5691+0.001,101 

ifband = 20
vnapwr = -20

f = np.linspace(fstop,fstart,points)
f1 = f[::-1]  
  
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_230609C/AMP_VS_Current_qubit_local_20230630/{time.strftime("%m%d")}_{time.strftime("%H%M%S")}_freq_{fstart}GHz_to_{fstop}GHz/'
if not os.path.exists(filename):
    os.makedirs(filename)

Globalnum = 4
filenameafter = 'AMP_VS_Current{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs_qubit_local_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))
filename1=filename+filenameafter

volt = np.linspace(-6, 6, 1201)
gs_qubit_local.Start_OutPut()
gs_qubit_local.setlevel_slow(volt[0],0.005)
time1=time.time()
S_data=[]
k=0
for i in tqdm(range(len(volt))):
    gs_qubit_local.setlevel_slow(volt[i],0.005)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    
    Phase = np.angle(np.array(S_data))
    S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 60.12))
    
    k += 1
    if k==100:
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_amp, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Amp', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_phase, corlor_scale='bwr', xlabel='Current[mA]',
                ylabel='Frequency[Ghz]', title='sweep_Phase', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
        k=0

# gs3.setlevel_slow(-0,0.01) 
Sdata=np.array(S_data)

# filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

about=' BOX?; setup: in:port 18,out: 20; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])

filename1=filename+filenameafter+'.npz'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)

heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_AMP",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(volt,f,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter+"_PHASE",filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
#%%
qubitnum=2
vna=keysight_vna('vna1',trace=21)
which='S21'
RTatt=60

fstart, fstop, points = 5.71-0.05, 5.71+0.05, 201 

ifband=20
vnapwr=-20

f=np.linspace(fstop,fstart,points)
f1=f[::-1]  
  
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

filename = f'/home/machine1/E/Data/YingHu/shaping2023/Energyspectrum/{time.strftime("%m%d")}_{time.strftime("%H%M%S")}_freq_{fstart}GHz_to_{fstop}GHz/'
if not os.path.exists(filename):
    os.makedirs(filename)

Globalnum=4
filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))
filename1=filename+filenameafter

volt=np.linspace(0.1,0.8,46)
gs7.Start_OutPut()
gs7.setlevel_slow(volt[0],0.005)
time1=time.time()
S_data=[]
for i in tqdm(range(len(volt))):
    gs7.setlevel_slow(volt[i],0.005)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    
    Phase = np.angle(np.array(S_data))
    S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 62.42))

    
    heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_amp, corlor_scale='bwr', xlabel='Current[mA]',
            ylabel='Frequency[Ghz]', title='sweep_Amp', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
    heatmap(np.linspace(volt[0], volt[i], i), f, S_dataplot_phase, corlor_scale='bwr', xlabel='Current[mA]',
            ylabel='Frequency[Ghz]', title='sweep_Phase', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)

    if i%(int(len(volt)/10))==0:
        usedtime=time.time()-time1
        time1=time.time()
        remaintime=(10-((i+1)*10/len(volt)))*usedtime
        print('scan_Fr: {:.0f}% completed.'.format(i/len(volt)*100))
        print('remain about {:.1f} mins.'.format(remaintime/60))
# gs3.setlevel_slow(-0,0.01) 
Sdata=np.array(S_data)

filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

about=' BOX?; setup: in:port 11,out: 13 ; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])

filename1=filename+filenameafter+'.npz'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+'.jpg',format='jpg',dpi=300)

imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+'.jpg',format='jpg',dpi=300)

heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
#%%单个比特效率
# qubitfre=4.671
# qubit_bandwith=0.3

# Currentvalue=0.36
# NoiseCuvalue=0.7



# def single_qubit_efficient()
# filename = f'/home/machine1/E/Data/YingHu/SinglePhoton/{time.strftime("%Y%m%d")}/'
filename = '/home/machine1/E/Data/YingHu/tanmo/Effciency&SmithCircles/'
if not os.path.exists(filename):
    os.makedirs(filename)   

gs_global.setCURRmode()
gs_global.getlevel()
gs_global.Start_OutPut()

qubitfre = 7.696
qubit_bandwith = 0.015

Currentvalue = 0.592
NoiseCuvalue = 0.9

vna = keysight_vna('vna1',trace=21)
which = 'S21'
# whichtrace=int(which[1]+which[2])
RTatt = 70
fstart = qubitfre - qubit_bandwith/2
fstop = qubitfre + qubit_bandwith/2
points = 201
ifband = 5
f = np.linspace(fstop,fstart,points)
f1 = f[::-1]    
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)

vnapwrstart = -28
vnapwrstop = 4
vnapwr_points = 17
vnapwr =np.linspace(vnapwrstart,vnapwrstop,vnapwr_points)


gs_global.setlevel_slow(Currentvalue,0.005)
time1 = time.time()
S_data = []
for i in tqdm(range(len(vnapwr))):
    vna.set_power(vnapwr[i])
    data = vna.get_data() 
    S_data.append(data)
    S_dataplot = np.rot90(np.array(S_data))
    S_dataplot_amp = 20*np.log10(np.abs(S_dataplot))
    if i!=0:
        plt.clf()
    imshow(S_dataplot_amp,figname='ReflectionVSPower_{}_signal'.format(qubitfre),extent=[vnapwr[0],vnapwr[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Power(dBm)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 

Sdata=np.array(S_data)
Sdata_signal=Sdata



about=' BOX?; setup: Reflection:setup: in:port31,out:20; \n\
    Attenuation: Vna={}dB,inline=52dB=42+10dB,outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,\n vnapower from {}dBm to {}dBm \ current value {}mA'.format(RTatt,round(fstart,3),round(fstop,3),points,ifband,vnapwrstart,vnapwrstop,Currentvalue)



filenameafter = 'reflectionPower_qubit_{}_Att_{}_signal'.format(qubitfre,RTatt)+str(time.strftime("%H%M%S"))
filename1 = filename+filenameafter+'.npz'
np.savez(filename1,fre = f,vnapwr = vnapwr, Sdata_signal = Sdata_signal, about = about,Magneticpoints = vnapwr_points,Frequencypoints = points,CenterFrequency = qubitfre,Frequencybandwidth = qubit_bandwith)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2 = filename+filenameafter
plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+'.jpg',format='jpg',dpi=300)
heatmap(vnapwr,f,S_dataplot_amp,corlor_scale = 'jet',xlabel = 'Power(dBm)',ylabel = 'Frequency[Ghz]',title = filenameafter,filename = filename,errobar = [],zmin = 0,zmax = 1,zauto = True)



gs_global.setlevel_slow(NoiseCuvalue,0.005)

time1 = time.time()
S_data=[]
for i in tqdm(range(len(vnapwr))):
    vna.set_power(vnapwr[i])
    data = vna.get_data() 
    S_data.append(data)
    S_dataplot = np.rot90(np.array(S_data))
    S_dataplot_amp = 20*np.log10(np.abs(S_dataplot))
    if i!=0:
        plt.clf()
    imshow(S_dataplot_amp,figname='ReflectionVSPower_{}_noise'.format(qubitfre),extent=[vnapwr[0],vnapwr[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Power(dBm)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
Sdata=np.array(S_data)
Sdata_noise=Sdata


about=' BOX?; setup: Reflection:setup: in:port31,out:20; \n\
    Attenuation: Vna={}dB,inline=52dB=42+10dB,outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,\n vnapower from {}dBm to {}dBm \ current value {}mA'.format(RTatt,round(fstart,3),round(fstop,3),points,ifband,vnapwrstart,vnapwrstop,NoiseCuvalue)

filenameafter = 'reflectionVSPower_qubit_{}_Att_{}_noise'.format(qubitfre,RTatt)+str(time.strftime("%H%M%S"))
filename1=filename+filenameafter+'.npz'
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
np.savez(filename1,fre=f,vnapwr=vnapwr, Sdata_noise=Sdata_noise,about=about,Magneticpoints=vnapwr_points,Frequencypoints=points,CenterFrequency=qubitfre,Frequencybandwidth=qubit_bandwith)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+'.jpg',format='jpg',dpi=300)
heatmap(vnapwr,f,S_dataplot_amp,corlor_scale='jet',xlabel='Power(dBm)',ylabel='Frequency[Ghz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)


fre1 = f

vnapwr1 = vnapwr
Amp1 = np.abs(np.array(Sdata_signal))
Phase1 = np.angle(np.array(Sdata_signal))

Magneticpoints = vnapwr_points;
Frequencypoints = points;
CenterFrequency = qubitfre;
Frequencybandwidth = qubit_bandwith;
Frequency0 = CenterFrequency-0.5*Frequencybandwidth;
Frequency1 = CenterFrequency+0.5*Frequencybandwidth;


fre2 = f
#volt2=data2['volt']
vnapwr2 = vnapwr
Amp2 = np.abs(np.array(Sdata_noise))
Phase2 = np.angle(np.array(Sdata_noise))

Z1 = Sdata_signal          #Convert the polor axis to real part and imaginary part
Z2 = Sdata_noise
Normalized_Z = Z1/Z2;

# from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as C
from scipy import misc
import scipy
import cmath
from scipy import ndimage,optimize
#import numpy as np  

def circle(w,y,r1,r2,omega):
    Z=1 - (r1/r2)*((1-1j*(w-CenterFrequency)/r2)/(1+((w-CenterFrequency)/r2)**2+omega/(r1*r2)))-re
    return Z


def f_2(xy0):
    r1=xy0[0]
    r2=xy0[1]
    omega=xy0[2]
    Z=circle(w,re,r1,r2,omega) 
    return abs(Z)

def yy(r1,r2,omega,w):    
    y2=1-(r1/r2)*((1-1j*(w-CenterFrequency)/r2)/(1+((w-CenterFrequency)/r2)**2+omega/(r1*r2)))
    return np.real(y2),np.imag(y2)



import matplotlib.pyplot as plt  
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator, MaxNLocator,LogLocator,FixedFormatter 
#import numpy as np  
w = np.linspace(Frequency0,Frequency1,Frequencypoints)
w1 = np.linspace(Frequency0,Frequency1,10*Frequencypoints)
#fig = plt.figure(figsize=(8,6),dpi=100)
fig = plt.figure(figsize=(8,8))
ax1 = fig.add_axes([0.15, 0.15, 0.8, 0.8])
for index in range(1,Magneticpoints):
    Real = np.real(Normalized_Z[index,:])
    Imaginary = np.imag(Normalized_Z[index,:])

    #最小二乘/
    re = Normalized_Z[index,:]
    xyc, _ =  optimize.leastsq(f_2,[0.003051,0.002467,1.291045335883221e-05])
    r1,r2,omega = xyc
    # Z=circle(w,y,r1 ,r2)
    for k in range(1000):
        xyc, _ =  optimize.leastsq(f_2,[r1,r2,omega])
        r1,r2,omega = xyc
    print(r1,r2,omega,r1/(2*r2))#圆心和半径
    real_fit, imag_fit = yy(r1,r2,omega,w1)
    plt.scatter(Real, Imaginary, linewidth=2)
    plt.plot(real_fit, imag_fit, linewidth=2)
    # plt.plot(real_fit*1.6-0.49,imag_fit*1.4,linewidth=2)

    plt.xlabel('Re(r)',fontsize=20)
    plt.ylabel('Im(r)',fontsize=20)

for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_size(fontsize=20)
    # ax1.xaxis.set_major_locator(FixedLocator([-1, -0.5, 0, 0.5, 1]))
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_size(fontsize=20)
    # ax1.yaxis.set_major_locator(FixedLocator([-0.8, -0.4, 0, 0.4, 0.8]))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
#ax1.set_xlabel(r'Current (mA)',fontsize=36)
ax1.set_xlabel(r'Re(r)',fontsize=28)
ax1.set_ylabel(r'Im(r)',fontsize=28)
plt.xlim(-1.1,1.1)
plt.ylim(-0.8,0.8)
# ax1.autoscale(tight=True)
about=' BOX?; setup: Smith circles \n\
    Attenuation: Vna={}dB,inline=52dB=42+10dB,outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,\n\
    vnapower from {}dBm to {}dBm \ current value {}mA'.format(RTatt,round(fstart,3),round(fstop,3),points,ifband,vnapwrstart,vnapwrstop,NoiseCuvalue)

filenameafter='reflectionVSPower_qubit_{}_Att_{}_smith'.format(qubitfre,RTatt)+str(time.strftime("%H%M%S"))
filename3=filename+filenameafter+'.npz'
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename3=filename+filenameafter
plt.savefig(filename3+'.pdf',format='pdf',dpi=300)
plt.savefig(filename3+'.jpg',format='jpg',dpi=300)

plt.show()
#%%scan_setting
gs3.setCURRmode()
gs3.getlevel()
gs3.Start_OutPut()

Currentvaluelist=np.linspace(-0.4,-0.6,11)
gs3.setlevel_slow(Currentvaluelist[0],0.005)
qubitfrelist=[]
for i in range(len(Currentvaluelist)):
# def single_qubit_efficient()
        
    #find resonant point of qubit
    Currentvalue_test=Currentvaluelist[i]
    gs3.setlevel_slow(Currentvalue_test,0.005)
    
    vna = keysight_vna('vna1', trace=21)
    which='S21'
    # whichtrace=int(which[1]+which[2])
    fstart_test=5.71-0.05
    fstop_test= 5.71+0.05
    points_test=501
    ifband_test=20
    vnapwr_test=-15
    f=np.linspace(fstart_test,fstop_test,points_test)
    
    vna.set_power(vnapwr_test)
    vna.set_startstopFre(fstart_test,fstop_test)
    vna.set_points_band(points_test,ifband_test)
    exec('test=vna.get_data_'+which+'()')
    frevalue=round(f[list(todB(test)).index(min(todB(test)))],3)
    
    qubitfrelist.append(frevalue)
print(qubitfrelist)
#%%
i=1
qubitfre=qubitfrelist[i-1]
qubit_bandwith=0.02
Currentvalue=Currentvaluelist[i-1]



for i in range(1):
    
    Currentvalue=Currentvalue+i*0.05
    print(Currentvalue)
    gs3.setlevel_slow(Currentvalue,0.005)
    
    vna=keysight_vna('vna1',trace=21)
    which='S21'
    # whichtrace=int(which[1]+which[2])
    RTatt=60
    fstart=qubitfre-qubit_bandwith
    fstop=qubitfre+qubit_bandwith
    points=401
    ifband=30
    vnapwr=-20
    f=np.linspace(fstart,fstop,points)
    
    vna.set_power(vnapwr)
    vna.set_startstopFre(fstart,fstop)
    vna.set_points_band(points,ifband)
    exec('test=vna.get_data_'+which+'()')
    plt.figure('test')
    plt.plot(f,todB(test))    
    plt.xlabel('Frequency(Ghz)')
    plt.ylabel('Amplitude(dB)')
    plt.title('VNA Result')
#%%
for i in range(len(Currentvaluelist)):
    # def single_qubit_efficient()
    filename = '/home/machine1/E/Data/YingHu/shaping2023/Effciency&SmithCircles/'
    if not os.path.exists(filename):
        os.makedirs(filename)
        

    
    qubitfre = qubitfrelist[i]
    qubit_bandwith = 0.06
    
    Currentvalue=Currentvaluelist[i]
    NoiseCuvalue=-0.3
    
    vna=keysight_vna('vna1',trace=21)
    which='S21'
    # whichtrace=int(which[1]+which[2])
    RTatt=60
    fstart=qubitfre - qubit_bandwith/2
    fstop=qubitfre + qubit_bandwith/2
    points=201
    ifband=5
    f=np.linspace(fstop,fstart,points)
    f1=f[::-1]    
    vna.set_startstopFre(fstart,fstop)
    vna.set_points_band(points,ifband)
    
    vnapwrstart=-29
    vnapwrstop=5
    vnapwr_points=10
    vnapwr=np.linspace(vnapwrstart,vnapwrstop,vnapwr_points)
    
    
    # exec('test=vna.get_data_'+which+'()')
    # plt.figure('test')                      
    # plt.plot(f1,todB(test))    
    # plt.xlabel('Frequency(Ghz)')
    # plt.ylabel('Amplitude(dB)')
    # plt.title('VNA Result')
    
    gs3.setlevel_slow(Currentvalue,0.005)
    time1=time.time()
    S_data=[]
    for i in tqdm(range(len(vnapwr))):
        vna.set_power(vnapwr[i])
        data = vna.get_data() 
        S_data.append(data)
        S_dataplot = np.rot90(np.array(S_data))
        S_dataplot_amp = 20*np.log10(np.abs(S_dataplot))
        if i!=0:
            plt.clf()
        imshow(S_dataplot_amp,figname='ReflectionVSPower_{}_signal'.format(qubitfre),extent=[vnapwr[0],vnapwr[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Power(dBm)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 

    Sdata=np.array(S_data)
    Sdata_signal=Sdata
    
    
    about=' BOX?; setup: Reflection:setup: in:port31,out:20; \n\
        Attenuation: Vna={}dB,inline=52dB=42+10dB,outline=76dB(RT)+36dB(4K) \n\
        fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,\n vnapower from {}dBm to {}dBm \ current value {}mA'.format(RTatt,round(fstart,3),round(fstop,3),points,ifband,vnapwrstart,vnapwrstop,Currentvalue)



    filenameafter = 'reflectionPower_qubit_{}_Att_{}_signal'.format(qubitfre,RTatt)+str(time.strftime("%H%M%S"))
    filename1 = filename+filenameafter+'.npz'
    np.savez(filename1,fre = f,vnapwr = vnapwr, Sdata_signal = Sdata_signal, about = about,Magneticpoints = vnapwr_points,Frequencypoints = points,CenterFrequency = qubitfre,Frequencybandwidth = qubit_bandwith)
    about_txt=open(filename+filenameafter+'.txt','a')
    about_txt.write(str(about))
    about_txt.close()

    filename2 = filename+filenameafter
    plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
    plt.savefig(filename2+'.jpg',format='jpg',dpi=300)
    heatmap(vnapwr,f,S_dataplot_amp,corlor_scale = 'jet',xlabel = 'Power(dBm)',ylabel = 'Frequency[Ghz]',title = filenameafter,filename = filename,errobar = [],zmin = 0,zmax = 1,zauto = True)
    
    
    
    gs3.setlevel_slow(NoiseCuvalue,0.005)

    time1 = time.time()
    S_data=[]
    for i in tqdm(range(len(vnapwr))):
        vna.set_power(vnapwr[i])
        data = vna.get_data() 
        S_data.append(data)
        S_dataplot = np.rot90(np.array(S_data))
        S_dataplot_amp = 20*np.log10(np.abs(S_dataplot))
        if i!=0:
            plt.clf()
        imshow(S_dataplot_amp,figname='ReflectionVSPower_{}_noise'.format(qubitfre),extent=[vnapwr[0],vnapwr[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Power(dBm)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    Sdata=np.array(S_data)
    Sdata_noise=Sdata


    about=' BOX?; setup: Reflection:setup: in:port31,out:20; \n\
        Attenuation: Vna={}dB,inline=52dB=42+10dB,outline=76dB(RT)+36dB(4K) \n\
        fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,\n vnapower from {}dBm to {}dBm \ current value {}mA'.format(RTatt,round(fstart,3),round(fstop,3),points,ifband,vnapwrstart,vnapwrstop,NoiseCuvalue)

    filenameafter = 'reflectionVSPower_qubit_{}_Att_{}_noise'.format(qubitfre,RTatt)+str(time.strftime("%H%M%S"))
    filename1=filename+filenameafter+'.npz'
    #filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
    np.savez(filename1,fre=f,vnapwr=vnapwr, Sdata_noise=Sdata_noise,about=about,Magneticpoints=vnapwr_points,Frequencypoints=points,CenterFrequency=qubitfre,Frequencybandwidth=qubit_bandwith)
    about_txt=open(filename+filenameafter+'.txt','a')
    about_txt.write(str(about))
    about_txt.close()

    filename2=filename+filenameafter
    plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
    plt.savefig(filename2+'.jpg',format='jpg',dpi=300)
    heatmap(vnapwr,f,S_dataplot_amp,corlor_scale='jet',xlabel='Power(dBm)',ylabel='Frequency[Ghz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
    
    
    
    fre1 = f

    vnapwr1 = vnapwr
    Amp1 = np.abs(np.array(Sdata_signal))
    Phase1 = np.angle(np.array(Sdata_signal))

    Magneticpoints = vnapwr_points;
    Frequencypoints = points;
    CenterFrequency = qubitfre;
    Frequencybandwidth = qubit_bandwith;
    Frequency0 = CenterFrequency-0.5*Frequencybandwidth;
    Frequency1 = CenterFrequency+0.5*Frequencybandwidth;


    fre2 = f
    #volt2=data2['volt']
    vnapwr2 = vnapwr
    Amp2 = np.abs(np.array(Sdata_noise))
    Phase2 = np.angle(np.array(Sdata_noise))

    Z1 = Amp1*np.exp(1j*Phase1)          #Convert the polor axis to real part and imaginary part
    Z2 = Amp2*np.exp(1j*Phase2)
    Normalized_Z = Z1/Z2;

    # from qutip import *
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.constants as C
    from scipy import misc
    import scipy
    import cmath
    from scipy import ndimage,optimize
    #import numpy as np  

    def circle(w,y,r1,r2,omega):
        Z=1 - (r1/r2)*((1-1j*(w-CenterFrequency)/r2)/(1+((w-CenterFrequency)/r2)**2+omega/(r1*r2)))-re
        return Z


    def f_2(xy0):
        r1=xy0[0]
        r2=xy0[1]
        omega=xy0[2]
        Z=circle(w,re,r1,r2,omega) 
        return abs(Z)

    def yy(r1,r2,omega,w):    
        y2=1-(r1/r2)*((1-1j*(w-CenterFrequency)/r2)/(1+((w-CenterFrequency)/r2)**2+omega/(r1*r2)))
        return np.real(y2),np.imag(y2)
    
    import matplotlib.pyplot as plt  
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator, MaxNLocator,LogLocator,FixedFormatter 
    #import numpy as np  
    w = np.linspace(Frequency0,Frequency1,Frequencypoints)
    w1 = np.linspace(Frequency0,Frequency1,10*Frequencypoints)
    #fig = plt.figure(figsize=(8,6),dpi=100)
    fig = plt.figure(figsize=(8,8))
    ax1 = fig.add_axes([0.15, 0.15, 0.8, 0.8])
    for index in range(1,Magneticpoints):
        Real = np.real(Normalized_Z[index,:])
        Imaginary = np.imag(Normalized_Z[index,:])

        #最小二乘/
        re = Normalized_Z[index,:]
        xyc, _ =  optimize.leastsq(f_2,[0.003051,0.002467,1.291045335883221e-05])
        r1,r2,omega = xyc
        # Z=circle(w,y,r1 ,r2)
        for k in range(1000):
            xyc, _ =  optimize.leastsq(f_2,[r1,r2,omega])
            r1,r2,omega = xyc
        print(r1,r2,omega,r1/(2*r2))#圆心和半径
        real_fit, imag_fit = yy(r1,r2,omega,w1)
        plt.scatter(Real, Imaginary, linewidth=2)
        plt.plot(real_fit, imag_fit, linewidth=2)
        # plt.plot(real_fit*1.6-0.49,imag_fit*1.4,linewidth=2)

        plt.xlabel('Re(r)',fontsize=20)
        plt.ylabel('Im(r)',fontsize=20)

    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_size(fontsize=20)
        # ax1.xaxis.set_major_locator(FixedLocator([-1, -0.5, 0, 0.5, 1]))
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_size(fontsize=20)
        # ax1.yaxis.set_major_locator(FixedLocator([-0.8, -0.4, 0, 0.4, 0.8]))
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
    #ax1.set_xlabel(r'Current (mA)',fontsize=36)
    ax1.set_xlabel(r'Re(r)',fontsize=28)
    ax1.set_ylabel(r'Im(r)',fontsize=28)
    plt.xlim(-1.1,1.1)
    plt.ylim(-0.9,0.9)
    # ax1.autoscale(tight=True)
    about=' BOX?; setup: Smith circles \n\
        Attenuation: Vna={}dB,inline=52dB=42+10dB,outline=76dB(RT)+36dB(4K) \n\
        fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,\n\
        vnapower from {}dBm to {}dBm \ current value {}mA'.format(RTatt,round(fstart,3),round(fstop,3),points,ifband,vnapwrstart,vnapwrstop,NoiseCuvalue)

    filenameafter='reflectionVSPower_qubit_{}_Att_{}_smith'.format(qubitfre,RTatt)+str(time.strftime("%H%M%S"))
    filename3=filename+filenameafter+'.npz'
    about_txt=open(filename+filenameafter+'.txt','a')
    about_txt.write(str(about))
    about_txt.close()

    filename3=filename+filenameafter
    plt.savefig(filename3+'.pdf',format='pdf',dpi=300)
    plt.savefig(filename3+'.jpg',format='jpg',dpi=300)

    plt.show()
#%%repeat sweep qubit energy spectrum 
qubitnum=1

gs6.setlevel_slow(-0.75,0.01)#QUBIT4
gs3.setlevel_slow(-0.25,0.01)#QUBIT6
gs1.setlevel_slow(-0.0,0.01)#QUBIT3
gs7.setlevel_slow(-1.25,0.01)#QUBIT1
gs9.setlevel_slow(-1,0.01)#QUBIT5

vna=keysight_vna('vna1',trace=21)
which='S21'
# whichtrace=int(which[1]+which[2])
RTatt=50
fstart=4
fstop=5.8
points=501
ifband=15
vnapwr=-30
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

Globalnum=3
volt=np.linspace(-1.25,0.5,101)
gs7.Start_OutPut()
gs7.setlevel_slow(volt[0],0.01)
time1=time.time()
S_data=[]
for i in range(len(volt)):
    gs7.setlevel_slow(volt[i],0.01)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    if i!=0:
        plt.clf()
    imshow(S_dataplot_amp,figname='Energyspetrum{}_qubit{}'.format(Globalnum,qubitnum),extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    
    if i%(int(len(volt)/10))==0:
        usedtime=time.time()-time1
        time1=time.time()
        remaintime=(10-((i+1)*10/len(volt)))*usedtime
        print('scan_Fr: {:.0f}% completed.'.format(i/len(volt)*100))
        print('remain about {:.0f} mins.'.format(round(remaintime/60,2)))
Sdata=np.array(S_data)

about=' BOX?; setup: in:port22,out:17; \n Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB,             outline=76dB(RT)+36dB(4K) \n fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])+'\n qubit1 current {}mA,\n qubit2 current xmA,\n qubit3 current {}mA,\n qubit4 current {}mA,\n qubit5 current {}mA,\n qubit6 current {}mA,\n'.format(round(float(gs7.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs1.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs6.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs9.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs3.Inst.query(':SOUR:LEV?'))*10e2,3))

filenameafter='Energyspetrum{}_qubit{}_vnapwr{}_RTatt{}_Cinch3_1_GS_07_'.format(Globalnum,qubitnum,vnapwr,RTatt)+str(time.strftime("%H%M%S"))
filename1=filename+filenameafter+'.npz'
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+'.jpg',format='jpg',dpi=300)
heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)











qubitnum=3

gs6.setlevel_slow(-0.75,0.01)#QUBIT4
gs3.setlevel_slow(-0.25,0.01)#QUBIT6
gs1.setlevel_slow(-2.2,0.01)#QUBIT3
gs7.setlevel_slow(0.35,0.01)#QUBIT1
gs9.setlevel_slow(-1,0.01)#QUBIT5

vna=keysight_vna('vna1',trace=21)
which='S21'
# whichtrace=int(which[1]+which[2])
RTatt=50
fstart=4
fstop=5.8
points=501
ifband=15
vnapwr=-30
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

Globalnum=3
volt=np.linspace(-2.2,0.2,101)
gs1.Start_OutPut()
gs1.setlevel_slow(volt[0],0.01)
time1=time.time()
S_data=[]
for i in range(len(volt)):
    gs1.setlevel_slow(volt[i],0.01)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    if i!=0:
        plt.clf()
    imshow(S_dataplot_amp,figname='Energyspetrum{}_qubit{}'.format(Globalnum,qubitnum),extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    
    if i%(int(len(volt)/10))==0:
        usedtime=time.time()-time1
        time1=time.time()
        remaintime=(10-((i+1)*10/len(volt)))*usedtime
        print('scan_Fr: {:.0f}% completed.'.format(i/len(volt)*100))
        print('remain about {:.0f} mins.'.format(round(remaintime/60,2)))
Sdata=np.array(S_data)

about=' BOX?; setup: in:port22,out:17; \n Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB,             outline=76dB(RT)+36dB(4K) \n fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])+'\n qubit1 current {}mA,\n qubit2 current xmA,\n qubit3 current {}mA,\n qubit4 current {}mA,\n qubit5 current {}mA,\n qubit6 current {}mA,\n'.format(round(float(gs7.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs1.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs6.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs9.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs3.Inst.query(':SOUR:LEV?'))*10e2,3))

filenameafter='Energyspetrum{}_qubit{}_vnapwr{}_RTatt{}_Cinch3_3_GS_01_'.format(Globalnum,qubitnum,vnapwr,RTatt)+str(time.strftime("%H%M%S"))
filename1=filename+filenameafter+'.npz'
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+'.jpg',format='jpg',dpi=300)
heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)












qubitnum=4

gs6.setlevel_slow(-1,0.01)#QUBIT4
gs3.setlevel_slow(-0.25,0.01)#QUBIT6
gs1.setlevel_slow(-0.0,0.01)#QUBIT3
gs7.setlevel_slow(0.35,0.01)#QUBIT1
gs9.setlevel_slow(-1,0.01)#QUBIT5

vna=keysight_vna('vna1',trace=21)
which='S21'
# whichtrace=int(which[1]+which[2])
RTatt=50
fstart=4
fstop=5.8
points=501
ifband=15
vnapwr=-30
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

Globalnum=3
volt=np.linspace(-1,1,101)
gs6.Start_OutPut()
gs6.setlevel_slow(volt[0],0.01)
time1=time.time()
S_data=[]
for i in range(len(volt)):
    gs6.setlevel_slow(volt[i],0.01) 
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    if i!=0:
        plt.clf()
    imshow(S_dataplot_amp,figname='Energyspetrum{}_qubit{}'.format(Globalnum,qubitnum),extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    
    if i%(int(len(volt)/10))==0:
        usedtime=time.time()-time1
        time1=time.time()
        remaintime=(10-((i+1)*10/len(volt)))*usedtime
        print('scan_Fr: {:.0f}% completed.'.format(i/len(volt)*100))
        print('remain about {:.0f} mins.'.format(round(remaintime/60,2)))
Sdata=np.array(S_data)

about=' BOX?; setup: in:port22,out:17; \n Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB,             outline=76dB(RT)+36dB(4K) \n fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])+'\n qubit1 current {}mA,\n qubit2 current xmA,\n qubit3 current {}mA,\n qubit4 current {}mA,\n qubit5 current {}mA,\n qubit6 current {}mA,\n'.format(round(float(gs7.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs1.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs6.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs9.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs3.Inst.query(':SOUR:LEV?'))*10e2,3))

filenameafter='Energyspetrum{}_qubit{}_vnapwr{}_RTatt{}_Cinch3_4_GS_06_'.format(Globalnum,qubitnum,vnapwr,RTatt)+str(time.strftime("%H%M%S"))
filename1=filename+filenameafter+'.npz'
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+'.jpg',format='jpg',dpi=300)
heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)







qubitnum=5

gs6.setlevel_slow(-0.75,0.01)#QUBIT4
gs3.setlevel_slow(-0.25,0.01)#QUBIT6
gs1.setlevel_slow(-0.0,0.01)#QUBIT3
gs7.setlevel_slow(0.35,0.01)#QUBIT1
gs9.setlevel_slow(-1,0.01)#QUBIT5

vna=keysight_vna('vna1',trace=21)
which='S21'
# whichtrace=int(which[1]+which[2])
RTatt=50
fstart=4
fstop=5.8
points=501
ifband=15
vnapwr=-30
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

Globalnum=3
volt=np.linspace(-1,1,101)
gs9.Start_OutPut()
gs9.setlevel_slow(volt[0],0.01)
time1=time.time()
S_data=[]
for i in range(len(volt)):
    gs9.setlevel_slow(volt[i],0.01)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    if i!=0:
        plt.clf()
    imshow(S_dataplot_amp,figname='Energyspetrum{}_qubit{}'.format(Globalnum,qubitnum),extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    
    if i%(int(len(volt)/10))==0:
        usedtime=time.time()-time1
        time1=time.time()
        remaintime=(10-((i+1)*10/len(volt)))*usedtime
        print('scan_Fr: {:.0f}% completed.'.format(i/len(volt)*100))
        print('remain about {:.0f} mins.'.format(round(remaintime/60,2)))
Sdata=np.array(S_data)

about=' BOX?; setup: in:port22,out:17; \n Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB,             outline=76dB(RT)+36dB(4K) \n fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])+'\n qubit1 current {}mA,\n qubit2 current xmA,\n qubit3 current {}mA,\n qubit4 current {}mA,\n qubit5 current {}mA,\n qubit6 current {}mA,\n'.format(round(float(gs7.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs1.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs6.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs9.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs3.Inst.query(':SOUR:LEV?'))*10e2,3))

filenameafter='Energyspetrum{}_qubit{}_vnapwr{}_RTatt{}_Cinch3_5_GS_09_'.format(Globalnum,qubitnum,vnapwr,RTatt)+str(time.strftime("%H%M%S"))
filename1=filename+filenameafter+'.npz'
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+'.jpg',format='jpg',dpi=300)
heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)









qubitnum=6

gs6.setlevel_slow(-0.75,0.01)#QUBIT4
gs3.setlevel_slow(-1.7,0.01)#QUBIT6
gs1.setlevel_slow(-0.0,0.01)#QUBIT3
gs7.setlevel_slow(0.35,0.01)#QUBIT1
gs9.setlevel_slow(-1,0.01)#QUBIT5

vna=keysight_vna('vna1',trace=21)
which='S21'
# whichtrace=int(which[1]+which[2])
RTatt=50
fstart=4
fstop=5.8
points=501
ifband=15
vnapwr=-30
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

Globalnum=3
volt=np.linspace(-1.7,-0.3,101)
gs3.Start_OutPut()
gs3.setlevel_slow(volt[0],0.01)
time1=time.time()
S_data=[]
for i in range(len(volt)):
    gs3.setlevel_slow(volt[i],0.01)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    if i!=0:
        plt.clf()
    imshow(S_dataplot_amp,figname='Energyspetrum{}_qubit{}'.format(Globalnum,qubitnum),extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    
    if i%(int(len(volt)/10))==0:
        usedtime=time.time()-time1
        time1=time.time()
        remaintime=(10-((i+1)*10/len(volt)))*usedtime
        print('scan_Fr: {:.0f}% completed.'.format(i/len(volt)*100))
        print('remain about {:.0f} mins.'.format(round(remaintime/60,2)))
gs3.setlevel_slow(-0.25,0.01) 
Sdata=np.array(S_data)

about=' BOX?; setup: in:port22,out:17; \n Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB,             outline=76dB(RT)+36dB(4K) \n fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])+'\n qubit1 current {}mA,\n qubit2 current xmA,\n qubit3 current {}mA,\n qubit4 current {}mA,\n qubit5 current {}mA,\n qubit6 current {}mA,\n'.format(round(float(gs7.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs1.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs6.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs9.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs3.Inst.query(':SOUR:LEV?'))*10e2,3))

filenameafter='Energyspetrum{}_qubit{}_vnapwr{}_RTatt{}_Cinch3_6_GS_03_'.format(Globalnum,qubitnum,vnapwr,RTatt)+str(time.strftime("%H%M%S"))
filename1=filename+filenameafter+'.npz'
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+'.jpg',format='jpg',dpi=300)
heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)





#%%linewidth with sweet point
qubitnum=4
vna=keysight_vna('vna1',trace=21)
which='S21'
# whichtrace=int(which[1]+which[2])
RTatt=50
fstart=5
fstop=5.25
points=501
ifband=5
vnapwr=-30
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)

data_amp=vna.get_data()
S_dataplot_amp=20*np.log10(np.abs(data_amp))


gs6.setlevel_slow(-1,0.02)#QUBIT4
gs3.setlevel_slow(-0.25,0.01)#QUBIT6
gs1.setlevel_slow(-0.0,0.01)#QUBIT3
gs7.setlevel_slow(0.35,0.01)#QUBIT1
gs9.setlevel_slow(-1,0.01)#QUBIT5


data_bg=vna.get_data()
S_dataplot_bg=20*np.log10(np.abs(data_bg))
S_dataplot_deal=20*np.log10(np.abs(data_amp)/np.abs(data_bg))

filenameafter='Energyspetrumlinewidth_qubit{}_vnapwr{}_RTatt{}_Cinch3_3_GS_01_'.format(qubitnum,vnapwr,RTatt)+str(time.strftime("%H%M%S"))
filename1=filename+filenameafter
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
np.savez(filename1,fre=f1,data_amp=data_amp,S_dataplot_amp=S_dataplot_amp,data_bg=data_bg,S_dataplot_deal=S_dataplot_deal)

plt.plot(f1,S_dataplot_amp)
plt.plot(f1,S_dataplot_bg)
plt.plot(f1,S_dataplot_deal)
plt.savefig(filename1+'.jpg',format='jpg',dpi=300)
plt.show()


# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   



#%%first-second level transition
which='S21'
# whichtrace=int(which[1]+which[2])
RTatt=40
fstart=4.6
fstop=5.4
points=501
ifband=30
vnapwr=-10
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')

Globalnum=9
volt=np.linspace(-1,0,101)
gs3.Start_OutPut()
gs3.setlevel_slow(volt[0],0.01)
time1=time.time()
S_data=[]
for i in range(len(volt)):
    gs3.setlevel_slow(volt[i],0.01)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    if i!=0:
        plt.clf()
    imshow(S_dataplot_amp,figname='Globalsweep{}'.format(Globalnum),extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    if (i+1)%(int(len(volt)/10))==0:
        usedtime=time.time()-time1
        time1=time.time()
        remaintime=(10-((i+1)*10/len(volt)))*usedtime
        print('scan_Fr: {:.0f}% completed.'.format(i/len(volt)*100))
        print('remain {:.0f} seconds.'.format(remaintime))
gs3.setlevel_slow(0,0.02) 
Sdata=np.array(S_data)
about=' BOX?; setup: in:port18,out:17; \n Attenuation: Vna=50dB=40+10dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n current from{}mA to{}mA'.format(fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])
filename1=filename+'Globalsweep{}_Att_{}_vnapwr_{}_Cinch3_1_GS_03_1to2transition'.format(Globalnum,RTatt,vnapwr)+str(time.strftime("%H%M%S"))+'.npz'
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
filename2=filename+'Globalsweep{}_Att_{}_vnapwr_{}_Cinch3_1_GS_03_1to2transition'.format(Globalnum,RTatt,vnapwr)+str(time.strftime("%H%M%S"))
plt.savefig(filename2+'.pdf',format='pdf',dpi=600)
plt.savefig(filename2+'.jpg',format='jpg',dpi=600)
heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Globalsweep{}_Att_{}_vnapwr_{}_Cinch3_1_GS_03_1to2transition'.format(Globalnum,RTatt,vnapwr)+str(time.strftime("%H%M%S")),filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)


#%% The nearest qubit-qubit coupling
#gs3.setlevel_slow(-0.84,0.01)#4.8GHZ
#gs4=GS_200('DC4')
#gs4.setlevel_slow(0.92,0.01)#4.62GHZ
#gs4=GS_200('DC4')
gs3.setlevel_slow(0,0.005)#-0.8--4.634GHZ

which='S21'
# whichtrace=int(which[1]+which[2])
RTatt=40
fstart=5.03 
fstop=5.5
points=401
ifband=30
vnapwr=-30
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')
Globalnum=7
volt=np.linspace(-1.12-0.05,-1.12+0.05,101)
gs6.Start_OutPut()
gs6.setlevel_slow(volt[0],0.005)
time1=time.time()
S_data=[]
for i in range(len(volt)):
    gs6.setlevel_slow(volt[i],0.01)
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    if i!=0:
        plt.clf()
    imshow(S_dataplot_amp,figname='Globalsweep{}'.format(Globalnum),extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    if i%(int(len(volt)/10))==0:
        usedtime=time.time()-time1
        time1=time.time()
        remaintime=(10-((i+1)*10/len(volt)))*usedtime
        print('scan_Fr: {:.0f}% completed.'.format((i+1)/len(volt)*100))
        print('{},remain {:.0f} seconds.'.format(time1,remaintime))
#gs3.setlevel_slow(0,0.02) 
Sdata=np.array(S_data)
about=' BOX?; setup: in:port18,out:17; \n Attenuation: Vna=70dB=40+30dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n current from{}mA to{}mA'.format(fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])
filename1=filename+'Globalsweep{}_Att_{}_Cinch3_3_GS_02(steady)_Cinch3_4_GS_06_'.format(Globalnum,RTatt)+str(time.strftime("%H%M%S"))+'.npz'
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)

filename2=filename+'Globalsweep{}_Att_{}_Cinch3_3_GS_02(steady)_Cinch3_4_GS_06_'.format(Globalnum,RTatt)+str(time.strftime("%H%M%S"))
plt.savefig(filename2+'.pdf',format='pdf',dpi=600)
plt.savefig(filename2+'.jpg',format='jpg',dpi=600)
heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Globalsweep{}_Att_{}_Cinch3_3_GS_02(steady)_Cinch3_4_GS_06_'.format(Globalnum,RTatt)+str(time.strftime("%H%M%S")),filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)

#%%  change vnapower measure cavity fre
cavityfre=6.729

which='S21'
# whichtrace=int(which[1]+which[2])
RTatt=16
fstart=cavityfre - 0.002
fstop=cavityfre + 0.002
points=201
ifband=30
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)

vnapwrstart=-10
vnapwrstop=5
vnapwr_points=16
vnapwr=np.linspace(vnapwrstart,vnapwrstop,vnapwr_points)

time1=time.time()
S_data=[]
for i in range(len(vnapwr)):
    vna.set_power(vnapwr[i])
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    if i!=0:
        plt.clf()
    imshow(S_dataplot_amp,figname='powersweep{}_'.format(cavityfre),extent=[vnapwr[0],vnapwr[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Power(dBm)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    if i%(int(len(vnapwr)/10))==0:
        print('scan_Fr: {:.0f}% completed.'.format((i+1)/len(vnapwr)*100))
time2=time.time()    
time_cost=time2-time1
print('time_cost',time_cost)   
Sdata=np.array(S_data)
about=' BOX?; setup: in:port26,out:29; \n Attenuation: Vna={}dB(RT), inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz, \n vnapower from {}dBm to {}dBm'.format(RTatt,fstart,fstop,points,ifband,vnapwrstart,vnapwrstop)
filename1=filename+'Powersweep{}_Att_{}_'.format(cavityfre,RTatt)+str(time.strftime("%H%M%S"))+'.npz'
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
np.savez(filename1,fre=f,vnapwr=vnapwr,Sdata=Sdata,about=about)

filename2=filename+'Powersweep{}_Att_{}_'.format(cavityfre,RTatt)+str(time.strftime("%H%M%S"))
plt.savefig(filename2+'.pdf',format='pdf',dpi=600)
plt.savefig(filename2+'.jpg',format='jpg',dpi=600)
heatmap(vnapwr,f,S_dataplot_amp,corlor_scale='jet',xlabel='Power(dBm)',ylabel='Frequency[Ghz]',title='Powersweep{}_Att_{}_'.format(cavityfre,RTatt)+str(time.strftime("%H%M%S")),filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
#%%sweep setting
qubitnum=4
# gs1.setlevel_slow(0,0.01)#QUBIT2
# gs7.setlevel_slow(-1.648,0.01)#QUBIT3
# gs6.setlevel_slow(-0.75,0.01)#QUBIT4
# gs9.setlevel_slow(-1,0.01)#QUBIT5
# gs3.setlevel_slow(-0.25,0.01)#QUBIT6

qubitfre=5.21
qubit_bandwith=0.5

Currentvalue=0.14
for i in range(3):

    Currentvalue=Currentvalue+i*0.05
    gs6.setlevel_slow(Currentvalue,0.01)
    
    vna=keysight_vna('vna1',trace=21)
    which='S21'
    # whichtrace=int(which[1]+which[2])
    RTatt=60
    fstart=qubitfre - qubit_bandwith/2
    fstop=qubitfre + qubit_bandwith/2
    points=501
    ifband=30
    vnapwr=-30
    f=np.linspace(fstart,fstop,points)
    
    vna.set_power(vnapwr)
    vna.set_startstopFre(fstart,fstop)
    vna.set_points_band(points,ifband)
    exec('test=vna.get_data_'+which+'()')
    plt.figure('test')
    plt.plot(f,todB(test))    
    plt.xlabel('Frequency(Ghz)')
    plt.ylabel('Amplitude(dB)')
    plt.title('VNA Result')
#%%  efficiency with Smith figture
qubitnum=3
# gs1.setlevel_slow(0,0.01)#QUBIT2
# gs7.setlevel_slow(-1.648,0.01)#QUBIT3
# gs6.setlevel_slow(-0.75,0.01)#QUBIT4
# gs9.setlevel_slow(-1,0.01)#QUBIT5
# gs3.setlevel_slow(-0.25,0.01)#QUBIT6

qubitfre=5.028
qubit_bandwith=0.2

Currentvalue=-1.24
NoiseCuvalue=0

vna=keysight_vna('vna1',trace=21)
which='S21'
# whichtrace=int(which[1]+which[2])
RTatt=60
fstart=qubitfre - qubit_bandwith/2
fstop=qubitfre + qubit_bandwith/2
points=301
ifband=5
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)

vnapwrstart=-30
vnapwrstop=6
vnapwr_points=19
vnapwr=np.linspace(vnapwrstart,vnapwrstop,vnapwr_points)

print("qubit signal measure running\n")
print("Controling bias current at qubit{}'s working point, {}GHz".format(qubitnum,qubitfre))
gs7.setlevel_slow(Currentvalue,0.01)
print("Finished!\n")
time1=time.time()
S_data=[]
for i in range(len(vnapwr)):
    vna.set_power(vnapwr[i])
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    if i!=0:
        plt.clf()
    imshow(S_dataplot_amp,figname='ReflectionVSPower_{}_signal'.format(qubitfre),extent=[vnapwr[0],vnapwr[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Power(dBm)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    if (i)%(int(len(vnapwr)/10))==0:
        usedtime=time.time()-time1
        time1=time.time()
        remaintime=(10-((i+1)*10/len(vnapwr)))*usedtime
        print('scan_Fr: {:.0f}% completed.'.format(i/len(vnapwr)*100))
        print('remain about {:.0f} mins.'.format(round(remaintime/60,1)))  
Sdata=np.array(S_data)
Sdata_signal=Sdata

print("Saving the working-point data\n Waiting the current to bias away from the working-point ")
about=' BOX?; setup: in:port22,out:17; \n Attenuation: Vna={}dB(RT), inline=52dB=42+10dB,             outline=76dB(RT)+36dB(4K) \n fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,\n vnapower from {}dBm to {}dBm \ current value {}mA'.format(RTatt,round(fstart,3),round(fstop,3),points,ifband,vnapwrstart,vnapwrstop,Currentvalue)+'\n qubit1 current 0.4mA(ezq),\n qubit2 current {}mA,\n qubit3 current {}mA,\n qubit4 current {}mA,\n qubit5 current {}mA,\n qubit6 current {}mA,\n'.format(round(float(gs1.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs7.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs6.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs9.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs3.Inst.query(':SOUR:LEV?'))*10e2,3))
print(str(about))
filenameafter='ReflectionVSPower_qubit{}_{}_Att_{}_signal'.format(qubitnum,qubitfre,RTatt)+str(time.strftime("%H%M%S"))
filename1=filename+filenameafter+'.npz'
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
np.savez(filename1,fre=f,vnapwr=vnapwr,Sdata=Sdata,about=about,Magneticpoints=vnapwr_points,Frequencypoints=points,CenterFrequency=qubitfre,Frequencybandwidth=qubit_bandwith)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+'.jpg',format='jpg',dpi=300)
heatmap(vnapwr,f,S_dataplot_amp,corlor_scale='jet',xlabel='Power(dBm)',ylabel='Frequency[Ghz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)




print("\n Noise background measure running\n")
print("Controling bias current away from qubit{}'s working point, {}GHz".format(qubitnum,qubitfre))
gs1.setlevel_slow(NoiseCuvalue,0.01)
print("Finished!\n")


time1=time.time()
S_data=[]
for i in range(len(vnapwr)):
    vna.set_power(vnapwr[i])
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    if i!=0:
        plt.clf()
    imshow(S_dataplot_amp,figname='ReflectionVSPower_{}_noise'.format(qubitfre),extent=[vnapwr[0],vnapwr[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Power(dBm)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    if (i)%(int(len(vnapwr)/10))==0:
        usedtime=time.time()-time1
        time1=time.time()
        remaintime=(10-((i+1)*10/len(vnapwr)))*usedtime
        print('scan_Fr: {:.0f}% completed.'.format(i/len(vnapwr)*100))
        print('remain about {:.0f} mins.'.format(round(remaintime/60,1)))  
Sdata=np.array(S_data)
Sdata_noise=Sdata

print("Saving the noise-background data\n Waiting for the Smith figture ... ")
about=' BOX?; setup: in:port22,out:17; \n Attenuation: Vna={}dB(RT), inline=60dB=42+8+10dB,             outline=76dB(RT)+36dB(4K) \n fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,\n vnapower from {}dBm to {}dBm \ noise_current value {}mA'.format(RTatt,round(fstart,3),round(fstop,3),points,ifband,vnapwrstart,vnapwrstop,Currentvalue)+'\n qubit1 current 0.4mA(ezq),\n qubit2 current {}mA,\n qubit3 current {}mA,\n qubit4 current {}mA,\n qubit5 current {}mA,\n qubit6 current {}mA,\n'.format(round(float(gs1.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs7.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs6.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs9.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs3.Inst.query(':SOUR:LEV?'))*10e2,3))
print(str(about))

filenameafter='ReflectionVSPower_qubit{}_{}_Att_{}_noise'.format(qubitnum,qubitfre,RTatt)+str(time.strftime("%H%M%S"))
filename1=filename+filenameafter+'.npz'
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
np.savez(filename1,fre=f,vnapwr=vnapwr,Sdata=Sdata,about=about,Magneticpoints=vnapwr_points,Frequencypoints=points,CenterFrequency=qubitfre,Frequencybandwidth=qubit_bandwith)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+'.jpg',format='jpg',dpi=300)
heatmap(vnapwr,f,S_dataplot_amp,corlor_scale='jet',xlabel='Power(dBm)',ylabel='Frequency[Ghz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)



fre1=f

vnapwr1=vnapwr
Amp1=np.abs(np.array(Sdata_signal))
Phase1=np.angle(np.array(Sdata_signal))

Magneticpoints=vnapwr_points;
Frequencypoints=points;
CenterFrequency=qubitfre;
Frequencybandwidth=qubit_bandwith;
Frequency0=CenterFrequency-0.5*Frequencybandwidth;
Frequency1=CenterFrequency+0.5*Frequencybandwidth;


fre2=f
#volt2=data2['volt']
vnapwr2=vnapwr
Amp2=np.abs(np.array(Sdata_noise))
Phase2=np.angle(np.array(Sdata_noise))

Z1=Amp1*np.exp(1j*Phase1)          #Convert the polor axis to real part and imaginary part
Z2=Amp2*np.exp(1j*Phase2)
Normalized_Z=Z1/Z2;

from matplotlib.ticker import ScalarFormatter, FormatStrFormatter  
#import numpy as np  

#fig = plt.figure(figsize=(8,6),dpi=100)
fig = plt.figure()
ax1 = fig.add_axes([0.18, 0.16, 0.76, 0.8])
for index in range(0,19):
    Real=np.real(Normalized_Z[index,:])
    Imaginary=np.imag(Normalized_Z[index,:])
    Power=-134+index*2
#    Label=unicode('Power'+'dBm').decode('utf8','ignore')
    plt.plot(Real,Imaginary,linewidth=2)
    
    
#labs = [l.get_label() for l in lns]
#ax1.legend(lns, labs, 'upper right',fontsize=20)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_size(fontsize=24)
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_size(fontsize=24)
#label_f1 = "-136dBm"
#plt.text(0.7, 0.8,label_f1,
#     horizontalalignment='left',
#     verticalalignment='center',fontsize=20,
#     transform = ax1.transAxes)
#plt.annotate('-116dBm',xy=(-0.5,0),xytext=(0.2,0.2),textcoords='offset points',arrowprops=dict(facecolor='green'),fontsize=20)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%0.3f'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
#ax1.set_xlabel(r'Current (mA)',fontsize=36)
ax1.set_xlabel(r'Re(r)',fontsize=28)
ax1.set_ylabel(r'Im(r)',fontsize=28)
# plt.xlim(-0.6,1.2)
# plt.ylim(-0.9,0.9)


ax1.autoscale(tight=True)
plt.savefig(filename2+'.png',format='png',dpi=300)


plt.show()
#%%phase
phase=np.angle(S_data)
fig=plt.figure()
ax2=fig.add_axes([0.2,0.2,0.75,0.75])
c=plt.imshow(phase.T,aspect='auto',interpolation='gaussian',origin='lower',
             extent=[volt.min(),volt.max(),f.min(),f.max()],
             vmax=phase.max(),vmin=phase.min())
c.set_cmap('bwr')
cbar=fig.colorbar(c,ticks=[-4,0,4])
cbar.set_label("Gain(dB)",rotation='270',labelpad=20,fontsize=20)

#%%sweep cavity energy specturm
qubitnum = 2
vna = keysight_vna('vna2',trace=21)
which = 'S21'
RTatt = 70

fstart, fstop, points = 6.56, 6.57, 101 

ifband = 30
vnapwr = np.linspace(-30,6,19)

f = np.linspace(fstop,fstart,points)
f1 = f[::-1]  
  

vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# exec('test=vna.get_data_'+which+'()')
# plt.figure('test')
# plt.plot(f[::-1],todB(test))    
# plt.xlabel('Frequency(Ghz)')
# plt.ylabel('Amplitude(dB)')   

filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_230609C/cavity_spectrum/{time.strftime("%m%d")}_{time.strftime("%H%M%S")}_freq_{fstart}GHz_to_{fstop}GHz/'
if not os.path.exists(filename):
    os.makedirs(filename)

Globalnum = 4
filenameafter = 'Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs_qubit_local_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))
filename1=filename+filenameafter


S_data=[]
for i in tqdm(range(len(vnapwr))):
    vna.set_power(vnapwr[i])
    data=vna.get_data() 
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    
    Phase = np.angle(np.array(S_data))
    S_dataplot_phase = np.rot90(phase_handle(Phase, f1, 52.65))

    
    heatmap(np.linspace(vnapwr[0], vnapwr[i], i), f, S_dataplot_amp, corlor_scale='bwr', xlabel='Current[mA]',
            ylabel='Frequency[Ghz]', title='sweep_Amp', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)
    heatmap(np.linspace(vnapwr[0], vnapwr[i], i), f, S_dataplot_phase, corlor_scale='bwr', xlabel='Current[mA]',
            ylabel='Frequency[Ghz]', title='sweep_Phase', filename=filename1, errobar=[], zmin=0, zmax=1, zauto=True)

# gs3.setlevel_slow(-0,0.01) 
Sdata=np.array(S_data)

filenameafter='Energyspetrum{}_qubit{}_with_qubit_f_{}_to_{}_vnapwr{}_RTatt{}_Cinch3_{}_gs1_'.format(Globalnum,qubitnum,fstart,fstop,vnapwr,RTatt,qubitnum)+str(time.strftime("%m%d%H%M%S"))

about=' BOX?; setup: in:port 7,out: 13; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    current from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,vnapwr[0],vnapwr[-1])

filename1=filename+filenameafter+'.npz'
np.savez(filename1,fre=f,vnapwr=vnapwr,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename+filenameafter+'.txt','a')
about_txt.write(str(about))
about_txt.close()

filename2=filename+filenameafter
imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[vnapwr[0],vnapwr[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+'.jpg',format='jpg',dpi=300)

imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[vnapwr[0],vnapwr[i],f[-1],f[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+'.jpg',format='jpg',dpi=300)

heatmap(vnapwr,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(vnapwr,f,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
#%%初始化网分 twotone
vna=keysight_vna('vna2',trace=43)
vna.two_tone_vna()
fstart=6.449#12.97#6.471#6.483#6.479#6.47
fstop=6.449#12.97#6.471#6.483#6.479#6.47
points=201
RTatt=60
ifband=20#每秒扫描次数
vnapwr=-10#-12
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
#%%初始化微波源 twotone
mw=MW('mw1')
mw.MW_setpower(-130)
mwpwr=-15

mwfre_start=4
mwfre_stop=10
points=points
mw.VnaTrig(mwfre_start,mwfre_stop,points)
mw.MW_setpower(mwpwr)
mw.start_Output()
fwr=np.linspace(mwfre_start,mwfre_stop,points)
fwr1=fwr[::-1]
# mw.Close_Output()
#%% twotone sweep z
twotonenum = 1#0.015#-0.016 #-0.0308#0.382;
# gs5.setlevel_slowV2(-0.2,0.0001)
# gs5.setlevel_slowV2(Global,0.0001)
# gs3.setlevel_slowV2(1.5,0.005)


# bg = vna.get_data()
# gs1.setlevel_slowV2(-0.24,0.005)
volt=np.linspace( 0.5, 2, 51)[::1]+0.4
alpha_volt = np.linspace(1, 5, 5    )
# gs6.Start_OutPut() #右z
# time.sleep(339)
gs3.setlevel_slowV2(volt[0],0.005)#-0.9,0.05
S_data=[]
i=0
mw.MW_setpower(mwpwr)
print(time.strftime("%H:%M:%S"))
for alpha in alpha_volt:
    filename = f'/home/machine1/E/Data/LSY_new/data/CELv1/twotone_energyspectrum_RightQubit_RZ(alpha={alpha}mA)/two_tone_qubit_freq_{mwfre_start}GHz_to_{mwfre_stop}GHz_cavity_freq_{fstart}GHz{time.strftime("%m%d")}/'
    if not os.path.exists(filename):
        os.makedirs(filename)
    name0=f"two_tone_qubit_freq_{mwfre_start}GHz_to_{mwfre_stop}GHz_cavity_freq_{fstart}GHz"
    filename1=filename+name0+time.strftime("%H%M%S")+'_'
    gs7.setlevel_slowV2(alpha, 5e-3)
    print("alpha: ", alpha)
    S_data=[]

    for i in tqdm(range(len(volt))):
        gs3.setlevel_slowV2(volt[i], 0.005)
        #gs2.setlevel_slow(round(-1.96-0.21/1.5*(volt[i]-2),3),0.001)
        # gs2.setlevel_slow(round(0.24-0.21/1.5*(Iz+1.2),3),0.001)
        # gs1.setlevel_slow(round(-0.24-0.27/1.5*(volt[i]-2),3),0.001)
        # for j in range(len(mwf)):
        # data=vna.get_data_S21() 
        data = vna.get_data_S43()
        S_data.append(data)
        S_dataplot=np.rot90(np.array(S_data))
        S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
        Phase=np.angle(np.array(S_data))
        S_dataplot_phase=np.rot90(phase_handle(Phase, f1, 50.03))
        amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
        heatmap(np.linspace(volt[0],volt[i],i),fwr1,S_dataplot_amp,
            corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Sweep_Amp',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True,auto_open=False)
        heatmap(np.linspace(volt[0],volt[i],i),fwr1,S_dataplot_phase,
                corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep_Phase',
                filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True,auto_open=False)
        # if i!=0:
        #     fg1.clf()
        #     fg2.clf()
        # fg1=imshow(amp,figname='Two_Tone_Amp',extent=[volt[0],volt[i],fwr[0],fwr[-1]],
        #            xround='%0.2f',yround='%0.3f',xlabel=r' volt (mA)',ylabel=r'Frequency (GHz)',
        #            cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
        # fg2=imshow(phase,figname='Two_Tone_Phase',extent=[volt[0],volt[i],fwr[0],fwr[-1]],
        #            xround='%0.2f',yround='%0.3f',xlabel=r'volt (mA)',ylabel=r'Frequency (GHz)',
        #            cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    
    about=f'''
    BOX?; setup: in:port 2,out: 1;

    Attenuation: Vna={RTatt-vnapwr} dB={RTatt}+{-vnapwr} dB, inline = 60 dB = 42+8+10 dB,\n\
    Amp: outline=76 dB(RT)+36 dB(4K)
        
    VNA Paras: fstart = {fstart} GHz,fstop = {fstop} GHz,points = {points},ifband = {ifband} Hz,vnapwr = {vnapwr} dB
    MW Paras: mwpower = {mwpwr} dB, fstart = {mwfre_start} GHz, fstop = {mwfre_stop} GHz, points = {points}
    GS current: {volt[0]}mA to {volt[-1]}mA, bias type: Local
    ''' 
    # filename1=filename+filenameafter+'.npz'
    np.savez(filename1+'.npz',fre=fwr1,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
    about_txt=open(filename1+'.txt','a')
    about_txt.write(str(about))
    about_txt.close()
    
    # warning !!!!
    volt = volt[::-1]
    
mw.MW_setpower(-130)
# mw.Close_Output()
# gs3.setlevel_slow(0,0.005)
# gs6.setlevel_slowV2(0,0.005)  
##%%SAVE for two tone
Sdata=np.array(S_data)


# filename1=filename+filenameafter+'.npz'
# np.savez(filename1+'.npz',fre=fwr1,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
# about_txt=open(filename1+'.txt','a')
# about_txt.write(str(about))
# about_txt.close()

filename2=filename1
# imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt[0],volt[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
# plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
# plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

# imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt[0],volt[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
# plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
# plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)
#%%
heatmap(volt,fwr1,S_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title="Sweep_amp",filename=filename1,errobar=[],zmin=-10,zmax=10,zauto=True)
#-np.mean(S_dataplot_amp,axis=0)[:,None].repeat(721,axis=1).T heatmap(volt,fwr1,S_dataplot_phase,corlor_scale='hot_r',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=filenameafter,filename=filename,errobar=[],zmin=-1.5,zmax=2.5 ,zauto=False)

        #%%初始化网分 twotone 2
vna=keysight_vna('vna2',trace=43)
vna.two_tone_vna()
fstart=6.449#12.97#6.471#6.483#6.479#6.47
fstop=6.449#12.97#6.471#6.483#6.479#6.47
points=1001
RTatt=60
ifband=30#每秒扫描次数
vnapwr=-20#-12
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
#%%初始化微波源 twotone
mw=MW('mw1')
mw.MW_setpower(-130)
mwpwr=-50#-35

mwfre_start=2
mwfre_stop=12
points=points
mw.VnaTrig(mwfre_start,mwfre_stop,points)
mw.MW_setpower(mwpwr)
mw.start_Output()

fwr=np.linspace(mwfre_start,mwfre_stop,points)
fwr1=fwr[::-1]
# mw.Close_Output()
#%% twotone sweep mwpower
twotonenum =1#0.015#-0.016 #-0.0308#0.382;
# gs5.setlevel_slowV2(-0.2,0.0001)
# gs5.setlevel_slowV2(Global,0.0001)
which='S43'
# gs3.setlevel_slowV2(1.5,0.005)
alpha=0#gs3.getlevel()*1000

filename = f'/home/machine1/E/Data/LSY_new/data/CELv1/twotone_energyspectrum(alpha={alpha}mA)_power/two_tone_qubit_freq_{mwfre_start}GHz_to_{mwfre_stop}GHz_cavity_freq_{fstart}GHz{time.strftime("%m%d")}/'
if not os.path.exists(filename):
    os.makedirs(filename)


name0=f"two_tone_qubit_freq_{mwfre_start}GHz_to_{mwfre_stop}GHz_cavity_freq_{fstart}GHz"
filename1=filename+name0+time.strftime("%H%M%S")+'_'
exec('bg=vna.get_data_'+which+'()')
# gs1.setlevel_slowV2(-0.24,0.005)
# volt=np.linspace(-2.4,-2.4,1)
# gs6.Start_OutPut() #右z
# time.sleep(339)
gs2.setlevel_slowV2(-2.4,0.005)#-0.9,0.05
S_data=[]
i=0
mw.MW_setpower(mwpwr)
mw.Close_Output()
mwpwr_list = np.linspace(-50,0,101)
print(time.strftime("%H:%M:%S"))
for i in tqdm(range(len(mwpwr_list))):

    mw.MW_setpower(mwpwr_list[i])
    mw.start_Output()
    exec('data=vna.get_data_'+which+'()')
    mw.Close_Output()
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    Phase=np.angle(np.array(S_data))
    S_dataplot_phase=np.rot90(phase_handle(Phase, f1, 50.11))
    amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
    heatmap(np.linspace(mwpwr_list[0],mwpwr_list[i],i),fwr1,S_dataplot_amp,
        corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Sweep_Amp',
        filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True, auto_open=False)
    heatmap(np.linspace(mwpwr_list[0],mwpwr_list[i],i),fwr1,S_dataplot_phase,
            corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep_Phase',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True, auto_open=False)

mw.MW_setpower(-130)
# mw.Close_Output()
# gs3.setlevel_slow(0,0.005)
# gs6.setlevel_slowV2(0,0.005)  
##%%SAVE for two tone
Sdata=np.array(S_data)
about=' BOX?; setup: in:port 2,out: 1; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    mwpower from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,mwpwr_list[0],mwpwr_list[-1])


# filename1=filename+filenameafter+'.npz'
np.savez(filename1+'.npz',fre=fwr1,mwpwr_list=mwpwr_list,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename1 +'.txt','a')
about_txt.write(str(about))
about_txt.close()

filenameafter = "sweep"

filename2=filename1
imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[mwpwr_list[0],mwpwr_list[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.2f',xlabel=r'power(dB)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[mwpwr_list[0],mwpwr_list[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.2f',xlabel=r'power(dB)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)

heatmap(mwpwr_list,fwr1,S_dataplot_amp,corlor_scale='jet',xlabel='power(dB)',ylabel='Frequency[GHz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(mwpwr_list,fwr1,S_dataplot_phase,corlor_scale='jet',xlabel='power(dB)',ylabel='Frequency[GHz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)

#%%初始化网分 twotone 3
vna=keysight_vna('vna2',trace=43)
vna.two_tone_vna()
fstart=6.450#12.97#6.471#6.483#6.479#6.47
fstop=6.450#12.97#6.471#6.483#6.479#6.47
points=1001
RTatt=60
ifband=20#每秒扫描次数
vnapwr=0#-12
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
#%%初始化微波源 twotone
mw=MW('mw1')
mw.MW_setpower(-130)
mwpwr=-25#-35

mwfre_start=2
mwfre_stop=12
points=points
mw.VnaTrig(mwfre_start,mwfre_stop,points)
mw.MW_setpower(mwpwr)
mw.start_Output()

fwr=np.linspace(mwfre_start,mwfre_stop,points)
fwr1=fwr[::-1]
# mw.Close_Output()
#%% twotone sweep mwpower
twotonenum =1#0.015#-0.016 #-0.0308#0.382;
# gs5.setlevel_slowV2(-0.2,0.0001)
# gs5.setlevel_slowV2(Global,0.0001)
which='S43'
# gs3.setlevel_slowV2(1.5,0.005)
alpha=0#gs3.getlevel()*1000

filename = f'/home/machine1/E/Data/LSY_new/data/CELv1/twotone_energyspectrum(alpha={alpha}mA)_vna_power/two_tone_qubit_freq_{mwfre_start}GHz_to_{mwfre_stop}GHz_cavity_freq_{fstart}GHz{time.strftime("%m%d")}/'
if not os.path.exists(filename):
    os.makedirs(filename)


name0=f"two_tone_qubit_freq_{mwfre_start}GHz_to_{mwfre_stop}GHz_cavity_freq_{fstart}GHz"
filename1=filename+name0+time.strftime("%H%M%S")+'_'
exec('bg=vna.get_data_'+which+'()')
# gs1.setlevel_slowV2(-0.24,0.005)
# volt=np.linspace(-2.4,-2.4,1)
# gs6.Start_OutPut() #右z
# time.sleep(339)
gs2.setlevel_slowV2(-2.4,0.005)#-0.9,0.05
S_data=[]
i=0
mw.MW_setpower(mwpwr)
mw.Close_Output()
vnapwr_list = np.linspace(-30,6,37)
print(time.strftime("%H:%M:%S"))
for i in tqdm(range(len(vnapwr_list))):

    vna.set_power(vnapwr_list[i])
    exec('data=vna.get_data_'+which+'()')
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    Phase=np.angle(np.array(S_data))
    S_dataplot_phase=np.rot90(phase_handle(Phase, f1, 50.11))
    amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
    heatmap(np.linspace(vnapwr_list[0],vnapwr_list[i],i),fwr1,S_dataplot_amp,
        corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Sweep_Amp',
        filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True, auto_open=False)
    heatmap(np.linspace(vnapwr_list[0],vnapwr_list[i],i),fwr1,S_dataplot_phase,
            corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep_Phase',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True, auto_open=False)

mw.MW_setpower(-130)
# mw.Close_Output()
# gs3.setlevel_slow(0,0.005)
# gs6.setlevel_slowV2(0,0.005)  
##%%SAVE for two tone
Sdata=np.array(S_data)
about=' BOX?; setup: in:port 2,out: 1; \n\
    Attenuation: Vna={}dB={}+{}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
    fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
    mwpower from{}mA to{}mA'.format(RTatt-vnapwr,RTatt,-vnapwr,fstart,fstop,points,ifband,vnapwr,vnapwr_list[0],vnapwr_list[-1])


# filename1=filename+filenameafter+'.npz'
np.savez(filename1+'.npz',fre=fwr1,vnapwr_list=vnapwr_list,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename1 +'.txt','a')
about_txt.write(str(about))
about_txt.close()

filenameafter = "sweep"

filename2=filename1
imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[vnapwr_list[0],vnapwr_list[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.2f',xlabel=r'power(dB)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[vnapwr_list[0],vnapwr_list[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.2f',xlabel=r'power(dB)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)

heatmap(vnapwr_list,fwr1,S_dataplot_amp,corlor_scale='jet',xlabel='power(dB)',ylabel='Frequency[GHz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(vnapwr_list,fwr1,S_dataplot_phase,corlor_scale='jet',xlabel='power(dB)',ylabel='Frequency[GHz]',title=filenameafter,filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)

#%%初始化网分 twotone 4
vna=keysight_vna('vna2',trace=43)
vna.two_tone_vna()
fstart=6.4504#12.97#6.471#6.483#6.479#6.47
fstop=6.4504#12.97#6.471#6.483#6.479#6.47
points=101
RTatt=60
ifband=20#每秒扫描次数
vnapwr=-5#-12
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
#%%初始化微波源 twotone
mw=MW('mw1')
mw.MW_setpower(-130)
mwpwr=-13 #35

mwfre_start=3
mwfre_stop=6.5
points=points
mw.VnaTrig(mwfre_start,mwfre_stop,points)
mw.MW_setpower(mwpwr)
mw.start_Output()

fwr=np.linspace(mwfre_start,mwfre_stop,points)
fwr1=fwr[::-1]
# mw.Close_Output()
#%% twotone sweep global
# gs7=GS_200('DC7') 
# gs7.setCURRmode()
# gs7.Start_OutPut()

twotonenum =1#0.015#-0.016 #-0.0308#0.382;
gs1.setlevel_slowV2(0,0.005)
gs2.setlevel_slowV2(-2.04,0.005)
which='S43'
# gs3.setlevel_slowV2(1.5,0.005)
alpha=0#gs3.getlevel()*1000
# filename = f'/home/machine1/E/Data/LSY_new/data/test_tanmo/twotone_energyspectrum(alpha={alpha}mA)_global/two_tone_qubit_freq_{mwfre_start}GHz_to_{mwfre_stop}GHz_cavity_freq_{fstart}GHz{time.strftime("%m%d")}/'
filename = f'/home/machine1/E/Data/LSY_new/data/CELv1/twotone_energyspectrum(alpha={alpha}mA)_global/two_tone_qubit_freq_{mwfre_start}GHz_to_{mwfre_stop}GHz_cavity_freq_{fstart}GHz{time.strftime("%m%d")}/'
if not os.path.exists(filename):
    os.makedirs(filename)


name0=f"two_tone_qubit_freq_{mwfre_start}GHz_to_{mwfre_stop}GHz_cavity_freq_{fstart}GHz"
filename1=filename+name0+time.strftime("%H%M%S")+'_'
# exec('bg=vna.get_data_'+which+'()')

S_data=[]
i=0
mw.MW_setpower(mwpwr)
mw.start_Output()
volt_list = np.linspace(-0.76, -0.72, 51)
print(time.strftime("%H:%M:%S"))
gs4.setlevel_slowV2(volt_list[0], 0.0001)
for i in tqdm(range(len(volt_list))):

    gs4.setlevel_slowV2(volt_list[i], 0.0001)
    exec('data=vna.get_data_'+which+'()')
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    Phase=np.angle(np.array(S_data))
    S_dataplot_phase=np.rot90(phase_handle(Phase, f1, 67.06))
    amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
    heatmap(np.linspace(volt_list[0],volt_list[i],i),fwr1,S_dataplot_amp,
        corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Sweep_Amp',
        filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True, auto_open=False)
    heatmap(np.linspace(volt_list[0],volt_list[i],i),fwr1,S_dataplot_phase,
            corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep_Phase',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True, auto_open=False)

mw.MW_setpower(-130)
# mw.Close_Output()
# gs3.setlevel_slow(0,0.005)
# gs6.setlevel_slowV2(0,0.005)  
##%%SAVE for two tone
Sdata=np.array(S_data)    
about=f'''
BOX?; setup: in:port 2,out: 1;

Attenuation: Vna={RTatt-vnapwr} dB={RTatt}+{-vnapwr} dB, inline = 60 dB = 42+8+10 dB,\n\
Amp: outline=76 dB(RT)+36 dB(4K)
    
VNA Paras: fstart = {fstart} GHz,fstop = {fstop} GHz,points = {points},ifband = {ifband} Hz,vnapwr = {vnapwr} dB
MW Paras: mwpower = {mwpwr} dB, fstart = {mwfre_start} GHz, fstop = {mwfre_stop} GHz, points = {points}
GS current: {volt_list[0]}mA to {volt_list[-1]}mA, bias type: Global
''' 


# filename1=filename+filenameafter+'.npz'
np.savez(filename1+'.npz',fre=fwr1,volt_list=volt_list,Sdata=Sdata,vnapwr=vnapwr,about=about)
about_txt=open(filename1 +'.txt','a')
about_txt.write(str(about))
about_txt.close()

filenameafter = "sweep"

filename2=filename1
imshow(S_dataplot_amp,figname=filename2+"_AMP",extent=[volt_list[0],volt_list[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.2f',xlabel=r'power(dB)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_AMP"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_AMP"+'.jpg',format='jpg',dpi=300)

imshow(S_dataplot_phase,figname=filename2+"_PHASE",extent=[volt_list[0],volt_list[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.2f',xlabel=r'power(dB)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase(rad)",color='bwr',interpolations='None') 
plt.savefig(filename2+"_PHASE"+'.pdf',format='pdf',dpi=300)
plt.savefig(filename2+"_PHASE"+'.jpg',format='jpg',dpi=300)
#%%
heatmap(volt_list,fwr1,S_dataplot_amp,corlor_scale='hot',xlabel='power(dB)',ylabel='Frequency[GHz]',title=filenameafter,filename=filename,errobar=[],zmin=-80,zmax=-65,zauto=False)
heatmap(volt_list,fwr1,S_dataplot_phase,corlor_scale='hot',xlabel='power(dB)',ylabel='Frequency[GHz]',title=filenameafter,filename=filename,errobar=[],zmin=-1,zmax=2,zauto=False)
    

#%%compensate


Currentvalue=Currentvalue_list[i]
gs_global.setlevel_slow(Currentvalue,0.0001)

def compensate_current(test,fwr,final_point,gs_qubit_local):
    test_amp=(20*np.log10(np.abs(test))).tolist()
    min_position = test_amp.index(min(test_amp))
    freq = fwr[min_position]
    
    delta0 = final_point - freq
    
    if delta0 > 0.01 :
        
        
        for i in range(10):
            qubit_local_current = round(float(gs_qubit_local.Inst.query(':SOUR:LEV?'))*1e3,5)
            gs_qubit_local.setlevel_slow(qubit_local_current + i * 0.005, Volt[3])
            data1 = vna.getdata()
            data1_amp=(20*np.log10(np.abs(data1))).tolist()
            min_position = data1_amp.index(min(data1_amp))
            freq1 = fwr[min_position]
         
        




if Fit_condition==True:
    fg1 =plt.figure('test')
    origindata = todB(test)
    fitdata = Lorentz_fitting_VNA(f,origindata)
    Width = round(half_value_width(f, fitdata)*1e3,3)
    plt.plot(f,origindata)
    plt.plot(f,fitdata)
    plt.xlabel('Frequency(Ghz)')
    plt.ylabel('Amplitude(dB)')
    plt.title('VNA Result_Width{}MHz'.format(Width))

    plt.savefig(filename+f'I_{Currentvalue}mA_Width_{Width}MHz_{time.strftime("%m%d%H%M%S")}'+'.png',format='png',dpi=600)
else:
    fg1 =plt.figure('test')
    plt.plot(f,todB(test))    
    plt.xlabel('Frequency(Ghz)')
    plt.ylabel('Amplitude(dB)')
    plt.title('VNA Result')
    plt.savefig(filename+f'I_{Currentvalue}mA_{time.strftime("%m%d%H%M%S")}'+'.png',format='png',dpi=600)

#%% two tone measure

def Vna_setting(vnadev,vnafre,vnapwr,vnapoint,vnaifband):
    vna=keysight_vna(f'{vnadev[0]}',trace=vnadev[1])
    vna.two_tone_vna()
    vna.set_power(vnapwr)
    vna.set_startstopFre(vnafre,vnafre)
    vna.set_points_band(vnapoint,vnaifband)
    return vna

def MW_setting(mwdev,mwpwr,mwfre,vnapoint):
    mw = MW(f'{mwdev}')
    mw.MW_setpower(-130)
    mw.VnaTrig(mwfre[0],mwfre[1],vnapoint)
    mw.start_Output()
    return mw

def two_tone_measure(DC_dev, phase_delay, Volt, filename, name0, mwfre, mwpwr,vnadev,mwdev):

        vna = Vna_setting(vnadev = vnadev,vnafre=vnafre,vnapwr=vnapwr,vnapoint=vnapoint,vnaifband=vnaifband)
        mw = MW_setting(vnapoint=vnapoint,mwdev = mwdev,mwpwr=mwpwr,mwfre=mwfre)
        which='S21'        
    
        filename1 = filename + name0 + time.strftime("%H%M%S") + '_'
        
        exec('bg=vna.get_data_'+which+'()')
    
        volt=np.linspace(Volt[0],Volt[1],Volt[2])
    
        DC_dev.setlevel_slow(volt[0],0.005)#-0.9,0.05
        S_data=[]
        i=0
        mw.MW_setpower(mwpwr)
        fwr=np.linspace(mwfre[0],mwfre[1],vnapoint)
        fwr1=fwr[::-1]
        f1=np.linspace(vnafre,vnafre,vnapoint)
        for i in tqdm(range(len(volt))):
            
            DC_dev.setlevel_slow(volt[i], Volt[3])
            data=vna.get_data() 
            S_data.append(data)
            S_dataplot=np.rot90(np.array(S_data))
            S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
            Phase=np.angle(np.array(S_data))
            S_dataplot_phase=np.rot90(phase_handle(Phase, f1, phase_delay))
            amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
            heatmap(np.linspace(volt[0],volt[i],i),fwr1,S_dataplot_amp,
                corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Amp',
                filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
            heatmap(np.linspace(volt[0],volt[i],i),fwr1,S_dataplot_phase,
                    corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Phase',
                    filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        mw.MW_setpower(-130)
    
        Sdata=np.array(S_data)
        about=f' BOX?; setup: in:port 9,out: 4; \n\
            Attenuation: Vna={RTatt-vnapwr}dB={RTatt}+{-vnapwr}dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
            fstart={mwfre[0]}Ghz,fstop={mwfre[-1]}Ghz,points={vnapoint},ifband={vnaifband}Hz,vnapwr={vnapwr}dB \n\
            current from{volt[0]}mA to{volt[-1]}mA'
    
        # filename1=filename+filenameafter+'.npz'
        np.savez(filename1+'.npz',fre=fwr1,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
        about_txt=open(filename1+'.txt','a')
        about_txt.write(str(about))
        about_txt.close()
    
        filename2 = filename1
        imshow(S_dataplot_amp,figname=filename2+"_Amp",extent=[volt[0],volt[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency(GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
        plt.savefig(filename2+"Amp"+'.pdf',format='pdf',dpi=300)
        plt.savefig(filename2+"Amp"+'.jpg',format='jpg',dpi=300)
        plt.clf()
        
        imshow(S_dataplot_phase,figname=filename2+"_Phase",extent=[volt[0],volt[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.2f',xlabel=r'Current(mA)',ylabel=r'Frequency(GHz)',cbarlabel=r"Phase",color='bwr',interpolations='None') 
        plt.savefig(filename2+"Phase"+'.pdf',format='pdf',dpi=300)
        plt.savefig(filename2+"Phase"+'.jpg',format='jpg',dpi=300)
        plt.clf()
        
        heatmap(volt,fwr1,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title="Amp",filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        heatmap(volt,fwr1,S_dataplot_phase,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title="Phase",filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    
    
    
#%% measure program

RTatt = 50                  #室温衰减
vnadev = ["vna1",21]        #VNA设备
vnafre = 6.5688             #VNA频率 -> Two tone cavity frequency
vnapwr_list = [-20]                 #VNA功率
vnapoint = 201             #VNA点数
vnaifband = 5               #VNA带宽

mwdev = "mw3"               #MW设备
mwpwr_list = [-30]                 #MW功率
mwfre = [4,6]            #MW扫频范围 -> qubit frquency scan 

# Volt = [-0.16, 0, 17, 0.002] # DC current [起始电流，终止电流，点数，调整速度]
Volt = [-6, -2, 41, 0.005] # DC current [起始电流，终止电流，点数，调整速度]
# Volt = [0.5, -0.5, 51, 0.005] # DC current [起始电流，终止电流，点数，调整速度]
DC_dev = gs_qubit_local    #bias current

gs_gmon_local.setlevel_slow(0,0.005)
# gs_qubit_local.setlevel_slow(0.13,0.001)
# two_tone_setting(vnadev = vnadev,vnafre=vnafre,vnapwr=vnapwr,vnapoint=vnapoint,vnaifband=vnaifband,mwdev=mwdev,mwpwr=mwpwr,mwfre=mwfre)
for mwpwr in mwpwr_list:
    for vnapwr in vnapwr_list:
        filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_230609C/twotone_energyspectrum_VS_qubit_bias/{time.strftime("%m%d")}_two_tone_qubit_freq_{mwfre[0]}GHz_to_{mwfre[1]}GHz_cavity_freq_{vnafre}GHz/'
        if not os.path.exists(filename):
            os.makedirs(filename)
        
        filenameafter=f'twotone_gs_qubit_local_gmon_{0}mA_RTatt{RTatt}_vnapwr{vnapwr}_BW{vnaifband}_mwpwr{mwpwr}_'
        
        two_tone_measure(DC_dev, phase_delay = 50, Volt=Volt, filename=filename, name0=filenameafter, mwfre=mwfre, mwpwr=mwpwr,vnadev = vnadev, mwdev = mwdev)
#%%

RTatt = 70                  #室温衰减
vnadev = ["vna2",21]        #VNA设备
vnafre = 6.845             #VNA频率 -> Two tone cavity frequency
vnapwr = -10                 #VNA功率
vnapoint = 201              #VNA点数
vnaifband = 30               #VNA带宽

mwdev = "mw1"               #MW设备
mwpwr = -25                 #MW功率
mwfre = [6.5, 7.5]            #MW扫频范围 -> qubit frquency scan 

Volt = [1, 0, 101, 0.005] # DC current [起始电流，终止电流，点数，调整速度]
# Volt = [0, -0.16, 17, 0.002] # DC current [起始电流，终止电流，点数，调整速度]
print(np.linspace(Volt[0],Volt[1],Volt[2]))

DC_dev = gs_global    #bias current

gs_qubit_local.setlevel_slow(0,0.001)
filename = f'/home/machine1/E/Data/YingHu/gmon_coupler_2309C_XS01/twotone_energyspectrum/{time.strftime("%m%d")}_two_tone_qubit_freq_{mwfre[0]}GHz_to_{mwfre[1]}GHz_cavity_freq_{vnafre}GHz/'
if not os.path.exists(filename):
    os.makedirs(filename)

filenameafter=f'twotone_gs_gmon_local_gs_qubit_-0.1mA_RTatt{RTatt}_vnapwr{vnapwr}_BW{vnaifband}_mwpwr{mwpwr}_'
two_tone_measure(DC_dev, phase_delay = 50, Volt=Volt, filename=filename, name0=filenameafter, mwfre=mwfre, mwpwr=mwpwr,vnadev = vnadev, mwdev= mwdev)
#%%
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_local_Couplequbit'
filename1=filename+'Globalsweep{}_Att_{}_vnapwr_{}_mwpwr{}_Cinch3_1_GS_03_'.format(twotonenum,RTatt,vnapwr)+str(time.strftime("%H%M%S"))+'.npz'
np.savez(filename1,fre=f,volt=volt[0:i+1],Sdata=Sdata,bg=bg,vnapwr=vnapwr,about=about,fwr=fwr,mwpwr=mwpwr)


# #filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
# np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr,about=about)
# filename2=filename+'Globalsweep{}_Att_{}_vnapwr_{}_Cinch3_1_GS_03_1'.format(twotonenum,RTatt,vnapwr)+str(time.strftime("%H%M%S"))
# plt.savefig(filename2+'.pdf',format='pdf',dpi=600)
# plt.savefig(filename2+'.jpg',format='jpg',dpi=600)
# heatmap(volt,f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Globalsweep{}_Att_{}_vnapwr_{}_Cinch3_1_GS_03_1to2transition'.format(twotonenum,RTatt,vnapwr)+str(time.strftime("%H%M%S")),filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)




#savetxt(filename2,about)
fg1=imshow(amp,figname='Two_Tone_Amp',
           extent=[volt[0],volt[i],fwr[0],fwr[-1]],xround='%0.2f',yround='%0.3f',
           xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',
           cbarlabel=r"Amplitude(dB)",color='jet',interpolations='None') 
plt.savefig(filename1+'sweep_Amp'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
fg2=imshow(phase,figname='Two_Tone_Phase',
           extent=[volt[0],volt[i],fwr[0],fwr[-1]],xround='%0.2f',yround='%0.3f',
           xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase(rad)",
           color='bwr',interpolations='None')
plt.savefig(filename1+'sweep_Phase'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
heatmap(np.linspace(volt[0],volt[i],i+1),fwr1,amp,
        corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Two_Tone_Amp',
        filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(np.linspace(volt[0],volt[i],i+1),fwr1,phase,
        corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Two_Tone_Phase',
        filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)

#%% SA-setup (resonance fluorescence)
sa=SA('sa1')
mw=MW('mw3')
Currentvalue=0.16
centerfre=5.22
fstart,fstop,points=centerfre-0.1,centerfre+0.1,501
sa.setFrePoints(fstart, fstop, points)
sa.setResBW(2) #单位Mhz
sa.setVideoBW(100) #单位hz
sa.setyRefLevel(-52)
sa.setyScaleDiv(1)
fre,spec=sa.getdata()
plt.figure()
plt.plot(fre,spec)
#%% resonance fluorescence-measure
def listmean(list_j):
    list_j=np.array(list_j)
    for i in range(list_j.shape[0]):
        if i==0:
            data=list_j[i]
        else:
            data=data+list_j[i,:]
    return data/list_j.shape[0]

which='singlephoton'
qubitnum=1

mw.MW_setFre(centerfre)

mwpower=np.linspace(-70,-10,31)
Powerpoints=31
loopnum=500
sadata_mwoff=[]
sadata_mwon=[]
sadata_i=[]
for i in range(len(mwpower)):
    mw.MW_setpower(mwpower[i])
    sadata_j=[]
    for j in range(loopnum):
        mw.Close_Output()
        fre,data_mwoff=sa.getdata()
        sadata_mwoff.append(data_mwoff)
        mw.start_Output()
        fre,data_mwon=sa.getdata()
        sadata_mwon.append(data_mwon)
        sadata_j.append(data_mwon-data_mwoff)
        if j%(int(loopnum/10))==0:
            print('scan_Fr: {:.0f}% completed.'.format((j+1)/loopnum*100))
        
    sadata_i.append(listmean(sadata_j))
    sadataplot=np.rot90(np.array(sadata_i))
    if i!=0:
        fg1.clf()
        fg2.clf()
    fg1=imshow(sadataplot,figname='resonance fluorescence_qubit_{}_center{}Ghz'.format(qubitnum,centerfre),Figsize=(10,8),extent=[mwpower[0],mwpower[i],fstart,fstop],
                xround='%0.2f',yround='%0.3f',xlabel=r'Power(dBm)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',
                interpolations='None',vmax=0.5) 
    fg2=imshow(sadataplot,figname=which+'resonance fluorescence_qubit_{}_center{}Ghz_Origin'.format(qubitnum,centerfre),Figsize=(10,8),extent=[mwpower[0],mwpower[i],fstart,fstop],
               xround='%0.2f',yround='%0.3f',xlabel=r'Power(dBm)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',
               interpolations='None') 
    print('{}of{}'.format(i,len(mwpower)))
mw.Close_Output()

about=' BOX?; setup: in:port22,out:17; \n  inline=52dB=42+10dB, outline=76dB(RT)+36dB(4K) \n centerfre={},fstart={}Ghz,fstop={}Ghz,Frequencypoints={},Powerpoints={},loopnum={} \ current value {}mA'.format(centerfre,fstart,fstop,points,Powerpoints,loopnum,Currentvalue)+'\n qubit1 current 0.4mA(ezq),\n qubit2 current {}mA,\n qubit3 current {}mA,\n qubit4 current {}mA,\n qubit5 current {}mA,\n qubit6 current {}mA,\n'.format(round(float(gs1.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs7.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs6.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs9.Inst.query(':SOUR:LEV?'))*10e2,3),round(float(gs3.Inst.query(':SOUR:LEV?'))*10e2,3))


filename2=filename+which+'resonance fluorescence_qubit_{}_center{}Ghz'.format(qubitnum,centerfre)+str(time.strftime("%H%M%S"))+'.txt'
savetxt(filename2,about)
filename1=filename+which+'resonance fluorescence_qubit_{}_center{}Ghz'.format(qubitnum,centerfre)+str(time.strftime("%H%M%S"))+'.npz'
np.savez(filename1,fre=fre,mwpower=mwpower,sadataplot=sadataplot,sadata_mwoff=sadata_mwoff,sadata_mwon=sadata_mwon,loopnum=loopnum)
fg1.savefig(filename+which+'resonance fluorescence_qubit_{}_center{}Ghz'.format(qubitnum,centerfre)+'.pdf',format='pdf',dpi=300)
fg1.savefig(filename+which+'resonance fluorescence_qubit_{}_center{}Ghz'.format(qubitnum,centerfre)+'.png') 
fg2.savefig(filename+which+'resonance fluorescence_qubit_{}_center{}Ghz_Origin'.format(qubitnum,centerfre)+'.pdf',format='pdf',dpi=300)
fg2.savefig(filename+which+'resonance fluorescence_qubit_{}_center{}Ghz_Origin'.format(qubitnum,centerfre)+'.png') 

about_txt=open(filename+'resonance fluorescence_qubit_{}_center{}Ghz'.format(qubitnum,centerfre)+str(time.strftime("%H%M%S"))+'.txt','a')
about_txt.write(str(about))
about_txt.close()
#%% pre anti cross fit(energy spectrum fit)
#from qutip import *
import numpy as np  
import matplotlib
import matplotlib.pyplot as plt  
#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.unicode'] = True
#matplotlib.rcParams['font.family']='sans-serif'
import numpy as np  
from scipy import misc
import scipy
from scipy import ndimage
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator, MaxNLocator,LogLocator,FixedFormatter

Initialfolder='E:\\Data\\'
folder='ZhouJianYong\\6qubits\\20220721\\'
filename='Energyspetrum3_qubit1_vnapwr-30_RTatt50_Cinch3_1_GS_07_005152'
data=np.load(Initialfolder+folder+filename+'.npz')
data.files
fre = data['fre']
volt=data['volt']
Sdata=data['Sdata']
Normalized_amplitude=0.1
Magneticpoints=101;
Frequencypoints=501;
CenterFrequency=4.9;
Frequencybandwidth=1.8;
Frequency0=CenterFrequency-0.5*Frequencybandwidth;
Frequency1=CenterFrequency+0.5*Frequencybandwidth;
period=1.82;   #mA
optimal_magnetic=-225#-1759.3;
Current0=-1.25; #mA
Current1=0.5  #mA
Field0=Current0*1000/period-optimal_magnetic; #in units of milli phi_0
Field1=Current1*1000/period-optimal_magnetic; #in units of milli phi_0


Amp=abs(Sdata)
Phase=np.angle(Sdata)
Z=Amp*np.exp(1j*Phase)
Z1=Z.max()
Normalized_Z=Z/Z1
Normalized_Amp=abs(Normalized_Z)
Normalized_Phase=np.angle(Normalized_Z)


Magneticfield_list = np.linspace(Field0, Field1, Magneticpoints);
Frequency_list = np.linspace(Frequency0, Frequency1, Frequencypoints);
Frequency_mat, Magneticfield_mat=np.meshgrid(Frequency_list, Magneticfield_list)
#fig = plt.figure(figsize=(8,6),dpi=100)
fig = plt.figure()
ax1 = fig.add_axes([0.15, 0.15, 0.8, 0.8])
c=plt.imshow(Normalized_Amp.T,aspect='auto',interpolation='gaussian',origin='lower',
             extent=[Field0,Field1,Frequency0,Frequency1],
            vmax=1, vmin=0)

c.set_cmap('hot')
cbar = fig.colorbar(c, ticks=[0, 0.5,1])

cbar.ax.set_yticklabels(['0','0.5','1.0'],fontsize=18)
cbar.set_label("$|r|$",rotation='270',labelpad=20,fontsize=20)


for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_size(fontsize=20)
    # ax1.xaxis.set_major_locator(FixedLocator([-20,-10,0,10,20]))
    # ax1.xaxis.set_minor_locator(MultipleLocator(2))
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_size(fontsize=20)
    tick.label.set_size(fontsize=20)
    # ax1.yaxis.set_major_locator(FixedLocator([4.8,5.3,5.8,6.3,6.8]))
    # ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
#label_f1 = "Probing power:-65dBm\n$\Delta=1.51\,$GHz\n$I_{p}=108\,$nA\n$g=45\,MHz$"
#plt.text(0.5, 0.8,label_f1,
#     horizontalalignment='right',
#     verticalalignment='center',fontsize=16,
#     transform = ax1.transAxes)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%0.0f'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
# ax1.set_xlabel(r'Current ($\mu$A)',fontsize=18)
ax1.set_xlabel(r'$\delta\Phi/\Phi_0\times10^3$',fontsize=18)
ax1.set_ylabel(r'Frequency (GHz)',fontsize=24)
# plt.xlim(1760,1770)
# plt.ylim(6.755,6.855)

#plt.xlim(225,226.8)
#ax1.autoscale(tight=True)
Magneticfield_list = np.linspace(Field0, Field1, Magneticpoints);
Frequency_list = np.linspace(Frequency0, Frequency1, Frequencypoints);
#
 

Ec=0.01
Ej=(5.36+Ec)**2/(8*Ec)
x=np.linspace(-300,300,641)
y=(8*Ec*Ej*np.cos(x*(np.pi)/1000))**(1/2)-Ec
plt.plot(x,y,'b--',linewidth=1.5)
# plt.xlim(-500,500)
#plt.ylim(3,18)
# ax2.autoscale(tight=True)

# theoretic calculation
#Field=np.linspace(-40,40,26);
#delta=4.81;
#Ip=0.212;
#biasenergy=2*Ip*np.array(Field);
#Energy=np.sqrt(biasenergy**2+delta**2);
#plt.plot(Field,Energy, color='k',marker='.',markersize=4,linewidth='0')
# plt.savefig(filename+'reflection_phase_fit.pdf',format='pdf',dpi=600)
plt.show()

#%% anti cross fit
from qutip import *
import numpy as np  
import matplotlib
#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.unicode'] = True
#matplotlib.rcParams['font.family']='sans-serif'
import numpy as np  
from scipy import misc
import scipy
from scipy import ndimage
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator, MaxNLocator,LogLocator,FixedFormatter

Initialfolder='E:\\Data\\'
folder='ZhouJianYong\\6qubits\\20220723\\'
filename='Energyspetrum4_qubit1withqubit2_vnapwr-30_RTatt50_Cinch3_2_GS_01_180019'
data=np.load(Initialfolder+folder+filename+'.npz')
data.files
fre = data['fre']
volt=data['volt']
Sdata=data['Sdata']

Normalized_amplitude=0.1
Magneticpoints=101;
Frequencypoints=601;
CenterFrequency=4.55;
Frequencybandwidth=0.3;
Frequency0=CenterFrequency-0.5*Frequencybandwidth;
Frequency1=CenterFrequency+0.5*Frequencybandwidth;
period=1.82;   #mA
optimal_magnetic=-223#-1759.3;
Current0=-0.2; #uA
Current1=0.15  #uA
Field0=Current0*1000/period-optimal_magnetic; #in units of milli phi_0
Field1=Current1*1000/period-optimal_magnetic; #in units of milli phi_0

Amp=abs(Sdata)
Phase=np.angle(Sdata)
Z=Amp*np.exp(1j*Phase)
Z1=Z.max()
Normalized_Z=Z/Z1
Normalized_Amp=abs(Normalized_Z)
Normalized_Phase=np.angle(Normalized_Z)


Magneticfield_list = np.linspace(Field0, Field1, Magneticpoints);
Frequency_list = np.linspace(Frequency0, Frequency1, Frequencypoints);

  
import matplotlib.pyplot as plt  
  
Frequency_mat, Magneticfield_mat=np.meshgrid(Frequency_list, Magneticfield_list)
#fig = plt.figure(figsize=(8,6),dpi=100)
fig = plt.figure()
ax1 = fig.add_axes([0.15, 0.15, 0.8, 0.8])
c=plt.imshow(Normalized_Amp.T,aspect='auto',interpolation='gaussian',origin='lower',
             extent=[Field0,Field1,Frequency0,Frequency1],
            vmax=1, vmin=0)

c.set_cmap('hot')
cbar = fig.colorbar(c, ticks=[0, 0.5,1])
#cbar.set_clim(0, 1)
cbar.ax.set_yticklabels(['0','0.5','1.0'],fontsize=18)
cbar.set_label("$|r|$",rotation='270',labelpad=20,fontsize=20)
################################################
# ##############caculate the transmission spectrum with J-C model
def compute(delta, wr0):#, g0_list,g1_list):
#    # Pre-compute operators for the hamiltonian
#
#    # qubit operators
    sz1 = tensor(sigmaz(), qeye(N))
    sx1 = tensor(sigmax(), qeye(N))

    # resonator
    a=tensor(qeye(2),destroy(N))

    # preallocate output array
    evals_mat = np.zeros((len(field_list), 2*N))
    for idx, field in enumerate(field_list):
        # evaluate the Hamiltonian
        #Ip=0.1
        #eps = 2*Ip*field
        g0 = 0.04#*np.exp(np.abs(field)/80)  # atom-resonator coupling strength
        Energy=(8*360*delta*np.cos(np.pi*field/1000))**(1/2)-delta
        H = 0.5*Energy*sz1+g0*(sp2*sm1+sm2*sp1)+wr0*sz2/2
        #H = 0.5*Energy*sz1+g0*(a+a.dag())*sx1+wr0*a.dag()*a
        #H = 0.5*delta*sx1+0.5*eps*sz1+g0*(a+a.dag())*sz1+wr0*a.dag()*a+wr1*b.dag()*b+g1*(b+b.dag())*sz1
        # find the energy eigenvalues and vectors of the composite system
        evals, evecs = H.eigenstates()
        evals_mat[idx, :] = evals
    return evals_mat

delta =0.01   # atom delta frequency
field_list=np.linspace(115,305,191) # magnetic field 
wr0 = 4.56    # resonator frequency
#g1_list=0
#Ip=0.2    # persistent current
N=3
#eps_list=2*Ip*field_list

#    run computation
evals_mat = compute(delta, wr0)#, g0_list,g1_list)

for n in [1,2,3,4]:
    plt.plot(field_list, (evals_mat[:, n]- evals_mat[:, 0]), color='w',linestyle='--',linewidth=2)
  
###############################################
##############calculate with Full Hamiltonian
#import FQ_3JJ_C_res_H as flo
#
#flux_shift=1
#fs=np.linspace(flux_shift+0.483,flux_shift+0.517,171)
#band=flo.band_GHz(90,3.1,1.02,6.1,1,fs,3.1e-15,6.013,6,5)
#plt.plot(1000*(fs-flux_shift-0.5),band[1,:]-band[0,:],color='w',linestyle='--',linewidth='2')
#plt.plot(1000*(fs-flux_shift-0.5),band[2,:]-band[0,:],color='w',linestyle='--',linewidth='2')
#plt.plot(1000*(fs-flux_shift-0.5),band[3,:]-band[0,:],color='w',linestyle='--',linewidth='2')
#plt.ylim(Frequency0,Frequency1)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_size(fontsize=16)
    ax1.xaxis.set_major_locator(FixedLocator([150,200,250,300]))
    # ax1.xaxis.set_minor_locator(MultipleLocator(2))
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_size(fontsize=20)
    # tick.label.set_size(fontsize=20)
    ax1.yaxis.set_major_locator(FixedLocator([4.4,4.5,4.6,4.7]))
    # ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
#label_f1 = "Probing power:-65dBm\n$\Delta=1.51\,$GHz\n$I_{p}=108\,$nA\n$g=45\,MHz$"
#plt.text(0.5, 0.8,label_f1,
#     horizontalalignment='right',
#     verticalalignment='center',fontsize=16,
#     transform = ax1.transAxes)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%0.0f'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
# ax1.set_xlabel(r'Current ($\mu$A)',fontsize=18)
ax1.set_xlabel(r'$\delta\Phi/\Phi_0\times10^3$',fontsize=18)
ax1.set_ylabel(r'Frequency (GHz)',fontsize=24)
# plt.xlim(-380,-378.5)
plt.ylim(4.4,4.7)
# plt.savefig(filename+'reflection_anticross_amp_fit.pdf',format='pdf',dpi=600)
#plt.xlim(225,226.8)
# ax1.autoscale(tight=True)
Magneticfield_list = np.linspace(Field0, Field1, Magneticpoints);
Frequency_list = np.linspace(Frequency0, Frequency1, Frequencypoints);
#
#  

##import numpy as np  
#
#


# theoretic calculation
#Field=np.linspace(-40,40,26);
#delta=4.81;
#Ip=0.212;
#biasenergy=2*Ip*np.array(Field);
#Energy=np.sqrt(biasenergy**2+delta**2);
#plt.plot(Field,Energy, color='k',marker='.',markersize=4,linewidth='0')
# plt.savefig(filename+'reflection_anticross_amp.pdf',format='pdf',dpi=600)
plt.show()

#%%qubit linewidth fit
import numpy as np
import matplotlib.pyplot as plt

data = np.load('E:\\Data\\ZhouJianYong\\6qubits\\20220719\\Energyspetrumlinewidth_qubit4_vnapwr-30_RTatt50_Cinch3_3_GS_01_210851.npz')
data.files
data_amp = data['data_amp']
data_bg = data['data_bg']
fre = data['fre']
NOR_power = (abs(data_amp)/abs(data_bg))**2


fit = np.loadtxt('C:\\Users\\peng\\Desktop\\powerfit.txt')

plt.figure()
plt.plot(fre,NOR_power,'k.',markersize=5)
plt.plot(fre,fit,'r-',lw=2)


#%%qubit linewidth fit


#from qutip import *
import numpy as np  
from scipy import misc
import scipy
from scipy import ndimage
from matplotlib import ticker
#from skimage.filter import tv_denoise

Initialfolder='E:\\Data\\'
folder='ZhouJianYong\\6qubits\\20220729\\'
filename1='ReflectionVSPower_qubit2_5.202_Att_60_signal110108'
filename2='ReflectionVSPower_qubit2_5.202_Att_60_noise111957'
data1=np.load(Initialfolder+folder+filename1+'.npz')
data2=np.load(Initialfolder+folder+filename2+'.npz')
data1.files
data2.files
Magneticpoints = data1['Magneticpoints']
Frequencypoints = data1['Frequencypoints']
CenterFrequency= data1['CenterFrequency']
Frequencybandwidth= data1['Frequencybandwidth']
Sdata1=data1['Sdata']
Sdata2=data2['Sdata']
Frequency0=CenterFrequency-0.5*Frequencybandwidth;
Frequency1=CenterFrequency+0.5*Frequencybandwidth;


Amp1=abs(Sdata1)
Phase1=np.angle(Sdata1)
Amp2=abs(Sdata2)
Phase2=np.angle(Sdata2)

from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator, MaxNLocator,LogLocator,FixedFormatter
import matplotlib.pyplot as plt  
#import numpy as np  
#class FixedOrderFormatter(ScalarFormatter):
#    """Formats axis ticks using scientific notation with a constant order of 
#    magnitude"""
#    def __init__(self, order_of_mag=0, useOffset=True, useMathText=False):
#        self._order_of_mag = order_of_mag
#        ScalarFormatter.__init__(self, useOffset=useOffset, 
#                                 useMathText=useMathText)
#    def _set_orderOfMagnitude(self, range):
#        """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
#        self.orderOfMagnitude = self._order_of_mag
Frequency=np.linspace(Frequency0,Frequency1,Frequencypoints);
Power1=Amp1**2;
# Z1=Power1.max()
Power2=Amp2**2;
# Normalized_Z=Z/Z1
column=1
#Average_Power=Power.sum(axis=0)/Magneticpoints
Normalized_Power=Power1[column]/Power2[column];
fig = plt.figure()
ax = fig.add_axes([0.18, 0.16, 0.76, 0.8])
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_size(fontsize=20)
    # ax.xaxis.set_major_locator(FixedLocator([5.20,5.24,5.28,5.32]))
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_size(fontsize=20)
    #ax.yaxis.set_major_locator(FixedLocator([0.92,0.94,0.96,0.98,1.00]))
 # Force the y-axis ticks to use 1e-9 as a base exponent 
#ax.yaxis.set_major_formatter(FixedOrderFormatter(-3))

# Make the x-axis ticks formatted to 0 decimal places 
ax.xaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
ax.set_xlabel(r'Frequency (GHz)',fontsize=20)
ax.set_ylabel(r'Normalized $|T|^2$',fontsize=20)

# theoretic calculation

# plt.plot(Frequency,Normalized_Power, color='k',marker='.',markersize=10,linewidth=0)
plt.plot(Frequency,Normalized_Power,'k.',markersize=18,linewidth=1)
ax.autoscale(tight=True)
#plt.xlim(5.958,5.998)
#plt.ylim(-0.01,1.01)
from scipy.optimize import leastsq # Levenberg-Marquadt Algorithm #
import pylab
# lorentz line for constructing the more complex function
Initial_width=0.001
#########################################################################
########################### DEFINING FUNCTIONS ##########################
y=Normalized_Power
def lorentzian(Frequency,p):
#    numerator =  (p[0]**2 )
    denominator = 4*( Frequency - (p[1]) )**2 + p[0]**2
    y = 2*p[2]/np.pi*(p[0]/denominator)
    return y

def residuals(p,y,Frequency):
    err = y - lorentzian(Frequency,p)
    return err

#########################################################################
######################## BACKGROUND SUBTRACTION #########################

# defining the 'background' part of the spectrum #
ind_bg_low = (Frequency > min(Frequency)) & (Frequency < (CenterFrequency-0.1*Frequencybandwidth))
ind_bg_high = (Frequency > (CenterFrequency+0.1*Frequencybandwidth)) & (Frequency < max(Frequency))

Frequency_bg = np.concatenate((Frequency[ind_bg_low],Frequency[ind_bg_high]))
y_bg = np.concatenate((y[ind_bg_low],y[ind_bg_high]))
#pylab.plot(x_bg,y_bg)

# fitting the background to a line # 
m, c = np.polyfit(Frequency_bg, y_bg, 1)

# removing fitted background # 
background = m*Frequency + c
y_bg_corr = y - background

p = [Initial_width,Frequency[y.argmin()],y.min()]  # [hwhm, peak center, intensity] #

# optimization # 
pbest = leastsq(residuals,p,args=(y_bg_corr,Frequency),full_output=1)
best_parameters = pbest[0]

# fit to data #
fit = lorentzian(Frequency,best_parameters)
pylab.plot(Frequency,fit+background,'r-',lw=2)
label_f1 = "$f_{0}=7.777\,$GHz\n $\kappa=3.38\,$MHz\n Probing power:-22dBm"
plt.text(0.54, 0.2,label_f1,
      horizontalalignment='left',
      verticalalignment='center',fontsize=17,
      transform = ax.transAxes)
plt.savefig(Initialfolder+folder+filename1+'linewidth at 7.777GHz.pdf',format='pdf',dpi=600)
plt.show()
#%%



