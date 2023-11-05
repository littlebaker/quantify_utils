  # -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 20:04:53 2021

@author: peng
"""
#%%
import time
from tqdm import tqdm, trange
from hunanu.Instrument.GS200_SPT  import GS_200
from hunanu.Instrument.VNA_keysightN5231b_4Portsand2Ports import keysight_vna
import matplotlib.pyplot as plt
import numpy as np
from def_plotly import heatmap
import matplotlib
from hunanu.Instrument.SMB_100A import MW
matplotlib.rcParams['text.usetex'] = True
# matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['font.family']='sans-serif'
from hunanu.Functions.function0 import savetxt
# from hunanu.Functions.function0 import savetxt,plot_removebg,todB,imshow_rxh,phase_handle
def phase_handle(Phase,f1,phase_delay):
    # phase_delay=67.86#62.03
    phase_offset=0
    phase_correct=np.zeros((Phase.shape[0],Phase.shape[1]))
    for i in range(len(f1)):
        phase_temp=Phase[:,i]+2*np.pi*f1[i]*phase_delay+phase_offset
        phase_correct[:,i]=np.angle(np.cos(phase_temp)+1j*np.sin(phase_temp))
    return phase_correct
def imshow(Amp,figname=None,Figsize=(10,8),extent=[0,1,0,1],xround='%0.0f',yround='%0.1f',xlabel=r'$\delta\Phi/\Phi_0\times10^3$',ylabel=r'Frequency (GHz)',cbarlabel=r"$|r|$",color='bwr',interpolations=None,vmax=None):
    if figname!=None:
        fig = plt.figure(figname,figsize=Figsize)
    else:
        fig = plt.figure(figsize=(10,8))
    ax1 = fig.add_axes([0.15, 0.15, 0.8, 0.8])
    if vmax==None:
        c=plt.imshow(Amp,aspect='auto',interpolation=interpolations,origin='upper',
                     extent=extent,
                    vmax=np.max(Amp), vmin=np.min(Amp))
    else:
        c=plt.imshow(Amp,aspect='auto',interpolation=interpolations,origin='upper',
             extent=extent,
            vmax=vmax, vmin=np.min(Amp))
    c.set_cmap(color)
#    c.set_cmap('jet')
    cbar = fig.colorbar(c)
    cbar.ax.tick_params(labelsize=18)
    cbar.set_label(cbarlabel,rotation='270',labelpad=20,fontsize=20)
    
    for tick in ax1.xaxis.get_major_ticks():
        tick.label.set_size(fontsize=20)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_size(fontsize=20)

#    ax1.xaxis.set_major_formatter(FormatStrFormatter(xround))
#    ax1.yaxis.set_major_formatter(FormatStrFormatter(yround))

    ax1.set_xlabel(xlabel,fontsize=24)
    ax1.set_ylabel(ylabel,fontsize=24)
    plt.pause(0.1)
    plt.show()
    return 
def plot22(f,Sdata21,Sdata43,Sdata41,Sdata23,name, show=True):
    fig=plt.figure('single',figsize=(8,8))
    ax1=fig.add_subplot(1,2,1);ax2=fig.add_subplot(1,2,2)
    plt.title(name)    
    ax1.plot(f[::-1],20*np.log10(np.abs(Sdata21)),label='S21',color='b',marker='+',markersize=2,linewidth=0.1)
    ax1.plot(f[::-1],20*np.log10(np.abs(Sdata43)),label='S43',color='g',marker='o',markersize=2,linewidth=0.1)
    ax1.legend(loc='upper right')
    
    ax2.plot(f[::-1],20*np.log10(np.abs(Sdata41)),label='S41',color='b',marker='+',markersize=2,linewidth=0.1)
    ax2.plot(f[::-1],20*np.log10(np.abs(Sdata23)),label='S23',color='g',marker='o',markersize=2,linewidth=0.1)
    ax2.legend(loc='upper right' )
    plt.pause(0.1)
    if show:
        plt.show()
    return fig
#%%VNA
vna=keysight_vna('vna2',trace=4)
fstart=6.45-0.025
fstop=6.45+0.025
points=51
ifband=20
vnapwr=-15
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# vna.close()
#%%Global preset
vna=keysight_vna('vna2',trace=41)
fstart=6.5
fstop=7
points=251
ifband=30
vnapwr=-5
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# gs4=GS_200('DC4')
# gs4.setCURRmode()
# gs4.Start_OutPut()
# gs5.setlevel_slow(-0.3274,0.0001)
# gs4.getlevel()
#%%SAVE_preset
# gs5=GS_200('DC5')
import time
import os
filename = f'/home/machine1/E/Data/LSY_new/data/CEL_v2/LQ_cavity_VS_Global current_2D/'
if not os.path.exists(filename):
    os.makedirs(filename)
name0='cav_fsart{}_fstop{}_vnapwr{}_BW{}'.format(fstart,fstop,vnapwr,ifband)
#%%扫global
name=name0+time.strftime("%H%M%S")+'_'
filename1=filename+name
# gs5.setCURRmode()
gs7.getlevel()
Globalnum=0

# time.sleep(500)
volt=np.linspace(-1,1,2001)
# gs5.Start_OutPut()
# gs5.setlevel_slowV2(-0.2,0.0001)
gs7.setlevel_slowV2(volt[0],0.0002)
S_data=[];S23_data=[];
S43_data=[];S41_data=[]
i=0
fig = plt.figure()
print(time.strftime("%H:%M:%S"))
for i in tqdm(range(len(volt))):
    gs7.setlevel_slowV2(volt[i],0.0001)
    # data=vna.get_data()
    Sdata21, Sdata43, Sdata41, Sdata23 = vna.get_data()
    # Sdata21, Sdata43, Sdata41, Sdata23 = vna.get_data_all()  
    # data2=vna.get_data_S43()
    S_data.append(Sdata21);
    # S23_data.append(Sdata23);S41_data.append(Sdata41)
    # S43_data.append(data2)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    Phase=np.angle(np.array(S_data))
    S_dataplot_phase=np.rot90(phase_handle(Phase,f1,50.17))
    
    S43_data.append(Sdata43)
    S43_dataplot=np.rot90(np.array(S43_data))
    S43_dataplot_amp=20*np.log10(np.abs(S43_dataplot))
    S43_Phase=np.angle(np.array(S43_data))
    S43_dataplot_phase=np.rot90(phase_handle(S43_Phase,f1,67.06))
    # S43_dataplot=np.rot90(np.array(S43_data))
    # S43_dataplot_amp=20*np.log10(np.abs(S43_dataplot))
    # Phase43=np.angle(np.array(S43_data))
    # S43_dataplot_phase=np.rot90(phase_handle(Phase43,f1,60))    
    # fg1=imshow(S_dataplot_amp,figname='JPAA80sweep_Amp{}'.format(Globalnum),extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    # fg2=imshow(S_dataplot_phase,figname='JPAA80sweep_Phase{}'.format(Globalnum),extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase(rad)",color='bwr',interpolations='None')
    
    heatmap(np.linspace(volt[0],volt[i],i+1),f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep21_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True,auto_open = False)
    heatmap(np.linspace(volt[0],volt[i],i+1),f,S_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep21_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True,auto_open = False)
    heatmap(np.linspace(volt[0],volt[i],i+1),f,S43_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep43_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True,auto_open = False)
    heatmap(np.linspace(volt[0],volt[i],i+1),f,S43_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep43_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True,auto_open = False)
    # if i!=0:
    #     fig.clf()
    # fig=plot22(f,Sdata21,Sdata43,Sdata41,Sdata23,str(round(volt[i],4))+'mA')
    # if i%(int(len(volt)/10))==0 and i!=0:
    #     print('scan_Fr: {:.0f}% completed.'.format((i+1)/len(volt)*100))
    #     print(time.strftime("%H:%M:%S"))

    # if i==0:
    #     t0=time.perf_counter()    
    # elif i==1:
    #     t1=time.perf_counter()
    #     T=(t1-t0)*len(volt)
    #     print('need {}min for one fig'.format(T/60))
print(time.strftime("%H:%M:%S"))
gs7.getlevel()
# gs5.Stop_OutPut()
print(time.strftime("%H:%M:%S"))

# filename1=filename+name
# Sdata=np.array(S_data)
# about=' BOX?; setup: in:port2,out:1; in:port24,out:4 \n Attenuation: RTatt = 60 dB,\
#     inline=60dB=42+10dB(add), outline=50dB(RT)+36dB(4K) \n\
#         fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n \
#             current from{}mA to{}mA'.format(fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])
# np.savez(filename1,fre=f[::-1],vnapwr=vnapwr,volt=volt[0:i],S21data=Sdata,S43data=S43_data,S23data=S23_data,S41data=S41_data,about=about)
#%%SAVE
def SaveFig(filename=None,name=None,HT='False',Num=2,FigData=[0,0],FigName=[0,0],Color=['',''],volt=[]):
    filename1=filename+name
    for k in range(Num):
        fg1=imshow(FigData[k],figname=FigName[k],
            extent=[volt[0],volt[-1],f[-1],f[0]],xround='%0.2f',yround='%0.3f',
            xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",
            color=Color[k],interpolations='None')
        plt.savefig(filename1+FigName[k]+'.png',format='png',dpi=600)
    if HT=='True':    
        heatmap(volt[0:i+1],f,FigData[k],
                corlor_scale=Color[k],xlabel='Current[mA]',ylabel='Frequency[Ghz]',
                title=name+FigName[k],filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
    return fg1
        
filename1=filename+name
Sdata=np.array(S_data)
about=' BOX?; setup: in:port9,out:1; \n Attenuation: Vna=40dB=40dB,\
    inline=60dB=42+8+10dB(add), outline=76dB(RT)+36dB(4K) \n\
        fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n \
            current from{}mA to{}mA'.format(fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_global.txt'
# np.savez(filename1,fre=f[::-1],vnapwr=vnapwr,volt=volt[0:i+1],S43data=S_data,about=about)
np.savez(filename1,fre=f[::-1],vnapwr=vnapwr,volt=volt[0:i],S21data=Sdata,S43data=S43_data,S23data=S23_data,S41data=S41_data,about=about)
# SaveFig(filename=filename,name=name,HT='False',Num=4,FigData=[S_dataplot_amp,S_dataplot_phase,S43_dataplot_amp,S43_dataplot_phase],FigName=['S21Fluxsweep_Amp','S21sweep_Phase','S43Fluxsweep_Amp','S43sweep_Phase'],Color=['jet','bwr','jet','bwr'],volt=volt)
SaveFig(filename=filename,name=name,HT='False',Num=2,FigData=[S_dataplot_amp,S_dataplot_phase],FigName=['S21Fluxsweep_Amp','S21sweep_Phase'],Color=['jet','bwr'],volt=volt)
# SaveFig(HT='False',Num=2,FigData=[S_dataplot_amp,S_dataplot_phase],FigName=['S43Fluxsweep_Amp','S43sweep_Phase'],Color=['jet','bwr'])

#%%DC-connect
gs2=GS_200('DC2')
gs2.setCURRmode()
gs2.Start_OutPut()
gs2.setRange_Lev(12)
gs4=GS_200('DC4')
gs4.setCURRmode()
gs4.Start_OutPut()
gs4.setRange_Lev(1.2)
# gs4.setlevel_slow(0,0.002)
gs3=GS_200('DC3')
gs3.setCURRmode()
gs3.Start_OutPut()
gs3.setRange_Lev(12)
# gs3.setlevel_slow(-1.5,0.002)
gs6=GS_200('DC6')
gs6.setCURRmode()
gs6.Start_OutPut()
gs6.setRange_Lev(12)
# gs6.setlevel_slow(0.1642,0.002)
gs7=GS_200('DC7')
gs7.setCURRmode()
gs7.setRange_Lev(1.2)
gs7.getlevel()
gs7.Start_OutPut()
gs1=GS_200('DC1')
gs1.setCURRmode()
gs1.setRange_Lev(12)
gs1.getlevel()
gs1.Start_OutPut()
#%%Preset 
gs4.setlevel_slowV2(2.5,0.005)
print(gs4.getlevel())
# time.sleep(100)
gs2.setlevel_slowV2(-3.205,0.005)
print(gs2.getlevel())
gs3.setlevel_slowV2(0.035,0.005)
print(gs3.getlevel())
gs1.setlevel_slowV2(-0.44,0.005)
print(gs1.getlevel())
gs6.setlevel_slowV2(-3.4,0.005)
print(gs6.getlevel())
gs7.setlevel_slowV2(3.05,0.005)
print(gs7.getlevel())
#%%扫偏置S21
# mw2=MW('mw2')
# mw2.MW_setpower(-130)
# mwpwr0=6#-45#-35
# fc0=6.470*2
# mw1.MW_setFre(fc0)
# mw1.MW_setpower(mwpwr0)
# mw1.start_Output()
# data2=np.load(filename+'2.5_2mA_2.npz')
# data2.files
# volt=np.array([round(j,3) for j in data2['Iz_list'][332:459]])
# fw01=data2['y01'][332:459]

# filename='E:\\Data\\RXH\\EP\\20211123\\'
# gs5.setlevel_slow(-0.2481,0.0001)
Global =-0.74#0.015#-0.016 #-0.0308#0.382;
# gs5.setlevel_slowV2(-0.2,0.0001)
# time.sleep(3600)
gs4.setlevel_slowV2(Global,0.0001)
# time.sleep(39)
# gs4.setlevel_slowV2(-1.45,0.005)#-0.4
alpha=0*1000
filename = f'/home/machine1/E/Data/LSY_new/data/CELv1/cavity_VS_RZ_loop_2D(Global={Global}mA,alpha={alpha}mA)/'
if not os.path.exists(filename):
    os.makedirs(filename)
name0='Global{}_cav_rightqubit_alpha{}_fsart{}_fstop{}_vnapwr{}_BW{}_'.format(Global,alpha,fstart,fstop,vnapwr,ifband)
name=name0+time.strftime("%H%M%S")+'_'
filename1=filename+name
# gs6.setlevel_slowV2(0,0.005)
volt=np.linspace(-3,3,121)
# volt=voltx
gs3.setlevel_slowV2(volt[0],0.005)
S_data=[]
S43_data=[]
i=0
print(time.strftime("%H:%M:%S"))
for i in tqdm(range(len(volt))):
    #tqdm
    # mw2.MW_setFre(fw01[i])
    gs3.setlevel_slowV2(volt[i],0.005)
    data=vna.get_data_S21() 
    if i!=0:
        fig.clf()
    fig=plot22(f,data,data,data,data,str(round(volt[i],4))+'mA')
    # Sdata21,Sdata43,Sdata41,Sdata23=vna.get_data_all()
    # data=Sdata21
    # data2=Sdata43
    S_data.append(data)
    # S43_data.append(data2)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    Phase=np.angle(np.array(S_data))
    S_dataplot_phase=np.rot90(phase_handle(Phase,f1,59.76)) 
    heatmap(np.linspace(volt[0],volt[i],i+1),f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=name+'sweep_Amp',filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
    heatmap(np.linspace(volt[0],volt[i],i+1),f,S_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=name+'sweep_Phase',filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
    # heatmap(np.linspace(volt[0],volt[i],i),f,S43_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep43_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    # heatmap(np.linspace(volt[0],volt[i],i),f,S43_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep43_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    if i==0:
        t0=time.perf_counter()    
    elif i==1:
        t1=time.perf_counter()
        T=(t1-t0)*len(volt)
        print('need {}min for one fig'.format(T/60))
# gs2.setlevel_slowV2(0,0.005)
# gs2.setlevel_slow(-1,0.005)#-1.86
print(time.strftime("%H:%M:%S"))
np.savez(filename1,fre=f[::-1],vnapwr=vnapwr,volt=volt[0:i+1],S21data=S_data,about=about)
SaveFig(filename=filename,name=name,HT='False',Num=2,FigData=[S_dataplot_amp,S_dataplot_phase],FigName=['S21Fluxsweep_Amp','S21sweep_Phase'],Color=['jet','bwr'],volt=volt[0:i])
# plt.plot(volt,fw01,'--')
# plt.ylim(6.39,6.53)
# bg=vna.get_data_S43()  
# fig=plot22(f,bg,bg,bg,bg,str(round(volt[i],4))+'mA')
# plt.savefig(filename+'left_cavity_bg'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
# np.savez(filename+'left_cavity_bg'+str(time.strftime("%H%M%S"))+'.npz',fre=f[::-1],vnapwr=vnapwr,bg=bg)
# plt.figure()
# volt0=-1.63
# gs2.setlevel_slow(volt0,0.005)
# data=vna.get_data_S21()
# plt.plot(f[::-1],20*np.log10(np.abs(data)),label=str(volt0))
# plt.legend()
plt.figure();i=34; plt.plot(f,S_dataplot_amp[:,i],'.-',label=str(volt[i])); plt.legend()
#%%扫偏置S43
# gs2=GS_200('DC2')
# gs2.setCURRmode()
# gs2.setRange_Lev(12)
# gs2.Start_OutPut()
# -0.75 ~ -0.73 51
# 3 - 6G 101
# 
# gs1=GS_200('DC1')
# gs1.setCURRmode()
# gs1.setRange_Lev(12)
# gs1.Start_OutPut()

vna=keysight_vna('vna2',trace=43)
fstart=6.5
fstop=6.8
points=151
ifband=30
vnapwr=-5
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# time.sleep(83)
Global =-0.756     #-0.0308#0.382;
gs4.setlevel_slowV2(Global,0.0001)
# gs4.getlevel()
volt=np.linspace(-4, 0, 401)
alpha_list=[0]
for alpha in alpha_list: 
    gs1.setlevel_slowV2(alpha,0.005) #alpha
    # alpha=gs1.getlevel()*1000
    
    filename = f'/home/machine1/E/Data/LSY_new/data/CELv1/cavity_VS_Z_loop_2D(Global={Global}mA,alpha={alpha}mA)/'
    if not os.path.exists(filename):
        os.makedirs(filename)
    name0='Global{}_cav_rightqubit_alpha{}_fsart{}_fstop{}_vnapwr{}_BW{}_'.format(Global,alpha,fstart,fstop,vnapwr,ifband)
    # time.sleep(300)
    name=name0+time.strftime("%H%M%S")+'_'
    filename1=filename+name
    # volt=np.linspace(-3, 0, 301)
    gs2.setlevel_slowV2(volt[0],0.005)
    S_data=[]
    S43_data=[]
    i=0
    print(time.strftime("%H:%M:%S"))
    for i in tqdm(range(len(volt))):
        gs2.setlevel_slow(volt[i],0.005)
        # gs2.setlevel_slow(round(-1.28-0.34/2*(volt[i]+0.78),3),0.001)
        data=vna.get_data_S43() 
        if i!=0:
            fig.clf()
            fg1.clf()
            fg2.clf()
        fig=plot22(f,data,data,data,data,str(round(volt[i],4))+'mA')
        # Sdata21,Sdata43,Sdata41,Sdata23=vna.get_data_all()
        # data=Sdata21
        # data2=Sdata43
        S_data.append(data)
        # S43_data.append(data2)
        S_dataplot=np.rot90(np.array(S_data))
        S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
        Phase=np.angle(np.array(S_data))
        S_dataplot_phase=np.rot90(phase_handle(Phase,f1,50.11))  
        fg1=imshow(S_dataplot_amp,figname='Amp_Global{}'.format(Global),Figsize=(8,6),extent=[volt[0],volt[i if i!=0 else 1],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
        fg2=imshow(S_dataplot_phase,figname='Phase_Global{}'.format(Global),Figsize=(8,6),extent=[volt[0],volt[i if i!=0 else 1],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase(rad)",color='bwr',interpolations='None')
        # if i>1:
        #     fg1=imshow(S_dataplot_amp,figname='Fluxsweep_Amp',extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
        #     fg2=imshow(S_dataplot_phase,figname='Fluxsweep_Phase',extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase(rad)",color='bwr',interpolations='None')
        # if i!=len(volt)-1 and i>1:
            # fg1.close()
            # fg2.close()   
        heatmap(np.linspace(volt[0],volt[i],i+1),f,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=name+'sweep_Amp',filename=filename,errobar=[],zmin=0,zmax=1,zauto=True, auto_open=False)
        heatmap(np.linspace(volt[0],volt[i],i+1),f,S_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=name+'sweep_Phase',filename=filename,errobar=[],zmin=0,zmax=1,zauto=True, auto_open=False)
        if i==0:
            t0=time.perf_counter()    
        elif i==1:
            t1=time.perf_counter()
            T=(t1-t0)*len(volt)
            print('need {}min for one fig'.format(T/60))
    # gs1.setlevel_slow(0,0.005)#-0.315
    # gs7.setlevel_slow(0,0.005)#-0.315
    print(time.strftime("%H:%M:%S"))
    np.savez(filename1,fre=f[::-1],vnapwr=vnapwr,volt=volt[0:i+1],S43data=S_data,about=about)
    SaveFig(filename=filename,name=name,HT='False',Num=2,FigData=[S_dataplot_amp,S_dataplot_phase],FigName=[f'S43Fluxsweep_Amp_{alpha}',f'S43sweep_Phase{alpha}'],Color=['jet','bwr'],volt=volt[0:i+1])
    volt=volt[::-1]
    print(alpha)
gs2.setlevel_slowV2(0,0.005)
#%%vna
vna=keysight_vna('vna2',trace=4)
fstart=6.45
fstop=6.55
points=501
ifband=10
vnapwr=-15
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# gs6.setlevel_slowV2(-3.5,0.005)
data,data2,Sdata41,Sdata23=vna.get_data_all()
fig=plot22(f,data,data2,Sdata41,Sdata23,str(round(Iz,4))+'mA')
#%%扫中间比特quick
Global =-0.36 #-0.0308#0.382;
gs5.setlevel_slowV2(Global,0.0001)
Iz_list=np.linspace(-3,3,201)
gs2.setlevel_slowV2(-1.0,0.005)
alpha_list=[0.00]
# time.sleep(59)
gs1.setlevel_slowV2(0,0.005)#中间alpha#-1.07
# gs2.setlevel_slowV2(Iz_list[0],0.005)#中间z
# gs3.setlevel_slowV2(alpha_list[0],0.005)#中间alpha#-1.07
gs6.setlevel_slowV2(Iz_list[0],0.005)#中间z
alpha=0#gs3.getlevel()*1000
# filename='E:\\Data\\RXH\\EP\\20220106\\'
name0='Global{}_middle_q_change_'.format(Global)
# gs4.getlevel()
i=0;j=0
ap=[str(round(alpha_list[i],2)) for i in range(len(alpha_list))]
for alpha in alpha_list:
    # gs3.setlevel_slowV2(alpha, 0.005)
    i=0
    name=name0+'alpha_'+ap[j]+'_'+time.strftime("%H%M%S")
    filename1=filename+name
    S_data=[]
    S23_data=[];S43_data=[];S41_data=[]
    print(time.strftime("%H:%M:%S"))
    # S43_data=[]
    for Iz in tqdm(Iz_list):
        gs6.setlevel_slow(Iz, 0.005) 
        # gs2.setlevel_slow(round(-1.16-0.06/0.4*(Iz+3.3),3),0.001)
        # gs1.setlevel_slow(round(-0.4-0.08/0.8*(Iz+3.3),3),0.001)
        data=vna.get_data_S21() 
        S_data.append(data);
        S_dataplot=np.rot90(np.array(S_data))
        S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
        Phase=np.angle(np.array(S_data))
        S_dataplot_phase=np.rot90(phase_handle(Phase,f1,60))
        heatmap(np.linspace(Iz_list[0],Iz_list[i],i+1),f,S_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title= name+'sweep_Phase',filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)    
        heatmap(np.linspace(Iz_list[0],Iz_list[i],i+1),f,S_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title= name+'sweep_Amp',filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
        if i!=0:
            fig.clf()
        fig=plot22(f,data,data,data,data,str(round(Iz,4))+'mA')
        if i==0:
            t0=time.perf_counter()    
        elif i==1:
            t1=time.perf_counter()
            T=(t1-t0)*len(Iz_list)
            print('need {}min for one fig'.format(T/60))  
        i=i+1
    about=' BOX?; setup: in:port9,out:1; z:gs2; alpha:gs4; \n Attenuation: Vna=60dB=30+30dB,\
     inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
         fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n \
             current from{}mA to{}mA, alpha={}mA'.format(fstart,fstop,points,ifband,vnapwr,Iz_list[0],Iz_list[i-1],alpha)
    np.savez(filename1,fre=f[::-1],vnapwr=vnapwr,volt=Iz_list[0:i],alpha=alpha,
             S21data=S_data,about=about)
    fg1=imshow(S_dataplot_amp,figname='Fluxsweep_Amp'+str(time.strftime("%H%M%S")),extent=[Iz_list[0],Iz_list[i-1],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='jet',interpolations='None') 
    plt.savefig(filename1+'_sweep_Amp'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg1)
    fg2=imshow(S_dataplot_phase,figname='Fluxsweep_Phase'+str(time.strftime("%H%M%S")),extent=[Iz_list[0],Iz_list[i-1],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase(rad)",color='bwr',interpolations='None')
    plt.savefig(filename1+'_sweep_Phase'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg2)
    print(j)
    j=j+1
    Iz_list=Iz_list[::-1]
# mw.Close_Output()
# gs2.setlevel_slowV2(0, 0.005)
# gs6.setlevel_slowV2(0, 0.005)
# gs3.setlevel_slowV2(0, 0.005)
# gs7.setlevel_slowV2(0, 0.005)
# gs1.setlevel_slowV2(0, 0.005)
# gs2.setlevel_slowV2(0, 0.005)
# gs4.setlevel_slowV2(0, 0.005)
#%%扫中间比特
# Global =-0.435 #-0.0308#0.382;
# gs5.setlevel_slow(Global,0.0001)
Iz_list=np.linspace(-3.5,-3.8,101)
# gs1.setlevel_slow(0,0.005)
alpha_list=[0.08,0.05]
gs6.setlevel_slowV2(Iz_list[0],0.005)#中间z
gs3.setlevel_slowV2(alpha_list[0],0.005)#中间alpha#-1.07
alpha=gs3.getlevel()*1000
# filename='E:\\Data\\RXH\\EP\\20220106\\'
name0='Global{}_ceter_q_change_'.format(Global)
# gs4.getlevel()
i=0;j=0
ap=[str(round(alpha_list[i],2)) for i in range(len(alpha_list))]
for alpha in alpha_list:
    gs3.setlevel_slowV2(alpha, 0.005)
    i=0
    name=name0+'alpha_'+ap[j]+'_'+time.strftime("%H%M%S")
    filename1=filename+name
    S_data=[]
    S23_data=[];S43_data=[];S41_data=[]
    print(time.strftime("%H:%M:%S"))
    # S43_data=[]
    for Iz in tqdm(Iz_list):
        gs6.setlevel_slow(Iz, 0.0005) 
        # gs2.setlevel_slow(round(-1.28-0.06/0.4*(Iz+3.58),3),0.001)
        gs2.setlevel_slow(round(-1.34-0.06/0.4*(Iz+3.4),3),0.001)
        gs1.setlevel_slow(round(-0.44-0.08/0.8*(Iz+3.4),3),0.001)
        # gs1.setlevel_slow(round(-0.44-0.08/0.8*(Iz+3.6),3),0.001)#-0.44
        # data=vna.get_data_S43() 
        data,data2,Sdata41,Sdata23=vna.get_data_all()
        S_data.append(data);
        S23_data.append(Sdata23);S41_data.append(Sdata41);S43_data.append(data2)
        S_dataplot=np.rot90(np.array(S_data))
        S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
        Phase=np.angle(np.array(S_data))
        S_dataplot_phase=np.rot90(phase_handle(Phase,f1,60))
        S43_dataplot=np.rot90(np.array(S43_data))
        S43_dataplot_amp=20*np.log10(np.abs(S43_dataplot))
        Phase43=np.angle(np.array(S43_data))
        S43_dataplot_phase=np.rot90(phase_handle(Phase43,f1,60))    
        S23_dataplot=np.rot90(np.array(S23_data)); S23_dataplot_amp=20*np.log10(np.abs(S23_dataplot))
        S41_dataplot=np.rot90(np.array(S41_data)); S41_dataplot_amp=20*np.log10(np.abs(S41_dataplot))
       
        # fg1=imshow(S_dataplot_amp,figname='JPAA80sweep_Amp{}'.format(Globalnum),extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
        # fg2=imshow(S_dataplot_phase,figname='JPAA80sweep_Phase{}'.format(Globalnum),extent=[volt[0],volt[i],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase(rad)",color='bwr',interpolations='None')
        heatmap(np.linspace(Iz_list[0],Iz_list[i],i+1),f,S_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title= name+'sweep_Amp',filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
        # heatmap(np.linspace(Iz_list[0],Iz_list[i],i+1),f,S_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=name+'sweep_Phase',filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
        heatmap(np.linspace(Iz_list[0],Iz_list[i],i+1),f,S43_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep43_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        heatmap(np.linspace(Iz_list[0],Iz_list[i],i+1),f,S23_dataplot_amp-S41_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Isolation',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        # if i%(int(len(volt)/10))==0 and i!=0:
        #     print('scan_Fr: {:.0f}% completed.'.format((i+1)/len(Iz_list)*100))
        #     print(time.strftime("%H:%M:%S"))
        if i!=0:
            fig.clf()
        fig=plot22(f,data,data2,Sdata41,Sdata23,str(round(Iz,4))+'mA')
        if i==0:
            t0=time.perf_counter()    
        elif i==1:
            t1=time.perf_counter()
            T=(t1-t0)*len(Iz_list)
            print('need {}min for one fig'.format(T/60))  
        i=i+1
    about=' BOX?; setup: in:port9,out:1; z:gs2; alpha:gs4; \n Attenuation: Vna=60dB=30+30dB,\
     inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
         fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n \
             current from{}mA to{}mA, alpha={}mA'.format(fstart,fstop,points,ifband,vnapwr,Iz_list[0],Iz_list[i-1],alpha)
    np.savez(filename1,fre=f[::-1],vnapwr=vnapwr,volt=Iz_list[0:i],alpha=alpha,
             S21data=S_data,S43data=S43_data,S23data=S23_data,S41data=S41_data,about=about)
    fg1=imshow(S_dataplot_amp,figname='S21Fluxsweep_Amp'+str(time.strftime("%H%M%S")),extent=[Iz_list[0],Iz_list[i-1],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='jet',interpolations='None') 
    plt.savefig(filename1+'_S21sweep_Amp'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg1)
    fg2=imshow(S43_dataplot_amp,figname='S43Fluxsweep_Amp'+str(time.strftime("%H%M%S")),extent=[Iz_list[0],Iz_list[i-1],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None')
    plt.savefig(filename1+'_S43sweep_Amp'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg2)
    S23_dataplot=np.rot90(np.array(S23_data)); S23_dataplot_amp=20*np.log10(np.abs(S23_dataplot))
    fg3=imshow(S23_dataplot_amp,figname='S23Fluxsweep_Amp',extent=[Iz_list[0],Iz_list[i-1],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='jet',interpolations='None') 
    plt.savefig(filename1+'_S23Fluxsweep_Amp'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg3)
    S41_dataplot=np.rot90(np.array(S41_data)); S41_dataplot_amp=20*np.log10(np.abs(S41_dataplot))
    fg4=imshow(S41_dataplot_amp,figname='S41Fluxsweep_Amp',extent=[Iz_list[0],Iz_list[i-1],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None')
    plt.savefig(filename1+'_S41Fluxsweep_Amp'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg4)
    # heatmap(Iz_list,f,S_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
    #         title=name+'S21sweep_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    # heatmap(Iz_list,f,S_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
    #         title=name+'S21sweep_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    # heatmap(Iz_list,f,S43_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
    #         title=name+'S43sweep_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    # heatmap(Iz_list,f,S43_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
    #         title=name+'S43sweep_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    print(j)
    j=j+1
    Iz_list=Iz_list[::-1]
# plt.figure()
# plt.plot(f,S43_dataplot_amp[:,109])
# time.sleep(3600)
    contents=f'current:\n\
        global gs5: {-0.36} mA;\n\
        left (alpha, z) gs4: {gs4.getlevel()*1000} mA, gs2: {gs2.getlevel()*1000} mA;\n\
        right (alpha, z) gs7: {gs7.getlevel()*1000} mA, gs1: {gs1.getlevel()*1000} mA;\n\
        middle (alpha, z) gs3: {gs3.getlevel()*1000} mA, gs6: {(Iz_list[0],Iz_list[i-1],i)} mA;\n-----------------\n\
    vna: att:{40} dB, Power: {vnapwr} dBm, BW: {ifband} Hz, f0: {fstart} GHz, f1: {fstop}GHz'
    savetxt(filename1+'.txt',contents)

gs6.setlevel_slowV2(0, 0.005)
gs3.setlevel_slowV2(0, 0.005)
gs7.setlevel_slowV2(0, 0.005)
gs1.setlevel_slowV2(0, 0.005)
gs2.setlevel_slowV2(0, 0.005)
gs4.setlevel_slowV2(0, 0.005)
#     # i=i-1
    # SaveFig(filename=filename,name=name,HT='False',Num=2,FigData=[S_dataplot_phase,S43_dataplot_phase],FigName=['S21Fluxsweep_Phase','S43sweep_Phase'],Color=['bwr','bwr'],volt=Iz_list)
#%%扫中间比特，谱分
from hunanu.Instrument.SA_keysightN9020B import keysightSA_N9020B as SA
sa=SA('sa2')
centerfre=6.477
fstart,fstop,points=centerfre-0.05,centerfre+0.05,201
fwr1=np.linspace(fstop,fstart,points)
sa.setFrePoints(fstart, fstop, points)
sa.setResBW(1) #单位Mhz
sa.setVideoBW(100) #单位hz
sa.setyRefLevel(-52)
sa.setyScaleDiv(1)
#####################################
def listmean(list_j):
    list_j = np.array(list_j)
    for i in range(list_j.shape[0]):
        if i == 0:
            data = list_j[i]
        else:
            data = data + list_j[i, :]
    return data / list_j.shape[0]

def spec_meas(sa,mw1,loopnum=60): 
    # mw1.MW_setFre(centerfre)
    # mw1.start_Output()
    sadata_j = []
    for j in range(loopnum):
        mw1.Close_Output()
        fre,data_mwoff= sa.getdata()
        mw1.start_Output()
        fre, data_mwon = sa.getdata()
        sadata_j.append((10**(data_mwon/10)*1e-3-10**(data_mwoff/10)*1e-3)/10**(dBline/10)/1e6)
    return listmean(sadata_j)
#####################################
dBline=100
about0=f'dBline{dBline}\n'+'[10**(data_mwon/10)*1e-3-10**(data_mwoff/10)*1e-3]/10**(dBline/10)/1e6'
which='S21'
gs4.setlevel_slowV2(-2.5,0.005)
lalpha=gs4.getlevel()*1000
gs2.setlevel_slowV2(0.24,0.005)
lz=gs2.getlevel()*1000
#####################################
centerfre=12.98
mwpwr0=10#-45#-35
mw1=MW('mw1')
mw1.MW_setFre(centerfre)
mw1.MW_setpower(mwpwr0)
#####################################
indesity = spec_meas(sa,mw1,loopnum=20)
plt.figure()
plt.plot(fwr1[::-1],indesity)
#%%扫中间比特，谱分,测量
####################################################
Global =0.365 #-0.0308#0.382;
gs5.setlevel_slowV2(Global,0.0001)
Iz_list=np.linspace(-1.2,-2,41)
# gs1.setlevel_slow(0,0.005)
alpha_list=[1.5]
# gs4.setlevel_slowV2(alpha_list[0],0.005)#中间alpha#-1.07
# gs2.setlevel_slowV2(Iz_list[0],0.005)#中间z
gs3.setlevel_slowV2(alpha_list[0],0.005)#中间alpha#-1.07
gs6.setlevel_slowV2(Iz_list[0],0.005)#中间z
alpha=gs3.getlevel()*1000
# filename='E:\\Data\\RXH\\EP\\20220106\\'
name0='Global{}_middle_q_change_lefta{}_lz{}'.format(Global,lalpha,lz)
# gs4.getlevel()
i=0;j=0
ap=[str(round(alpha_list[i],2)) for i in range(len(alpha_list))]
for alpha in alpha_list:
    gs3.setlevel_slowV2(alpha, 0.005)
    i=0
    name=name0+'alpha_'+ap[j]+'_'+time.strftime("%H%M%S")
    filename1=filename+name
    S_data=[]
    S23_data=[];S43_data=[];S41_data=[]
    print(time.strftime("%H:%M:%S"))
    # S43_data=[]
    for Iz in tqdm(Iz_list):
        gs6.setlevel_slow(Iz, 0.005) 
        gs2.setlevel_slow(round(0.24-0.21/1.5*(Iz+1.2),3),0.001)
        data=spec_meas(sa,mw1,loopnum=20)
        S_data.append(data);
        S_dataplot_amp=np.rot90(np.array(S_data))
        heatmap(np.linspace(Iz_list[0],Iz_list[i],i+1),fwr1,S_dataplot_amp,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title= name+'sweep_Amp',filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
        if i!=0:
            fig.clf()
        fig=plt.figure('singleplot')
        plt.plot(fwr1[::-1],data)
        plt.title(str(round(Iz,4))+'mA')
        plt.pause(0.1)
        plt.show()
        if i==0:
            t0=time.perf_counter()    
        elif i==1:
            t1=time.perf_counter()
            T=(t1-t0)*len(Iz_list)
            print('need {}min for one fig'.format(T/60))  
        i=i+1
    about=' BOX?; setup: in:port9,out:1; z:gs2; alpha:gs4; \n Attenuation: Vna=60dB=30+30dB,\
     inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
         fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n \
             current from{}mA to{}mA, alpha={}mA'.format(fstart,fstop,points,ifband,vnapwr,Iz_list[0],Iz_list[i-1],alpha)
    np.savez(filename1,fre=fwr1[::-1],vnapwr=vnapwr,volt=Iz_list[0:i],alpha=alpha,
             S21data=S_data,about=about, about0=about0)
    fg1=imshow(S_dataplot_amp,figname='Fluxsweep_Amp'+str(time.strftime("%H%M%S")),extent=[Iz_list[0],Iz_list[i-1],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='jet',interpolations='None') 
    plt.savefig(filename1+'_sweep_Amp'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    # plt.close(fg1)
    print(j)
    j=j+1
    Iz_list=Iz_list[::-1]
#%%
# gs6.setlevel_slowV2(0, 0.005)
# gs3.setlevel_slowV2(0, 0.005)
# gs7.setlevel_slowV2(0, 0.005)
# gs1.setlevel_slowV2(0, 0.005)
gs2.setlevel_slowV2(0, 0.005)
gs4.setlevel_slowV2(0, 0.005)
gs5.setlevel_slowV2(0, 0.0002)
#%%扫一幅图
# gs6.setlevel_slow(0,0.002)#-0.15#0.1642
# gs3.setlevel_slow(0,0.002)#左alpha
# gs2.setlevel_slow(-4.6,0.005)#中间z#-5.04
# gs4.setlevel_slow(-2.2,0.005)#中间alpha#-1.07
# gs1.setlevel_slow(0.8,0.002)
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
volt=np.linspace(1.32,1.32,1)
# gs1.setlevel_slow(0,0.005)
# gs2.setlevel_slow(0,0.005)

fig=plt.figure(figsize=(8,8))
ax1=fig.add_subplot(2,2,1)
ax2=fig.add_subplot(2,2,2)
ax3=fig.add_subplot(2,1,2)
colorl=['b','g','r','c','m','y','k','b','purple',
       'pink','brown']
i=0
Sdata21,Sdata43,Sdata41,Sdata23=vna.get_data_all()

# index=np.argmin(abs(f1-6.477))   
# f1=deldata(f1,index); Sdata21=deldata(Sdata21,index); Sdata41=deldata(Sdata41,index)
# f=f1[::-1]; Sdata23=deldata(Sdata23,index);Sdata43=deldata(Sdata43,index)
ax1.plot(f[::-1],NorAmp(Sdata21),label='S21',color=colorl[i],marker='+',markersize=2,linewidth=0.1)
ax1.plot(f[::-1],NorAmp(Sdata43),label='S43',color=colorl[i+1],marker='o',markersize=2,linewidth=0.1)
if i==0:
    ax1.legend(loc='upper right')

ax2.plot(f[::-1],Phase_handle(np.angle(Sdata21),f[::-1],60,0),label='S21',color=colorl[i],marker='+',markersize=2,linewidth=0.1)
ax2.plot(f[::-1],Phase_handle(np.angle(Sdata43),f[::-1],65,0),label='S43',color=colorl[i+1],marker='o',markersize=2,linewidth=0.1)
if i==0:
    ax2.legend(loc='upper right')

ax3.plot(f[::-1],NorAmp(Sdata41),label='S41',color=colorl[i],marker='+',markersize=2,linewidth=0.1)
ax3.plot(f[::-1],NorAmp(Sdata23),label='S23',color=colorl[i+1],marker='o',markersize=2,linewidth=0.1)
if i==0:
    ax3.legend(loc='upper right' )

# ax=fig.add_subplot(2,1,2)
# ax.plot(f[::-1],todB(Sdata23),label='S23')
# plt.legend(loc='upper right')

plt.savefig(filename+name0+'left_cavity'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
np.savez(filename+name0+'two_cavity'+str(time.strftime("%H%M%S"))+'.npz',fre=f[::-1],vnapwr=vnapwr,Sdata21=Sdata21,Sdata43=Sdata43,Sdata41=Sdata41,Sdata23=Sdata23,about=about)
#%%###########################################################################
Sdata21=vna.get_data() 
S_dataplot_amp=20*np.log10(np.abs(Sdata21))
#%%
Phase=np.angle(np.array(Sdata21))
S_dataplot_phase=phase_handle(Phase[None,:],f1,50.11)
#%%
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

import os
filename = f'/home/machine1/E/Data/LSY_new/data/CELv1/'
if not os.path.exists(filename):
    os.makedirs(filename)
plt.savefig(filename+'left_cavity_spectrum'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)

#%%
fig=plt.figure(figsize=(8,8))
ax1=fig.add_subplot(2,1,1);ax2=fig.add_subplot(2,1,2)
l_cut=70;r_cut=200
Amp=NorAmp(Sdata21)
A0,A,x0,kappa=fitlorentz(f1[l_cut:r_cut],Amp[l_cut:r_cut])
ax1.plot(f1,Amp,'g')
ax1.plot(f1[l_cut:r_cut],lorentz(f1[l_cut:r_cut],A0,A,x0,kappa),'r--')
ax1.text(0.1,0.9,str(round(kappa*1000,3))+'\n'+str(round(x0,4)),transform=ax1.transAxes,fontsize=20);
l_cut=120;r_cut=300
Amp=NorAmp(Sdata43)
A0,A,x0,kappa=fitlorentz(f1[l_cut:r_cut],Amp[l_cut:r_cut])
ax2.plot(f1,Amp,'g')
ax2.plot(f1[l_cut:r_cut],lorentz(f1[l_cut:r_cut],A0,A,x0,kappa),'r--')
ax2.text(0.1,0.9,str(round(kappa*1000,3))+'\n'+str(round(x0,4)),transform=ax2.transAxes,fontsize=20);
plt.savefig(filename+name0+'fit_two_cavity'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
#%%取点画图
volt=np.linspace(0.2,-1,121)
gs5.setlevel_slow(0.2,0.0001)
S21_data=[]
S43_data=[]
i=0
print(time.strftime("%H:%M:%S"))
for i in range(len(volt)):
    gs5.setlevel_slow(volt[i],0.0001)
    Sdata21=vna.get_data_S21() 
    # Sdata21,Sdata43,Sdata41,Sdata23=vna.get_data_all()
    data_amp21=np.mean(np.abs(Sdata21))
    data_pha21=np.mean(np.angle(Sdata21))
    # data_amp43=np.mean(np.abs(Sdata43))
    # data_pha43=np.mean(np.angle(Sdata43))
    S21_data.append((data_amp21,data_pha21))
    # S43_data.append((data_amp43,data_pha43))
    if i==0:
        t0=time.perf_counter()    
    elif i==1:
        t1=time.perf_counter()
        T=(t1-t0)*len(volt)
        print('need {}min for one fig'.format(T/60))
        print(time.strftime("%H:%M:%S"))
S21_data=np.array(S21_data)
# S43_data=np.array(S43_data)
fig=plt.figure()
ax=fig.add_subplot(2,2,1)
ax.plot(volt[0:i],S21_data[:,0],label='S21')
ax=fig.add_subplot(2,2,2)
ax.plot(volt[0:i],S21_data[:,1],label='S21Phase')

#%%改一边比特扫中间比特preset
def plot4(f,Sdata21,Sdata43,Sdata41,Sdata23,name):
    fig=plt.figure(figsize=(8,8))
    ax1=fig.add_subplot(2,2,1);ax2=fig.add_subplot(2,2,2);ax3=fig.add_subplot(2,1,2)
    plt.title(name)
    ax1.plot(f[::-1],NorAmp(Sdata21),label='S21',color='b',marker='+',markersize=2,linewidth=0.1)
    ax1.plot(f[::-1],NorAmp(Sdata43),label='S43',color='g',marker='o',markersize=2,linewidth=0.1)
    ax1.legend(loc='upper right')

    ax2.plot(f[::-1],Phase_handle(np.angle(Sdata21),f[::-1],60,0),label='S21',color='b',marker='+',markersize=2,linewidth=0.1)
    ax2.plot(f[::-1],Phase_handle(np.angle(Sdata43),f[::-1],65,0),label='S43',color='g',marker='o',markersize=2,linewidth=0.1)
    ax2.legend(loc='upper right')
    
    ax3.plot(f[::-1],NorAmp(Sdata41),label='S41',color='b',marker='+',markersize=2,linewidth=0.1)
    ax3.plot(f[::-1],NorAmp(Sdata23),label='S23',color='g',marker='o',markersize=2,linewidth=0.1)
    ax3.legend(loc='upper right' )
    return fig
def find_peaks(xData,yData,xstar,xstop):
    index0=int((xstar-xData[0])/(xData[1]-xData[0]))
    index1=int((xstop-xData[0])/(xData[1]-xData[0]))
    sort=np.argsort(yData[index0:index1])
    x0=sort[0]
    x0_data=xData[x0+index0]
    y0_data=yData[x0+index0]
    print(x0_data,y0_data)
    return x0_data,y0_data
gs7.setlevel_slow(2.9,0.005)#右alpha
gs3.setlevel_slow(-3.3,0.005)#中间alpha
gs2.setlevel_slow(-1.87,0.005)#左z
gs4.setlevel_slow(-0.2,0.005)#左alpha
# gs5.setlevel_slow(-0.2481,0.0001)
filename='E:\\Data\\RXH\\EP\\20211222\\'
name0='ceter_r_l_q_calphar{}_ralphal{}_lalpha{}_zl{}_BW{}_'.format(-3.3,2.45,-0.2,-1.87,ifband)
gs1.setlevel_slow(-0.49,0.005)#右z
gs6.setlevel_slow(-0.73,0.005)#中间z
Sdata21,Sdata43,Sdata41,Sdata23=vna.get_data_all()
fg=plot4(f,Sdata21,Sdata43,Sdata41,Sdata23,'-0.60 mA')
x0_data,y0_data=find_peaks(f1,Sdata21,6.46,6.48);x0_data,y0_data=find_peaks(f1,Sdata43,6.46,6.48)
#%%改一边比特扫中间比特
Izr_list=np.linspace(-0.22,-0.22,1)
Iz_list=np.linspace(-0.4,-0.8,41)
gs1.setlevel_slow(Izr_list[0],0.005)#右z
gs6.setlevel_slow(Iz_list[0],0.005)#中间z
# gs4.getlevel()
i=0;j=0
# ap=[str(round(alpha_list[i],1)) for i in range(len(alpha_list))]
for Izr in Izr_list:
    gs1.setlevel_slow(Izr, 0.005)
    i=0
    name='global_0.39_'+name0+'Izr_{}'.format(round(Izr,3))+'_'+time.strftime("%H%M%S")
    filename1=filename+name
    S_data=[]
    print(time.strftime("%H:%M:%S"))
    S_data=[];S23_data=[];
    S43_data=[];S41_data=[]
    for Iz in Iz_list:
        gs6.setlevel_slow(Iz, 0.005)
        data,data2,Sdata41,Sdata23=vna.get_data_all()
        S_data.append(data);
        S43_data.append(data2);
        S23_data.append(Sdata23); S41_data.append(Sdata41)
        S_dataplot=np.rot90(np.array(S_data))
        S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
        Phase=np.angle(np.array(S_data))
        S_dataplot_phase=np.rot90(phase_handle(Phase,f1,60))
        S43_dataplot=np.rot90(np.array(S43_data))
        S43_dataplot_amp=20*np.log10(np.abs(S43_dataplot))
        Phase43=np.angle(np.array(S43_data))
        S43_dataplot_phase=np.rot90(phase_handle(Phase43,f1,60.02)) 
        heatmap(np.linspace(Iz_list[0], Iz,i+1),f,S_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=name+'sweep_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        heatmap(np.linspace(Iz_list[0], Iz,i+1),f,S_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title=name+'sweep_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        heatmap(np.linspace(Iz_list[0], Iz,i+1),f,S43_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep43_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        heatmap(np.linspace(Iz_list[0], Iz,i+1),f,S43_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep43_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        if i%(int(len(Iz_list)/10))==0 and i!=0:
            print('scan_Fr: {:.0f}% completed.'.format((i+1)/len(Iz_list)*100))
            print(time.strftime("%H:%M:%S"))
        if i==0:
            t0=time.perf_counter()    
        elif i==1:
            t1=time.perf_counter()
            T=(t1-t0)*len(Iz_list)
            print('need {}min for one fig'.format(T/60))  
        i=i+1
    about=' BOX?; setup: in:port9,out:1; zr:gs1; zm:gs6; \n \
        Attenuation: Vna=40dB=40dB,inline=60dB=42+8+10dB(add), outline=76dB(RT)+36dB(4K) \n\
         fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n \
             current from{}mA to{}mA, Izr={}mA'.format(fstart,fstop,points,ifband,vnapwr,Iz_list[0],Iz_list[i-1],Izr)
    np.savez(filename1,fre=f[::-1],vnapwr=vnapwr,volt=Iz_list[0:i],Izr=Izr,S21data=S_data,
             S43data=S43_data,S41data=S41_data,S23data=S23_data,about=about)
    fg1=imshow(S_dataplot_amp,figname='S21Fluxsweep_Amp'+str(time.strftime("%H%M%S")),extent=[Iz_list[0],Iz_list[i-1],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='jet',interpolations='None') 
    plt.savefig(filename1+'_S21sweep_Amp'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg1)
    fg2=imshow(S_dataplot_phase,figname='S21Fluxsweep_Phase'+str(time.strftime("%H%M%S")),extent=[Iz_list[0],Iz_list[i-1],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase(rad)",color='bwr',interpolations='None')
    plt.savefig(filename1+'_S21sweep_Phase'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg2)
    fg3=imshow(S43_dataplot_amp,figname='S43Fluxsweep_Amp',extent=[Iz_list[0],Iz_list[i-1],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='jet',interpolations='None') 
    plt.savefig(filename1+'_S43sweep_Amp'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg3)
    fg4=imshow(S43_dataplot_phase,figname='S43Fluxsweep_Phase',extent=[Iz_list[0],Iz_list[i-1],f[-1],f[0]],xround='%0.2f',yround='%0.3f',xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase(rad)",color='bwr',interpolations='None')
    plt.savefig(filename1+'_S43sweep_Phase'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg4)
    # heatmap(Iz_list,f,S_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
    #         title=name+'S21sweep_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    # heatmap(Iz_list,f,S_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
    #         title=name+'S21sweep_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    # heatmap(Iz_list,f,S43_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
    #         title=name+'S43sweep_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    # heatmap(Iz_list,f,S43_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
    #         title=name+'S43sweep_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    print(j)
    j=j+1
    Iz_list=Iz_list[::-1]
#%%SA-setup
from hunanu.Instrument.SA_keysightN9020B import keysightSA_N9020B as SA
sa=SA('sa2')
# mw=MW('mw1')
centerfre=6.470
fstart,fstop,points=centerfre-0.02,centerfre+0.02,801
sa.setFrePoints(fstart, fstop, points)
sa.setResBW(3) #单位Mhz #需要小于线宽
sa.setVideoBW(100) #单位hz
sa.setyRefLevel(-52)
sa.setyScaleDiv(1)
# fre,spec=sa.getdata()
# plt.figure()
# plt.plot(fre,spec)
Iz_list=np.linspace(-0.75,-0.45,241)
gs6.setlevel_slowV2(Iz_list[0],0.005)#中间z
which='port4'
name0='Global{}_ceter_q_alpha{}_'.format(Global,alpha)
name=name0+which+'_'+time.strftime("%H%M%S")
spec_data=[]
print(time.strftime("%H:%M:%S"))
i=0
for Iz in tqdm(Iz_list):
    gs6.setlevel_slowV2(Iz, 0.005)  
    fre,spec=sa.getdata();fre=fre/1e9
    spec_data.append(spec)
    spec_dataplot=np.rot90(np.array(spec_data))
    heatmap(np.linspace(Iz_list[0],Iz_list[i],i+1),fre[::-1],spec_dataplot,corlor_scale='jet',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title= name+'Spectrum',filename=filename,errobar=[],zmin=0,zmax=1,zauto=True)
    if i!=0:
        fig.clf()
    fig=plt.figure('single',figsize=(8,8))
    plt.title(str(round(Iz,4))+'mA')
    plt.plot(fre,spec)
    plt.pause(0.1)
    plt.show()    
    if i==0:
        t0=time.perf_counter()    
    elif i==1:
        t1=time.perf_counter()
        T=(t1-t0)*len(Iz_list)
        print('need {}min for one fig'.format(T/60))  
    i=i+1
about=' BOX?; setup: in:port9,3,out:1,4; z:gs6; alpha:gs3; \n Attenuation: Vna=40dB=30+10dB,\
 inline=60dB=42+8+10dB(add), outline=76dB(RT)+36dB(4K) \n\
     fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n \
         current from{}mA to{}mA, alpha={}mA'.format(fstart,fstop,points,ifband,vnapwr,Iz_list[0],Iz_list[i-1],alpha)
np.savez(filename1,fre=fre,volt=Iz_list[0:i],alpha=alpha,spec_data=spec_data,about=about)
fg1=imshow(spec_dataplot,figname='S21Fluxsweep_Amp'+str(time.strftime("%H%M%S")),
           extent=[Iz_list[0],Iz_list[i-1],fre[0],fre[-1]],xround='%0.2f',yround='%0.3f',
           xlabel=r'Current (mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Power (dBm)",color='jet',interpolations='None') 
plt.savefig(filename1+'Spectrum'+'.png',format='png',dpi=600)
# gs6.setlevel_slow(0, 0.005)
data=np.load('E:\\Data\\RXH\\EP\\20220112\\Global0.24_ceter_q_change_alpha_4.0_175100.npz')
data.files
spec_data=data['spec_data']
fre=data['fre']
volt=data['volt']
fig=plt.figure('single',figsize=(8,8))
# plt.title(str(round(Iz,4))+'mA')
plt.plot(fre,spec_data[0,:])
plt.plot(fre,spec_data[97,:])
plt.plot(fre,spec_data[110,:])
plt.plot(np.linspace(6.47,6.47,201),np.linspace(spec_data.min(),spec_data.max(),201))
