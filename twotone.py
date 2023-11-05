# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 13:05:49 2021

@author: peng
"""
import time
from tqdm import tqdm
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
# from hunanu.Functions.function0 import savetxt,plot_removebg,todB,imshow,phase_handle
# from hunanu.Functions.function0 import savetxt,plot_removebg,todB,imshow_rxh,phase_handle
def todB(S21):
    dBS21=10*np.log10(np.abs(S21))
    return(dBS21)
# vna=keysight_vna('vna2',trace=4)#trace=4 端口四条线， trace=24 两端口四条线 trace=21 两端口一条线 S21

def imshow(Amp,figname=None,Figsize=(10,8),extent=[0,1,0,1],xround='%0.0f',yround='%0.1f',xlabel=r'$\delta\Phi/\Phi_0\times10^3$',ylabel=r'Frequency (GHz)',cbarlabel=r"$|r|$",color='bwr',interpolations=None,vmax=None,vmin=None):
    if figname!=None:
        fig = plt.figure(figname,figsize=Figsize)
    else:
        fig = plt.figure(figsize=(10,8))
    ax1 = fig.add_axes([0.15, 0.15, 0.8, 0.8])
    if vmax==None:
        vmax=np.max(Amp)
    if vmin==None:
        vmin=np.min(Amp)
    c=plt.imshow(Amp,aspect='auto',interpolation=interpolations,origin='upper',
         extent=extent,
        vmax=vmax, vmin=vmin)        
    c.set_cmap(color)
    # c.set_clim(-2,1)
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
    return fig
def plot_removebg(S_dataplot_amp,S_dataplot_phase):
    ave1=np.mean(S_dataplot_amp,axis=0)#对每一列取平均，算出每个电流下反射的平均值
    amp=[]
    for i in range(S_dataplot_amp.shape[1]):
        amp.append((S_dataplot_amp[:,i]-ave1[i]))#shape[1]指的是列数
    amp=np.array(amp)
    
    ave2=np.mean(S_dataplot_phase,axis=0)
    phase=[]
    for i in range(S_dataplot_phase.shape[1]):
        phase.append((S_dataplot_phase[:,i]-ave2[i]))
    phase=np.array(phase)
    return amp.T,phase.T
def phase_handle(Phase,f1,phase_delay=0):
    # phase_delay=0 
    phase_offset=0
    phase_correct=np.zeros((Phase.shape[0],Phase.shape[1]))
    for i in range(len(f1)):
        # phase_temp=Phase[:,i]+2*np.pi*f1[i]*phase_delay+phase_offset
        phase_temp=Phase[:,i]+2*np.pi*((f1[i]-f1[0])/(len(f1)-1)*i+f1[0])*phase_delay+phase_offset
        phase_correct[:,i]=np.angle(np.cos(phase_temp)+1j*np.sin(phase_temp))
    return phase_correct
#%%初始化电流源
# gs2=GS_200('DC2')
# gs2.setCURRmode()
gs3=GS_200('DC3')
gs3.setCURRmode()
gs3.setRange_Lev(12)
gs3.getlevel()
gs1=GS_200('DC1')
gs1.setCURRmode()
gs1.Start_OutPut()
gs1.setRange_Lev(12)
gs1.getlevel()
gs7=GS_200('DC7')
gs7.setCURRmode()
gs7.Start_OutPut()
gs7.setRange_Lev(1.2)
gs6=GS_200('DC6')
gs6.setCURRmode()
gs6.Start_OutPut()
gs6.setRange_Lev(12)
gs2=GS_200('DC2')
gs2.setCURRmode()
gs2.Start_OutPut()
gs2.setRange_Lev(12)
#%%初始化网分
vna=keysight_vna('vna2',trace=43)
vna.two_tone_vna()
fstart=6.471#12.97#6.471#6.483#6.479#6.47
fstop=6.471#12.97#6.471#6.483#6.479#6.47
points=1001
ifband=30#每秒扫描次数
vnapwr=-8#-8
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)

#%%初始化微波源
mw=MW('mw1')
mw.MW_setpower(-130)
mwpwr=-20#-20
vary=8#8
fc0=10#10
mwfre_start=fc0-vary
mwfre_stop=fc0+vary
points=points
mw.VnaTrig(mwfre_start,mwfre_stop,points)
mw.MW_setpower(mwpwr)
mw.start_Output()
fwr=np.linspace(mwfre_start,mwfre_stop,points)
fwr1=fwr[::-1]
# mw.Close_Output()
#%%
filename='E:\\Data\\RXH\\EP\\20220719\\'
#%% twotone sweep z
Global =-0.36#0.015#-0.016 #-0.0308#0.382;
# gs5.setlevel_slowV2(-0.2,0.0001)
# gs5.setlevel_slowV2(Global,0.0001)
which='S21'
# gs3.setlevel_slowV2(1.5,0.005)
alpha=0#gs3.getlevel()*1000
name0='Global_'+str(Global)+which+'twotone_middlequbit_alpha{}_fsart{}_fstop{}_vnapwr{}_BW{}_mwpwr{}_'.format(alpha,fstart,fstop,vnapwr,ifband,mwpwr)
filename1=filename+name0+time.strftime("%H%M%S")+'_'
mw.MW_setpower(-130)
exec('bg=vna.get_data_'+which+'()')
# gs1.setlevel_slowV2(-0.24,0.005)
volt=np.linspace(1.3,2.3,51)
# gs6.Start_OutPut() #右z
# time.sleep(339)
gs6.setlevel_slowV2(volt[0],0.005)#-0.9,0.05
S_data=[]
i=0
mw.MW_setpower(mwpwr)
print(time.strftime("%H:%M:%S"))
for i in tqdm(range(len(volt))):
    gs6.setlevel_slow(volt[i],0.005);#gs2.setlevel_slow(round(-1.96-0.21/1.5*(volt[i]-2),3),0.001)
    # gs2.setlevel_slow(round(0.24-0.21/1.5*(Iz+1.2),3),0.001)
    # gs1.setlevel_slow(round(-0.24-0.27/1.5*(volt[i]-2),3),0.001)
    # for j in range(len(mwf)):
    # data=vna.get_data_S21() 
    exec('data=vna.get_data_'+which+'()')
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    Phase=np.angle(np.array(S_data))
    S_dataplot_phase=np.rot90(phase_handle(Phase,fwr))
    amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
    heatmap(np.linspace(volt[0],volt[i],i),fwr1,S_dataplot_amp,
        corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Sweep_Amp',
        filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    heatmap(np.linspace(volt[0],volt[i],i),fwr1,S_dataplot_phase,
            corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep_Phase',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    # if i!=0:
    #     fg1.clf()
    #     fg2.clf()
    # fg1=imshow(amp,figname='Two_Tone_Amp',extent=[volt[0],volt[i],fwr[0],fwr[-1]],
    #            xround='%0.2f',yround='%0.3f',xlabel=r' volt (mA)',ylabel=r'Frequency (GHz)',
    #            cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    # fg2=imshow(phase,figname='Two_Tone_Phase',extent=[volt[0],volt[i],fwr[0],fwr[-1]],
    #            xround='%0.2f',yround='%0.3f',xlabel=r'volt (mA)',ylabel=r'Frequency (GHz)',
    #            cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    if i==0:
        t0=time.perf_counter()    
    elif i==1:
        t1=time.perf_counter()
        T=(t1-t0)*len(volt)
        print('need {}min for one fig'.format(T/60))  
mw.MW_setpower(-130)
# mw.Close_Output()
# gs3.setlevel_slow(0,0.005)
# gs6.setlevel_slowV2(0,0.005)  
##%%SAVE for two tone
Sdata=np.array(S_data)
about=' BOX?; setup: in:port9,3,out:1,4;\n\
    Attenuation: Attenuation: Vna=40dB=30+10dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
        fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
            current from{}mA to{}mA'.format(fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_local_Couplequbit'
np.savez(filename1,fre=f,volt=volt[0:i+1],Sdata=Sdata,bg=bg,vnapwr=vnapwr,about=about,fwr=fwr,mwpwr=mwpwr)
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
#%%AC stark
# mw=MW('mw3')
mw.MW_setpower(-130)
mwpwr=-40#-45#-35
vary=0.5
fc0=5.6
mwfre_start=fc0-vary
mwfre_stop=fc0+vary
points=points
mw.VnaTrig(mwfre_start,mwfre_stop,points)
mw.MW_setpower(mwpwr)
mw.start_Output()
fwr=np.linspace(mwfre_start,mwfre_stop,points)
fwr1=fwr[::-1]
# mw.Close_Output()
# Global =0.24#0.015#-0.016 #-0.0308#0.382;
# gs5.setlevel_slowV2(Global,0.0001)
which='S21';att='30dB'
gs4.setlevel_slowV2(0,0.005)
alpha=gs4.getlevel()*1000
volt=1.35
gs2.setlevel_slowV2(volt,0.005)#-0.9,0.05
name0='Global_'+str(Global)+which+'Acstark_leftqubit_alpha{}_att{}_f{}_mwpwr{}_volt{}'.format(alpha,att,fstart,mwpwr,volt)
filename1=filename+name0+time.strftime("%H%M%S")+'_'
vnapwr=np.linspace(-30,2,65)
S_data=[]
i=0
print(time.strftime("%H:%M:%S"))
for i in tqdm(range(len(vnapwr))):
    vna.set_power(vnapwr[i])
    # mw.MW_setpower(mwpwr[i])
    exec('data=vna.get_data_'+which+'()')
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    Phase=np.angle(np.array(S_data))
    S_dataplot_phase=np.rot90(phase_handle(Phase,fwr))
    amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
    heatmap(np.linspace(vnapwr[0],vnapwr[i],i),fwr1,S_dataplot_amp,
        corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='AcSweep_Amp',
        filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    heatmap(np.linspace(vnapwr[0],vnapwr[i],i),fwr1,S_dataplot_phase,
            corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='AcSweep_Phase',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    # if i!=0:
    #     fg1.clf()
    #     fg2.clf()
    # fg1=imshow(amp,figname='Two_Tone_Amp',extent=[volt[0],volt[i],fwr[0],fwr[-1]],
    #            xround='%0.2f',yround='%0.3f',xlabel=r' volt (mA)',ylabel=r'Frequency (GHz)',
    #            cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    # fg2=imshow(phase,figname='Two_Tone_Phase',extent=[volt[0],volt[i],fwr[0],fwr[-1]],
    #            xround='%0.2f',yround='%0.3f',xlabel=r'volt (mA)',ylabel=r'Frequency (GHz)',
    #            cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    if i==0:
        t0=time.perf_counter()    
    elif i==1:
        t1=time.perf_counter()
        T=(t1-t0)*len(vnapwr)
        print('need {}min for one fig'.format(T/60))  
mw.MW_setpower(-130)
mw.Close_Output()
# gs2.setlevel_slow(0,0.005)
# gs6.setlevel_slow(0,0.005)  
##%%SAVE for AC
Sdata=np.array(S_data)
about=' BOX?; setup: in:port9,out:1;\n\
    Attenuation: Attenuation: Vna=40dB=30+10dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
        fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,volt={}mA \n\
            current from{}dBm to{}dBm'.format(fstart,fstop,points,ifband,volt,vnapwr[0],vnapwr[-1])
#filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_local_Couplequbit'
np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,vnapwr=vnapwr[0:i+1],about=about,fwr=fwr,mwpwr=mwpwr)
#savetxt(filename2,about)
fg1=imshow(amp,figname='Ac_stark_Amp',
           extent=[vnapwr[0],vnapwr[i],fwr[0],fwr[-1]],xround='%0.2f',yround='%0.3f',
           xlabel=r'Power(dBm)',ylabel=r'Frequency (GHz)',
           cbarlabel=r"Amplitude(dB)",color='jet',interpolations='None') 
plt.savefig(filename1+'sweep_Amp'+'.png',format='png',dpi=600)
fg2=imshow(phase,figname='Ac_stark_Phase',
           extent=[vnapwr[0],vnapwr[i],fwr[0],fwr[-1]],xround='%0.2f',yround='%0.3f',
           xlabel=r'Power(dBm)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase(rad)",
           color='bwr',interpolations='None')
plt.savefig(filename1+'sweep_Phase'+'.png',format='png',dpi=600)
heatmap(np.linspace(vnapwr[0],vnapwr[i],i+1),fwr1,amp,
        corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='Ac_stark_Amp',
        filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
heatmap(np.linspace(vnapwr[0],vnapwr[i],i+1),fwr1,phase,
        corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='Ac_stark_Phase',
        filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
#%%改alpha测two_toneS21
# mw2=MW('mw2')
# mw2.MW_setpower(-130)
# mwpwr0=-130#-45#-35
# fc0=6.474
# mw2.MW_setFre(fc0)
# mw2.MW_setpower(mwpwr0)
# mw2.start_Output()

# data2=np.load(filename+'2.5mA.npz')
# data2.files
# Iz_list=data2['Iz_list']
# fw01=data2['y01']#[0:38][::-1]
which='S21'
name0='Global_'+str(Global)+which+'twotone_leftqubit_Pump{}_fsart{}_fstop{}_vnapwr{}_BW{}_mwpwr{}_'.format(mwpwr0,fstart,fstop,vnapwr,ifband,mwpwr)
filename1=filename+name0+time.strftime("%H%M%S")+'_'
Iz_list=np.linspace(1.75,2.75,51)
# time.sleep(140)
alpha_list=[-1.45]
gs2.setlevel_slowV2(Iz_list[0],0.005)#左z
gs4.setlevel_slowV2(alpha_list[0],0.005)#左alpha#-1.07
# gs2.getlevel()
i=0;j=0
ap=[str(round(alpha_list[i],3)) for i in range(len(alpha_list))]
for alpha in alpha_list:
    gs4.setlevel_slow(alpha, 0.005)
    i=0
    name=name0+'alpha_'+ap[j]+'_'+time.strftime("%H%M%S")
    filename1=filename+name
    S_data=[]
    print(time.strftime("%H:%M:%S"))
    # S43_data=[]
    mw.MW_setpower(-130)
    exec('bg=vna.get_data_'+which+'()')
    mw.MW_setpower(mwpwr)
    mw.start_Output()    
    for Iz in Iz_list:
        # mw2.MW_setFre(fw01[i])
        gs2.setlevel_slow(Iz, 0.005)
        # data2=vna.get_data_S43() 
        exec('data=vna.get_data_'+which+'()')
        S_data.append(data/bg)
        # S43_data.append(data2)
        S_dataplot=np.rot90(np.array(S_data))
        S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
        Phase=np.angle(np.array(S_data))
        S_dataplot_phase=np.rot90(phase_handle(Phase,fwr,phase_delay=0))
        amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
        heatmap(np.linspace(Iz_list[0], Iz,i+1),fwr1,S_dataplot_amp,
                corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
                title=name+'sweep_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        heatmap(np.linspace(Iz_list[0], Iz,i+1),fwr1,S_dataplot_phase,
                corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
                title=name+'sweep_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
            # heatmap(np.linspace(Iz_list[0], Iz,i),f,S43_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep43_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
            # heatmap(np.linspace(Iz_list[0], Iz,i),f,S43_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep43_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
       
        # gs4.setlevel_slowV2(1.45,0.005)#左alpha#-1.07
        # gs2.setlevel_slow(-2.015, 0.005)
        # exec('data=vna.get_data_'+which+'()')
        fig=plt.figure('single')
        plt.plot(fwr1[::-1],np.angle(data))
        plt.pause(0.1)
        plt.show()
        if i!=(len(Iz_list)-1):
            fig.clf()
        if i%(int((len(Iz_list)-1)/10))==0 and i!=0:
            print('scan_Fr: {:.0f}% completed.'.format((i+1)/len(Iz_list)*100))
            print(time.strftime("%H:%M:%S"))
        if i==0:
            t0=time.perf_counter()    
        elif i==1:
            t1=time.perf_counter()
            T=(t1-t0)*len(Iz_list)
            print('need {}min for one fig'.format(T/60))  
        i=i+1
    about=' BOX?; setup: in:port9,out:1;\n\
        Attenuation: Attenuation: Vna=30dB=30+0dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
            fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
                current from{}mA to{}mA'.format(fstart,fstop,points,ifband,vnapwr,Iz_list[0],Iz_list[-1])
    np.savez(filename1,fre=f,volt=Iz_list[0:i],Sdata=S_data,bg=bg,vnapwr=vnapwr,about=about,fwr=fwr,mwpwr=mwpwr)
    fg1=imshow((amp.T-amp[:,0].T).T,figname='S21sweep_Amp'+str(time.strftime("%H%M%S")),
               extent=[Iz_list[0],Iz_list[i-1],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',
               xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",
               color='jet',interpolations='None') 
    plt.savefig(filename1+'_S21sweep_Amp'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg1)
    fg2=imshow(phase,figname='S21sweep_Phase'+str(time.strftime("%H%M%S")),
               extent=[Iz_list[0],Iz_list[i-1],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',
               xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='jet',interpolations='None')     
    plt.savefig(filename1+'_S21sweep_Phase'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg2)
    heatmap(Iz_list[0:i],fwr1,amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
            title=name+'S21sweep_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    heatmap(Iz_list[0:i],fwr1,(phase.T).T,corlor_scale='hot',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
            title=name+'S21sweep_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    print(j)
    j=j+1
    Iz_list=Iz_list[::-1]
    mw.MW_setpower(-130)
    mw.Close_Output()
    # mw2.Close_Output()
# gs2.setlevel_slowV2(-3.205,0.005)#左z
# gs4.setlevel_slow(0,0.005)#左alpha#-1.07
#%%
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
    c.set_clim(-2,1)
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
    return fig

fg1=imshow((phase.T-phase[:,0].T).T,figname='S21sweep_Phase'+str(time.strftime("%H%M%S")),
           extent=[Iz_list[0],Iz_list[i-1],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',
           xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='jet',interpolations='None')

def find_peaks(xData,yData,xstar,xstop,form='max'):
    index0=int((xstar-xData[0])/(xData[1]-xData[0]))
    index1=int((xstop-xData[0])/(xData[1]-xData[0]))
    sort=np.argsort(yData[index0:index1])
    if form=='max':
        x0=sort[-1]
    else:
        x0=sort[0]
    x0_data=xData[x0+index0]
    y0_data=yData[x0+index0]
    print(x0_data,y0_data)
    return x0_data,y0_data
x0=[];y0=[]
for i in range(len(Iz_list)):
    if i>15 and i<36:
        x0_data,y0_data=find_peaks(fwr1,phase[:,i]-phase[:,50],fwr1[725],fwr1[900])
    else:
        x0_data,y0_data=find_peaks(fwr1,phase[:,i]-phase[:,50],fwr1[715],fwr1[725],form='min')
    x0.append(x0_data);y0.append(y0_data)
plt.figure()
plt.plot(Iz_list,x0)
xxx=Iz_list[16:21];yyy=np.array(x0[16:21])
z1 = np.polyfit(xxx, yyy, 1) # 用7次多项式拟合，可改变多项式阶数；
p1 = np.poly1d(z1) #得到多项式系数，按照阶数从高到低排列
yvals=p1(xxx)
plt.plot(xxx,yvals)
xxx=np.linspace(Iz_list[1],Iz_list[21],81)
yvals=p1(xxx+0.02)
plt.plot(xxx,yvals)
ym=np.linspace(6.472,6.472,len(xxx))
plt.plot(xxx,ym)
plt.figure()
plt.plot(fwr1,phase[:,17]-phase[:,50])

voltx=xxx;fw01=yvals
#%%change mw pwr
# mw=MW('mw3')
# mw.MW_setpower(-130)
# mwpwr=-45#-45#-35
# vary=0.25
# fc0=5
# mwfre_start=fc0-vary
# mwfre_stop=fc0+vary
# points=points
# mw.VnaTrig(mwfre_start,mwfre_stop,points)
# mw.MW_setpower(mwpwr)
# mw.start_Output()
# fwr=np.linspace(mwfre_start,mwfre_stop,points)
# fwr1=fwr[::-1]

# mw.Close_Output()
# Global =0.24#0.015#-0.0V216 #-0.0308#0.382;
# gs5.setlevel_slowV2(Global,0.0001)
which='S21';att='40dB'
gs4.setlevel_slowV2(0,0.005)
alpha=gs4.getlevel()*1000
volt=-1.08
gs2.setlevel_slowV2(volt,0.005)#-0.9,0.05
name0='Global_'+str(Global)+which+'Acstark_leftqubit_alpha{}_att{}_f{}_volt{}'.format(alpha,att,fstart,volt)
filename1=filename+name0+time.strftime("%H%M%S")+'_'
mwpwr=np.linspace(-20,0,21)
S_data=[]
i=0
print(time.strftime("%H:%M:%S"))
for i in tqdm(range(len(mwpwr))):
    # vna.set_power(vnapwr[i])
    mw.MW_setpower(mwpwr[i])
    exec('data=vna.get_data_'+which+'()')
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    Phase=np.angle(np.array(S_data))
    S_dataplot_phase=np.rot90(phase_handle(Phase,fwr))
    amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
    heatmap(np.linspace(mwpwr[0],mwpwr[i],i),fwr1,S_dataplot_amp,
        corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='AcSweep_Amp',
        filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    heatmap(np.linspace(mwpwr[0],mwpwr[i],i),fwr1,S_dataplot_phase,
            corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='AcSweep_Phase',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    # if i!=0:
    #     fg1.clf()
    #     fg2.clf()
    # fg1=imshow(amp,figname='Two_Tone_Amp',extent=[volt[0],volt[i],fwr[0],fwr[-1]],
    #            xround='%0.2f',yround='%0.3f',xlabel=r' volt (mA)',ylabel=r'Frequency (GHz)',
    #            cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    # fg2=imshow(phase,figname='Two_Tone_Phase',extent=[volt[0],volt[i],fwr[0],fwr[-1]],
    #            xround='%0.2f',yround='%0.3f',xlabel=r'volt (mA)',ylabel=r'Frequency (GHz)',
    #            cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    if i==0:
        t0=time.perf_counter()    
    elif i==1:
        t1=time.perf_counter()
        T=(t1-t0)*len(mwpwr)
        print('need {}min for one fig'.format(T/60))  
mw.MW_setpower(-130)
mw.Close_Output()
#%%sweep with spectrum analyzer_preset
mw1=MW('mw1')
mw1.MW_setpower(-130)
mwpwr0=-30#-45#-35
fc0=6.449
mw1.MW_setFre(fc0)
mw1.MW_setpower(mwpwr0)
mw1.start_Output()

# mw=MW('mw1')
# mw.MW_setpower(-130)
# mwpwr0=-50#-45#-35
# fc0=6.449
# mw.MW_setFre(fc0)
# mw.MW_setpower(mwpwr0)
# mw.start_Output()
#%%
dBline=100
mwpwr=-17#-45#-35
from hunanu.Instrument.SA_keysightN9020B import keysightSA_N9020B as SA
sa=SA('sa2')
centerfre=6.449
fstart,fstop,points=centerfre-0.00025,centerfre+0.00025,1001
fwr1=np.linspace(fstop,fstart,points)
sa.setFrePoints(fstart, fstop, points)
sa.setResBW(0.005) #单位Mhz
sa.setVideoBW(100) #单位hz
sa.setyRefLevel(-52)
sa.setyScaleDiv(1)
# gs4.setlevel_slowV2(1.5,0.005)
alpha=gs7.getlevel()*1000
centerfre=6.477*2
mw1.MW_setFre(centerfre)
mw1.MW_setpower(-130)
fre, data_mwoff = sa.getdata()

mw1.MW_setpower(mwpwr)
mw1.start_Output()
# mw=MW('mw1')
fre,spec=sa.getdata()
plt.figure()
plt.plot(fre,spec)
net_spect = (10**(spec/10)*1e-3-10**(data_mwoff/10)*1e-3)/10**(dBline/10);#dB
indesity = net_spect/1e6;
plt.figure()
plt.plot(fre,indesity)
sadata_j = []
for j in range(loopnum):
    mw1.Close_Output()
    fre,data_mwoff= sa.getdata()
    mw1.start_Output()
    fre, data_mwon = sa.getdata()
    sadata_j.append((10**(data_mwon/10)*1e-3-10**(data_mwoff/10)*1e-3)/10**(dBline/10)/1e6)
    # if j % (int(loopnum / 10)) == 0:
    #     print('scan_Fr: {:.0f}% completed.'.format((j + 1) / loopnum * 100))
fig=plt.figure('singleplot')
plt.plot(fwr1[::-1],listmean(sadata_j))
#%%sweep with spectrum analyzer
alpha = -1.5
local = 1.7
mwpower=np.array([-10])
mwfreq = np.linspace(12, 14, 51)
filename = f'/home/machine1/E/Data/LSY_new/CELv1/CEL_spectrum_RightQubit_RZ(alpha={alpha}mA)(local={local})/two_tone_qubit_freq_{mwfreq[0]}Ghz_to_{mwfreq[-1]}GHz_{mwpower[0]}db_to_{mwpower[-1]}db_{time.strftime("%m%d")}/'
if not os.path.exists(filename):
    os.makedirs(filename)

def listmean(list_j):
    list_j = np.array(list_j)
    for i in range(list_j.shape[0]):
        if i == 0:
            data = list_j[i]
        else:
            data = data + list_j[i, :]
    return data / list_j.shape[0]

#mwpower = np.power(10, np.linspace(-3, 0, 31))
# mwpower = np.linspace(-30,10,41)
loopnum = 20
dBline=100
about=f'dBline{dBline}'+'[10**(data_mwon/10)*1e-3-10**(data_mwoff/10)*1e-3]/10**(dBline/10)/1e6'
which='S21'
def spec_volt(mwpower,volt=1.63,centerfre=6.486*2,Mode='Power'): 
    # volt=1.63
    gs3.setlevel_slowV2(volt,0.001)#-0.9,0.05
    sadata_i = []; sadata_off = []; sadata_on = []
    if Mode=='Power':
         name0=which+'adder_Global_'+str(Global)+'spectrum_leftqubit_alpha{}_volt{}_att{}_pf{}_avg{}_'.format(alpha,volt,att,centerfre,loopnum)
         mw1.MW_setFre(centerfre)
         mw1.start_Output()
         mwpower=mwpower
    elif Mode=='Fre':
         name0=which+'adder_Global_'+str(Global)+'spectrum_leftqubit_alpha{}_volt{}_att{}_ppower{}_avg{}_'.format(alpha,volt,att,mwpower,loopnum)
         mw1.MW_setFre(centerfre[0])
         mw1.MW_setpower(mwpower)
         mw1.start_Output()
         mwpower=centerfre        
    filename1=filename+name0+time.strftime("%H%M%S")+'_'     
    for i in range(len(mwpower)):
        # start_time = time.time()
        if Mode=='Power':
            mw1.MW_setpower(mwpower[i])
            xlabel='Power(dBm)'
        elif Mode=='Fre':
            mw1.MW_setFre(mwpower[i])
            xlabel='Pump Fre(MHz)'
        sadata_j = []; sadata_j_on = []; sadata_j_off = []
        for j in range(loopnum):
            mw1.Close_Output()
            fre,data_mwoff= sa.getdata()
            mw1.start_Output()
            fre, data_mwon = sa.getdata()
            sadata_j.append((10**(data_mwon/10)*1e-3-10**(data_mwoff/10)*1e-3)/10**(dBline/10)/1e6)
            sadata_j_on.append(data_mwon); sadata_j_off.append(data_mwoff);
            # if j % (int(loopnum / 10)) == 0:
            #     print('scan_Fr: {:.0f}% completed.'.format((j + 1) / loopnum * 100))
        sadata_i.append(listmean(sadata_j))
        sadata_on.append(listmean(sadata_j_on))
        sadata_off.append(listmean(sadata_j_off))
        sadataplot = np.rot90(np.array(sadata_i))
        sadataplot_on = np.rot90(np.array(sadata_on))
        sadataplot_off = np.rot90(np.array(sadata_off))        
        fig=plt.figure('singleplot')
        plt.plot(fwr1[::-1],listmean(sadata_j))
        plt.pause(0.1)
        plt.show()
        if i!=(len(mwpower)-1):
            fig.clf()
        if i==0:
            t0=time.perf_counter()    
        elif i==1:
            t1=time.perf_counter()
            T=(t1-t0)*len(mwpower)
            print('need {}min for one fig'.format(T/60))  
        # plt.plot(listmean(sadata_j))
        heatmap(np.linspace(mwpower[0],mwpower[i],i+1),fwr1,sadataplot_on,
                corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='Amp_on',
                filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        heatmap(np.linspace(mwpower[0],mwpower[i],i+1),fwr1,sadataplot,
            corlor_scale='bwr',xlabel=xlabel,ylabel='Frequency(GHz)',title='Sweep_Amp',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        # print("time consume: ", time.time()-start_time)
    fg1=imshow(sadataplot,figname='sweep_Amp'+str(time.strftime("%H%M%S")),
            extent=[mwpower[0],mwpower[-1],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',
            xlabel=xlabel,ylabel='Frequency (GHz)',cbarlabel=r"S (W/Hz)",
            color='jet',interpolations='None') 
    plt.savefig(filename1+'_sweep_Amp'+'.png',format='png',dpi=600)
    # plt.close(fg1)
    fg2=imshow(sadataplot_on,figname='Amp_on'+str(time.strftime("%H%M%S")),
            extent=[mwpower[0],mwpower[-1],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',
            xlabel=xlabel,ylabel='Frequency(GHz)',cbarlabel=r"Power(dBm)",
            color='jet',interpolations='None') 
    plt.savefig(filename1+'_Amp_on'+'.png',format='png',dpi=600)
    # plt.close(fg2)
    mw1.Close_Output()
    np.savez(filename1,volt=volt,fre=fwr1,mwpower=mwpower,sadataplot=sadataplot,sadataplot_on=sadataplot_on, sadataplot_off=sadataplot_off,vnapwr=vnapwr,centerfre=centerfre,about=about)

# centerfre=6.488*2
# vlist=np.linspace(0.24,0.25,1)
# frec=np.linspace(6.475,6.479,21)*2
# mwp=[-7]
# for k in range(len(mwp)):
#     spec_volt(mwpower=mwp[k],volt=-2.57,centerfre=frec,Mode='Fre')#volt-2.625
#     print(mwp[k])


freq = [6.449]
for k in range(len(freq)):
    spec_volt(-10,volt=1.7,centerfre=mwfreq, Mode="Fre")
    print(frec[k])
    
# for k in range(len(frec)):
#     spec_volt(volt=-2.575,centerfre=frec[k])
#     print(frec[k])   
# gs4.setlevel_slowV2(0,0.005)
# gs2.setlevel_slowV2(0,0.005)

# data=np.load(filename+'S21Global_-0.36spectrum_leftqubit_alpha1.5_volt-2.57_att0dB_pf12.954_avg20_192257_.npz')
# data.files
# sadataplot_on=data['sadataplot_on'];
# sadataplot=data['sadataplot'];
# plt.figure()
# plt.plot(fwr1,sadataplot[:,10])
#%%Sweep with seed light
## mw用来加种子光； mw1用来加pump
# mwpower = np.linspace(-30,10,41)
loopnum = 20
dBline=100
about=f'dBline{dBline}'+'[10**(data_mwon/10)*1e-3-10**(data_mwoff/10)*1e-3]/10**(dBline/10)/1e6'
which='S21'
mwpwr0=-4#-45#-35
fc0=12.954
mw1.MW_setFre(fc0)
mw1.MW_setpower(mwpwr0)
mw1.Close_Output()
def spec_volt(mwpower,volt=1.63,centerfre=6.486*2,Mode='Power'): 
    # volt=1.63
    gs2.setlevel_slowV2(volt,0.001)#-0.9,0.05
    sadata_i = []; sadata_off = []; sadata_on = []
    if Mode=='Power':
         name0='pf{}_pp{}'.format(fc0,mwpwr0)+which+'Global_'+str(Global)+'spectrum_leftqubit_alpha{}_volt{}_att{}_sf{}_avg{}_'.format(alpha,volt,att,centerfre,loopnum)
         mw.MW_setFre(centerfre)
         mw.start_Output()
         mwpower=mwpower
    elif Mode=='Fre':
         name0='without pump'+which+'Global_'+str(Global)+'spectrum_leftqubit_alpha{}_volt{}_att{}_spower{}_avg{}_'.format(alpha,volt,att,mwpower,loopnum)
         mw.MW_setFre(centerfre[0])  #'pf{}_pp{}'.format(fc0,mwpwr0)
         mw.MW_setpower(mwpower)
         mw.start_Output()
         mwpower=centerfre        
    filename1=filename+name0+time.strftime("%H%M%S")+'_'     
    for i in range(len(mwpower)):
        # start_time = time.time()
        if Mode=='Power':
            mw.MW_setpower(mwpower[i])
            xlabel='Power(dBm)'
        elif Mode=='Fre':
            mw.MW_setFre(mwpower[i])
            xlabel='Fre(GHz)'
        sadata_j = []; sadata_j_on = []; sadata_j_off = []
        for j in range(loopnum):
            # mw1.Close_Output()
            mw.Close_Output()
            fre,data_mwoff= sa.getdata()
            # mw1.start_Output()
            mw.start_Output()
            fre, data_mwon = sa.getdata()
            sadata_j.append((10**(data_mwon/10)*1e-3-10**(data_mwoff/10)*1e-3)/10**(dBline/10)/1e6)
            sadata_j_on.append(data_mwon); sadata_j_off.append(data_mwoff);
            # if j % (int(loopnum / 10)) == 0:
            #     print('scan_Fr: {:.0f}% completed.'.format((j + 1) / loopnum * 100))
        sadata_i.append(listmean(sadata_j))
        sadata_on.append(listmean(sadata_j_on))
        sadata_off.append(listmean(sadata_j_off))
        sadataplot = np.rot90(np.array(sadata_i))
        sadataplot_on = np.rot90(np.array(sadata_on))
        sadataplot_off = np.rot90(np.array(sadata_off))
        fig=plt.figure('singleplot')
        plt.plot(fwr1[::-1],listmean(sadata_j))
        plt.pause(0.1)
        plt.show()
        if i!=(len(mwpower)-1):
            fig.clf()
        if i==0:
            t0=time.perf_counter()    
        elif i==1:
            t1=time.perf_counter()
            T=(t1-t0)*len(mwpower)
            print('need {}min for one fig'.format(T/60))  
        # plt.plot(listmean(sadata_j))
        heatmap(np.linspace(mwpower[0],mwpower[i],i+1),fwr1,sadataplot_on,
                corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='Amp_on',
                filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        heatmap(np.linspace(mwpower[0],mwpower[i],i+1),fwr1,sadataplot,
            corlor_scale='bwr',xlabel=xlabel,ylabel='Frequency(GHz)',title='Sweep_Amp',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        # print("time consume: ", time.time()-start_time)
    fg1=imshow(sadataplot,figname='sweep_Amp'+str(time.strftime("%H%M%S")),
            extent=[mwpower[0],mwpower[-1],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',
            xlabel=xlabel,ylabel='Frequency (GHz)',cbarlabel=r"S (W/Hz)",
            color='jet',interpolations='None') 
    plt.savefig(filename1+'_S21sweep_Amp'+'.png',format='png',dpi=600)
    # plt.close(fg1)
    fg2=imshow(sadataplot_on,figname='Amp_on'+str(time.strftime("%H%M%S")),
            extent=[mwpower[0],mwpower[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',
            xlabel=xlabel,ylabel='Frequency(GHz)',cbarlabel=r"Power(dBm)",
            color='jet',interpolations='None') 
    plt.savefig(filename1+'_Amp_on'+'.png',format='png',dpi=600)
    # plt.close(fg2)
    mw1.Close_Output()
    mw.Close_Output()
    np.savez(filename1,volt=volt,fre=fwr1,mwpower=mwpower,sadataplot=sadataplot,sadataplot_on=sadataplot_on, sadataplot_off=sadataplot_off,
             vnapwr=vnapwr,centerfre=centerfre,about=about)

# centerfre=6.488*2
# vlist=np.linspace(0.24,0.25,1)
frec=np.linspace(6.450,6.495,46)
# time.sleep(15)
mwp=[-60]
for k in range(len(mwp)):
    spec_volt(mwpower=mwp[k],volt=-2.57,centerfre=frec,Mode='Fre')#volt-2.625
    print(mwp[k])
# mwpower=np.linspace(-70,-30,21)
# frec=[6.47]
# for k in range(len(frec)):
#     spec_volt(mwpower,volt=-2.625,centerfre=frec[k])
#     print(frec[k])
#%%sweep with vna
vna=keysight_vna('vna2',trace=21)
fstart=6.44
fstop=6.50
points=301
ifband=10
vnapwr=-25
fwr1=np.linspace(fstop,fstart,points)   
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
#########################################
sadata_i = []
sadata_i_gain = []
gs4.setlevel_slowV2(2.5,0.005)
alpha=gs4.getlevel()*1000
volt=-3.205#0.135
# time.sleep(119)
gs2.setlevel_slowV2(volt,0.001)
# mw1=MW('mw3')
centerfre=6.4730*2#12.954
mwpwr=-45#-45#-35
mw1.MW_setFre(centerfre)
mw1.MW_setpower(mwpwr)
mw1.start_Output()
mw1.MW_setpower(-130)
#########################################
data_mwoff = vna.get_data_S21()
mwpower = np.linspace(-25,15,41)#np.linspace(-35,15,26)
att=40
name0='Global_'+str(Global)+which+'vna_sweepGain_leftqubit_alpha{}_volt{}_vatt{}_f{}_'.format(alpha,volt,att,centerfre)
filename1=filename+name0+time.strftime("%H%M%S")+'_'
for i in range(len(mwpower)):
    mw1.MW_setpower(mwpower[i])
    data_mwon = vna.get_data_S21()
    sadata_i.append(20*np.log10(np.abs(data_mwon)))
    sadata_i_gain.append(20*np.log10(np.abs(data_mwon/data_mwoff)))
    sadataplot = np.rot90(np.array(sadata_i))
    sadataplot_gain = np.rot90(np.array(sadata_i_gain))
    if i==0:
        t0=time.perf_counter()  
    elif i==1:
        t1=time.perf_counter()
        T=(t1-t0)*len(mwpower)
        print('need {}min for one fig'.format(T/60))  
    fig=plt.figure('singleplot')
    plt.plot(fwr1[::-1],10*np.log10(np.abs(data_mwon)))
    plt.pause(0.1)
    plt.show()
    if i!=(len(mwpower)-1):
        fig.clf()
    heatmap(np.linspace(mwpower[0],mwpower[i],i+1),fwr1,sadataplot_gain,
            corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='Amp_gain',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    heatmap(np.linspace(mwpower[0],mwpower[i],i+1),fwr1,sadataplot,
        corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[GHz]',title='Sweep_Amp',
        filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    # print("time consume: ", time.time()-start_time)
fg1=imshow(sadataplot,figname='sweep_Amp'+str(time.strftime("%H%M%S")),
        extent=[mwpower[0],mwpower[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',
        xlabel='Power(dBm)',ylabel='Frequency(GHz)',cbarlabel=r"Amp(dBm)",
        color='jet',interpolations='None') 
plt.savefig(filename1+'_S21sweep_Amp'+'.png',format='png',dpi=600)
plt.close(fg1)
fg2=imshow(sadataplot_gain,figname='Amp_gain'+str(time.strftime("%H%M%S")),
        extent=[mwpower[0],mwpower[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',
        xlabel='Power(dBm)',ylabel='Frequency(GHz)',cbarlabel=r"Gain(dB)",
        color='jet',interpolations='None') 
plt.savefig(filename1+'_Gain'+'.png',format='png',dpi=600)
plt.close(fg2)
mw1.Close_Output()
np.savez(filename1,volt=volt,fre=fwr1,mwpwr=mwpower,sadataplot=sadataplot,sadataplot_gain=sadataplot_gain,data_mwoff=data_mwoff,vnapwr=vnapwr,about=about)
#%%
mwpwr=2#-45#-35
centerfre=12.946#12.954
mw.MW_setFre(centerfre)
mw.MW_setpower(-130)
mw.start_Output()
sadata_i = []
sadata_i_gain = []
data_mwoff = vna.get_data_S21()
mw.MW_setpower(mwpwr)
mwfre = np.linspace(12.5,13.5,51)#np.linspace(-35,15,26)
name0='Global_'+str(Global)+which+'vna_sweepGain_leftqubit_alpha{}_volt{}_att{}_f{}_'.format(alpha,volt,att,centerfre)
filename1=filename+name0+time.strftime("%H%M%S")+'_'
for i in range(len(mwfre)):
    mw.MW_setFre(mwfre[i])
    mw.MW_setpower(mwpwr)
    data_mwon = vna.get_data_S21()
    sadata_i.append(20*np.log10(np.abs(data_mwon)))
    sadata_i_gain.append(20*np.log10(np.abs(data_mwon/data_mwoff)))
    sadataplot = np.rot90(np.array(sadata_i))
    sadataplot_gain = np.rot90(np.array(sadata_i_gain))
    if i==0:
        t0=time.perf_counter()  
    elif i==1:
        t1=time.perf_counter()
        T=(t1-t0)*len(mwfre)
        print('need {}min for one fig'.format(T/60))  
    fig=plt.figure('singleplot')
    plt.title(str(mwfre[i]))
    plt.plot(fwr1[::-1],10*np.log10(np.abs(data_mwon)))
    plt.pause(0.1)
    plt.show()
    mw.MW_setpower(-130)
    time.sleep(0.001)
    if i!=(len(mwpower)-1):
        fig.clf()
    heatmap(np.linspace(mwfre[0],mwfre[i],i+1),fwr1,sadataplot_gain,
            corlor_scale='bwr',xlabel='PFre[dBm]',ylabel='Frequency[Ghz]',title='Amp_gain',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    heatmap(np.linspace(mwfre[0],mwfre[i],i+1),fwr1,sadataplot,
        corlor_scale='bwr',xlabel='PFre[dBm]',ylabel='Frequency[GHz]',title='Sweep_Amp',
        filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    # print("time consume: ", time.time()-start_time)
fg1=imshow(sadataplot,figname='sweep_Amp'+str(time.strftime("%H%M%S")),
        extent=[mwfre[0],mwfre[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',
        xlabel='PFre(dBm)',ylabel='Frequency(GHz)',cbarlabel=r"Amp(dBm)",
        color='jet',interpolations='None') 
plt.savefig(filename1+'_S21sweep_Amp'+'.png',format='png',dpi=600)
plt.close(fg1)
fg2=imshow(sadataplot_gain,figname='Amp_gain'+str(time.strftime("%H%M%S")),
        extent=[mwfre[0],mwfre[i],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',
        xlabel='PFre(dBm)',ylabel='Frequency(GHz)',cbarlabel=r"Gain(dB)",
        color='jet',interpolations='None') 
plt.savefig(filename1+'_Gain'+'.png',format='png',dpi=600)
plt.close(fg2)
mw.Close_Output()
np.savez(filename1,volt=volt,fre=fwr1,mwpwr=mwpower,sadataplot=sadataplot,sadataplot_gain=sadataplot_gain,data_mwoff=data_mwoff,vnapwr=vnapwr,about=about)

#%%改alpha测two_toneS43
vna=keysight_vna('vna2',trace=43)
vna.two_tone_vna()
fstart=6.472#6.483#6.479#6.47
fstop=6.472#6.483#6.  79#6.47
points=1001
ifband=30#每秒扫描次数
vnapwr=-5
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)

##%%初始化微波源
mw=MW('mw3')
# mw.start_Output()
mw.MW_setpower(-130)
mwpwr=-40#-45#-35
vary=3
fc0=9
mwfre_start=fc0-vary
mwfre_stop=fc0+vary
points=points
mw.VnaTrig(mwfre_start,mwfre_stop,points)
mw.MW_setpower(mwpwr)
# mw.start_Output()
fwr=np.linspace(mwfre_start,mwfre_stop,points)
fwr1=fwr[::-1]

which='S43'
name0='Global_'+str(Global)+which+'_'
filename1=filename+name0+time.strftime("%H%M%S")+'_'


Iz_list=np.linspace(0.2,-0.8,51)
alpha_list=[3.5]
gs1.setlevel_slowV2(Iz_list[0],0.005)#右z 
gs7.setlevel_slowV2(alpha_list[0],0.005)#右alpha#-1.07
# gs4.getlevel()
i=0;j=0
ap=[str(round(alpha_list[i],1)) for i in range(len(alpha_list))]
for alpha in alpha_list:
    gs7.setlevel_slow(alpha, 0.005)
    i=0
    name=name0+'alpha_'+ap[j]+'twotone_rightqubit_fsart{}_fstop{}_vnapwr{}_BW{}_mwpwr{}_'.format(fstart,fstop,vnapwr,ifband,mwpwr)+time.strftime("%H%M%S")
    filename1=filename+name
    S_data=[]
    print(time.strftime("%H:%M:%S"))
    mw.MW_setpower(-130)
    exec('bg=vna.get_data_'+which+'()')
    mw.MW_setpower(mwpwr)
    mw.start_Output()    
    for Iz in Iz_list:
        gs1.setlevel_slow(Iz, 0.005)
        # data2=vna.get_data_S43() 
        exec('data=vna.get_data_'+which+'()')
        S_data.append(data)
        # S43_data.append(data2)
        S_dataplot=np.rot90(np.array(S_data))
        S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
        Phase=np.angle(np.array(S_data))
        S_dataplot_phase=np.rot90(phase_handle(Phase,fwr))
        amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
        heatmap(np.linspace(Iz_list[0], Iz,i+1),fwr1,S_dataplot_amp,
                corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
                title='sweep_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        heatmap(np.linspace(Iz_list[0], Iz,i+1),fwr1,S_dataplot_phase,
                corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
                title='sweep_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
            # heatmap(np.linspace(Iz_list[0], Iz,i),f,S43_dataplot_amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep43_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
            # heatmap(np.linspace(Iz_list[0], Iz,i),f,S43_dataplot_phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep43_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        if i%(int((len(Iz_list)-1)/10))==0 and i!=0:
            print('scan_Fr: {:.0f}% completed.'.format((i+1)/len(Iz_list)*100))
            print(time.strftime("%H:%M:%S"))
        if i==0:
            t0=time.perf_counter()    
        elif i==1:
            t1=time.perf_counter()
            T=(t1-t0)*len(Iz_list)
            print('need {}min for one fig'.format(T/60))  
        i=i+1
    about=' BOX?; setup: in:port9,out:1;\n\
        Attenuation: Attenuation: Vna=40dB=30+10dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
            fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
                current from{}mA to{}mA'.format(fstart,fstop,points,ifband,vnapwr,Iz_list[0],Iz_list[-1])
    np.savez(filename1,fre=f,volt=Iz_list,Sdata=S_data,bg=bg,vnapwr=vnapwr,about=about,fwr=fwr,mwpwr=mwpwr)
    fg1=imshow(amp,figname='S43sweep_Amp'+str(time.strftime("%H%M%S")),
               extent=[Iz_list[0],Iz_list[-1],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',
               xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Amplitude(dB)",color='jet',interpolations='None') 
    plt.savefig(filename1+'_S43sweep_Amp'+'.png',format='png',dpi=600)
    plt.close(fg1)
    fg2=imshow(phase,figname='S43sweep_Phase'+str(time.strftime("%H%M%S")),
               extent=[Iz_list[0],Iz_list[-1],fwr1[-1],fwr1[0]],xround='%0.2f',yround='%0.3f',
               xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase(rad)",color='bwr',interpolations='None')
    plt.savefig(filename1+'_S43sweep_Phase'+'.png',format='png',dpi=600)
    plt.close(fg2)
    heatmap(Iz_list,fwr1,amp,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
            title='S43sweep_Amp',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    heatmap(Iz_list,fwr1,phase,corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',
            title='S43sweep_Phase',filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    print(j)
    j=j+1
    Iz_list=Iz_list[::-1]
    mw.MW_setpower(-130)
    mw.Close_Output()
# gs3.setlevel_slow(0,0.005)#右z
# gs1.setlevel_slow(0,0.01)#右alpha#-1.07
#%%AC stsrk S43
vna=keysight_vna('vna2',trace=43)
vna.two_tone_vna()
fstart=6.48#6.483#6.479#6.47
fstop=6.48#6.483#6.479#6.47
points=501
ifband=30#每秒扫描次数
vnapwr=-10
f=np.linspace(fstop,fstart,points)
f1=f[::-1]    
vna.set_power(vnapwr)
vna.set_startstopFre(fstart,fstop)
vna.set_points_band(points,ifband)
# mw=MW('mw3')
mw.MW_setpower(-130)
mwpwr=-40#-45#-35
vary=0.25
fc0=7.5
mwfre_start=fc0-vary
mwfre_stop=fc0+vary
points=points
mw.VnaTrig(mwfre_start,mwfre_stop,points)
mw.MW_setpower(mwpwr)
mw.start_Output()
fwr=np.linspace(mwfre_start,mwfre_stop,points)
fwr1=fwr[::-1]
# mw.Close_Output()
Global =0.24#0.015#-0.016 #-0.0308#0.382;
# gs5.setlevel_slowV2(Global,0.0001)
which='S43';att='20dB'
gs7.setlevel_slow(4.2,0.005)
alpha=gs7.getlevel()*1000
volt=-1
gs1.setlevel_slow(volt,0.005)#-0.9,0.05
name0='Global_'+str(Global)+which+'Acstark_leftqubit_alpha{}_att{}_f{}_mwpwr{}_volt{}_'.format(alpha,att,fstart,mwpwr,volt)
filename1=filename+name0+time.strftime("%H%M%S")+'_'
vnapwr=np.linspace(-30,6,73)
S_data=[]
i=0
print(time.strftime("%H:%M:%S"))
for i in tqdm(range(len(vnapwr))):
    vna.set_power(vnapwr[i])
    exec('data=vna.get_data_'+which+'()')
    S_data.append(data)
    S_dataplot=np.rot90(np.array(S_data))
    S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
    Phase=np.angle(np.array(S_data))
    S_dataplot_phase=np.rot90(phase_handle(Phase,fwr))
    amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
    heatmap(np.linspace(vnapwr[0],vnapwr[i],i),fwr1,S_dataplot_amp,
        corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='AcSweep_Amp',
        filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    heatmap(np.linspace(vnapwr[0],vnapwr[i],i),fwr1,S_dataplot_phase,
            corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='AcSweep_Phase',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    # if i!=0:
    #     fg1.clf()
    #     fg2.clf()
    # fg1=imshow(amp,figname='Two_Tone_Amp',extent=[volt[0],volt[i],fwr[0],fwr[-1]],
    #            xround='%0.2f',yround='%0.3f',xlabel=r' volt (mA)',ylabel=r'Frequency (GHz)',
    #            cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    # fg2=imshow(phase,figname='Two_Tone_Phase',extent=[volt[0],volt[i],fwr[0],fwr[-1]],
    #            xround='%0.2f',yround='%0.3f',xlabel=r'volt (mA)',ylabel=r'Frequency (GHz)',
    #            cbarlabel=r"Amplitude(dB)",color='bwr',interpolations='None') 
    if i==0:
        t0=time.perf_counter()    
    elif i==1:
        t1=time.perf_counter()
        T=(t1-t0)*len(vnapwr)
        print('need {}min for one fig'.format(T/60))  
mw.MW_setpower(-130)
mw.Close_Output()
# gs2.setlevel_slow(0,0.005)
# gs6.setlevel_slow(0,0.005)  
#%%mwpwr,vnapwr
def two_tone_mwpwr(mwpwr,vnapwr):
    vna.set_power(vnapwr)
    mw.MW_setpower(mwpwr)
    mw.start_Output()
    filename='E:\\Data\\RXH\\EP\\20211128\\'
    name0='twotone_leftqubit_fsart{}_fstop{}_vnapwr{}_BW{}_mwpwr{}_'.format(fstart,fstop,vnapwr,ifband,mwpwr)
    filename1=filename+name0+time.strftime("%H%M%S")+'_'
    bg=vna.get_data_S21() 
    volt=np.linspace(2.15,2.65,26)
    gs3.setlevel_slow(2.15,0.005)#-0.9,0.05
    gs3.Start_OutPut() #Couple Qubit
    S_data=[]
    i=0
    print(time.strftime("%H:%M:%S"))
    for i in range(len(volt)):
        gs3.setlevel_slow(volt[i],0.005)
        # for j in range(len(mwf)):
        data=vna.get_data_S21() 
        S_data.append(data)
        S_dataplot=np.rot90(np.array(S_data))
        S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
        Phase=np.angle(np.array(S_data))
        S_dataplot_phase=np.rot90(phase_handle(Phase,fwr))
        amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
        heatmap(np.linspace(volt[0],volt[i],i),fwr1,S_dataplot_amp,
            corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Sweep_Amp',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        heatmap(np.linspace(volt[0],volt[i],i),fwr1,S_dataplot_phase,
                corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep_Phase',
                filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        if i%(int(len(volt)/10))==0:
            print('scan_Fr: {:.0f}% completed.'.format((i+1)/len(volt)*100))  
            print(time.strftime("%H:%M:%S"))
    
            # 
        if i==0:
            t0=time.perf_counter()    
        elif i==1:
            t1=time.perf_counter()
            T=(t1-t0)*len(volt)
            print('need {}min for one fig'.format(T/60))  
    mw.MW_setpower(-130)
    mw.Close_Output()
    Sdata=np.array(S_data)
    about=' BOX?; setup: in:port9,out:1;\n\
        Attenuation: Attenuation: Vna=40dB=30+10dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
            fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
                current from{}mA to{}mA'.format(fstart,fstop,points,ifband,vnapwr,volt[0],volt[-1])
    #filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_local_Couplequbit'
    np.savez(filename1,fre=f,volt=volt,Sdata=Sdata,bg=bg,vnapwr=vnapwr,about=about,fwr=fwr,mwpwr=mwpwr)
    #savetxt(filename2,about)
    fg1=imshow(amp,figname='Two_Tone_Amp',
               extent=[volt[0],volt[i],fwr[0],fwr[-1]],xround='%0.2f',yround='%0.3f',
               xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',
               cbarlabel=r"Amplitude(dB)",color='jet',interpolations='None') 
    plt.savefig(filename1+'S21sweep_Amp'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg1)
    fg2=imshow(phase,figname='Two_Tone_Phase',
               extent=[volt[0],volt[i],fwr[0],fwr[-1]],xround='%0.2f',yround='%0.3f',
               xlabel=r'Current(mA)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase(rad)",
               color='bwr',interpolations='None')
    plt.savefig(filename1+'S21sweep_Phase'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg2)
    heatmap(np.linspace(volt[0],volt[i],i),fwr1,amp,
            corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Two_Tone_Amp',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    heatmap(np.linspace(volt[0],volt[i],i),fwr1,phase,
            corlor_scale='bwr',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Two_Tone_Phase',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    pass
mwpwrlist=[-60]
vnapwrlist=[0]
for vnapwr in vnapwrlist:
    for mwpwr in mwpwrlist:
        two_tone_mwpwr(mwpwr,vnapwr)
#%%mwpwrlist
def two_tone_mwpwrlist(mwpwrlist,vnapwr):
    vna.set_power(vnapwr)
    # filename='E:\\Data\\RXH\\EP\\20211129\\'
    name0='twotone_middlequbit_sweetpoint_fsart{}_fstop{}_vnapwr{}_BW{}_'.format(fstart,fstop,vnapwr,ifband)
    filename1=filename+name0+time.strftime("%H%M%S")+'_'
    mw.start_Output()
    bg=vna.get_data_S21() 
    volt=1.74
    # gs2.Start_OutPut() #left Qubit
    gs6.setlevel_slow(volt,0.005)#-0.9,0.05
    S_data=[]
    i=0
    for mwpwr in mwpwrlist:
        i=i+1
        mw.MW_setpower(mwpwr)
        time.sleep(1)
        data=vna.get_data_S21() 
        S_data.append(data)
        S_dataplot=np.rot90(np.array(S_data))
        S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
        Phase=np.angle(np.array(S_data))
        S_dataplot_phase=np.rot90(phase_handle(Phase,fwr))
        amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
        heatmap(np.linspace(mwpwrlist[0],mwpwr,i),fwr1,S_dataplot_amp,
            corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='Sweep_Amp',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        heatmap(np.linspace(mwpwrlist[0],mwpwr,i),fwr1,S_dataplot_phase,
                corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='sweep_Phase',
                filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
        if i%(int(len(mwpwrlist)/10))==0:
            print('scan_Fr: {:.0f}% completed.'.format((i+1)/len(mwpwrlist)*100))  
            print(time.strftime("%H:%M:%S"))
    
            # 
        if i==1:
            t0=time.perf_counter()    
        elif i==2:
            t1=time.perf_counter()
            T=(t1-t0)*len(mwpwrlist)
            print('need {}min for one fig'.format(T/60))  

    mw.MW_setpower(-130)
    mw.Close_Output()
    Sdata=np.array(S_data)
    about=' BOX?; setup: in:port9,out:1;\n\
        Attenuation: Attenuation: Vna=40dB=30+10dB, inline=60dB=42+8+10dB, outline=76dB(RT)+36dB(4K) \n\
            fstart={}Ghz,fstop={}Ghz,points={},ifband={}Hz,vnapwr={}dB \n\
                mwpwr from{}dBm to{}dBm'.format(fstart,fstop,points,ifband,vnapwr,mwpwrlist[0],mwpwrlist[-1])
    #filename2=filename+str(time.strftime("%H%M%S"))+'PTv4_3bit_local_Couplequbit'
    np.savez(filename1,fre=f,volt=volt,bg=bg,Sdata=Sdata,vnapwr=vnapwr,about=about,fwr=fwr,mwpwrlist=mwpwrlist)
    #savetxt(filename2,about)
    fg1=imshow(amp,figname='Two_Tone_Amp'+str(time.strftime("%H%M%S")),
               extent=[mwpwrlist[0],mwpwrlist[-1],fwr[0],fwr[-1]],xround='%0.2f',yround='%0.3f',
               xlabel=r'Power(dBm)',ylabel=r'Frequency (GHz)',
               cbarlabel=r"Amplitude(dB)",color='jet',interpolations='None') 
    plt.savefig(filename1+'S21sweep_Amp'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg1)
    fg2=imshow(phase,figname='Two_Tone_Phase'+str(time.strftime("%H%M%S")),
               extent=[mwpwrlist[0],mwpwrlist[-1],fwr[0],fwr[-1]],xround='%0.2f',yround='%0.3f',
               xlabel=r'Power(dBm)',ylabel=r'Frequency (GHz)',cbarlabel=r"Phase(rad)",
               color='bwr',interpolations='None')
    plt.savefig(filename1+'S21sweep_Phase'+str(time.strftime("%H%M%S"))+'.png',format='png',dpi=600)
    plt.close(fg2)
    heatmap(np.linspace(mwpwrlist[0],mwpwr,i),fwr1,amp,
            corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='Two_Tone_Amp',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    heatmap(np.linspace(mwpwrlist[0],mwpwr,i),fwr1,phase,
            corlor_scale='bwr',xlabel='Power[dBm]',ylabel='Frequency[Ghz]',title='Two_Tone_Phase',
            filename=filename1,errobar=[],zmin=0,zmax=1,zauto=True)
    pass
mw=MW('mw3')
mw.MW_setpower(-130)
mwpwr=-130#-45#-35
vary=0.2
fc0=10.922
mwfre_start=fc0-vary
mwfre_stop=fc0+vary
points=points
mw.VnaTrig(mwfre_start,mwfre_stop,points)
mw.MW_setpower(mwpwr)
mw.start_Output()
fwr=np.linspace(mwfre_start,mwfre_stop,points)
fwr1=fwr[::-1]
mwpwrlist=np.linspace(-30,-70,41)
vnapwrlist=[0]
for vnapwr in vnapwrlist:
    two_tone_mwpwrlist(mwpwrlist,vnapwr)