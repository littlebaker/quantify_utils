import numpy as np
from def_plotly import heatmap
from hunanu.Functions.function0 import savetxt,plot_removebg,todB,imshow,phase_handle,Lorentz_fitting_VNA,half_value_width


data = np.load("/home/machine1/E/Data/LSY_new/data/CELv1/twotone_energyspectrum(alpha=0mA)/two_tone_qubit_freq_2GHz_to_12GHz_cavity_freq_6.449GHz1022/two_tone_qubit_freq_2GHz_to_12GHz_cavity_freq_6.449GHz154947_.npz")

data.files

freq = data["fre"]
volt = data["volt"]
S_data = data["Sdata"]


f1 = np.array([6.449]*len(freq))
fwr1 = freq[::-1]

S_dataplot=np.rot90(np.array(S_data))
S_dataplot_amp=20*np.log10(np.abs(S_dataplot))
Phase=np.angle(np.array(S_data))
S_dataplot_phase=np.rot90(phase_handle(Phase, f1, 50.11))
amp,phase=plot_removebg(S_dataplot_amp,S_dataplot_phase)
heatmap(np.linspace(volt[0],volt[-1],len(volt)),fwr1,S_dataplot_amp,
    corlor_scale='hot',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='Sweep_Amp',
    errobar=[],zmin=-78,zmax=-65,zauto=False)
heatmap(np.linspace(volt[0],volt[-1],len(volt)),fwr1,S_dataplot_phase,
        corlor_scale='hot',xlabel='Current[mA]',ylabel='Frequency[Ghz]',title='sweep_Phase',
        errobar=[],zmin=-1,zmax=5,zauto=False)
