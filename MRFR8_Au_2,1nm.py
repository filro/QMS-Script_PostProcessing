# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 10:30:40 2022

@author: filro
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 09:36:38 2022

@author: filro
"""
# =========================================================================== #
#                            Pre-processing section                           #
# =========================================================================== #

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from collections import OrderedDict

from scipy import interpolate 

linestyles_dict = OrderedDict(
    [('solid',               (0, ())),
    
     ('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely dotted',      (0, (1, 1))),

     ('loosely dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely dashed',      (0, (5, 1))),

     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])

# Data import
df = pd.read_csv("C:/Users/filro/OneDrive - Danmarks Tekniske Universitet/Skrivebord/Data/QMS/MRFR8_Au_2nm_ClusterSource_CO_Oxidation_1bar_TCcontrolled_2to1_O2toCO.txt",
                 sep="\t", skiprows=20)

print(df)

# =========================================================================== #
#                              Processing section                             #
# =========================================================================== #

#_____find plateau______#
Index_tmin = min(df.index[df['M40-x'] > 0.2].tolist())
Index_tmax = max(df.index[df['M40-x'] < 1.2].tolist())
Baseline_Ar = np.mean(df["M40-y"][Index_tmin:Index_tmax])

FM_CO2 = 11.23 #A*s/mol
FM_Ar = 9.22 #A*s/mol
FM_H2 = 9.50 #A*s/mol
FM_CO = 10.45 #A*s/mol
FM_O2 = 7.27 #A*s/mol

r = 0.5 #cm
Area_microreactor = math.pi*r**2 #cm^2
MW_Au = 196.97 #ng/nmol
LoadingAu_Au_W = 134 #ng/cm^2
LoadingAu_Au_mol = LoadingAu_Au_W/MW_Au #nmol/cm^2
Actual_LoadingAu_Au = Area_microreactor*LoadingAu_Au_mol #nmol

AreaSingleParticle = math.pi*(2.1/2)**2 #nm^2
TotalAreaCovered_Microreactor = 5/100*Area_microreactor*1e14 #5 % of the area in nm2
N_of_particles = TotalAreaCovered_Microreactor/AreaSingleParticle
AreaParticle = 4*math.pi*(3.7/2)**2 #nm^2
ActiveArea_Au = N_of_particles*AreaParticle/1e14 #cm2 = 0.48762362162861

#_____________________________Ar baselined data_______________________________#
df["M28-y-new"] = df["M28-y"]/df["M40-y"]*Baseline_Ar*1e9/FM_CO*60
df["M32-y-new"] = df["M32-y"]/df["M40-y"]*Baseline_Ar*1e9/FM_O2*60
df["M40-y-new"] = df["M40-y"]/df["M40-y"]*Baseline_Ar*1e9/FM_Ar*60
df["M44-y-new"] = df["M44-y"]/df["M40-y"]*Baseline_Ar*1e9/FM_CO2*60

df["M44-y-TOF"] = df["M44-y-new"]/ActiveArea_Au

#____________INTERPOLATION_________________#
# f44 = interpolate.interp1d(df["M44-x"],df["M44-y"],fill_value="extrapolate")
# df["M44-y-new"] = f44(df["M40-x"])/df["M40-y"]*Baseline_Ar
# df["M44-y-new"] = f44(df["M40-x"])
#___________END INTERPOLATION______________#
# df["M45-y-new"] = df["M45-y"]/df["M40-y"]*Baseline_Ar

#_______________________Averaging and trend - CO2 SIGNAL______________________#

# index of signals at different temperatures

T92min = min(df.index[df['M28-x'] > 2.3].tolist())
T92max = min(df.index[df['M28-x'] > 2.45].tolist())

T107min = min(df.index[df['M28-x'] > 3.3].tolist())
T107max = min(df.index[df['M28-x'] > 3.45].tolist())

T122min = min(df.index[df['M28-x'] > 4.3].tolist())
T122max = min(df.index[df['M28-x'] > 4.45].tolist())

T137min = min(df.index[df['M28-x'] > 5.3].tolist())
T137max = min(df.index[df['M28-x'] > 5.45].tolist())

T152min = min(df.index[df['M28-x'] > 6.3].tolist())
T152max = min(df.index[df['M28-x'] > 6.45].tolist())

T167min = min(df.index[df['M28-x'] > 7.3].tolist())
T167max = min(df.index[df['M28-x'] > 7.45].tolist())

T182min = min(df.index[df['M28-x'] > 8.3].tolist())
T182max = min(df.index[df['M28-x'] > 8.45].tolist())

T197min = min(df.index[df['M28-x'] > 9.3].tolist())
T197max = min(df.index[df['M28-x'] > 9.45].tolist())

T213min = min(df.index[df['M28-x'] > 10.3].tolist())
T213max = min(df.index[df['M28-x'] > 10.45].tolist())

T228min = min(df.index[df['M28-x'] > 11.3].tolist())
T228max = min(df.index[df['M28-x'] > 11.45].tolist())

T243min = min(df.index[df['M28-x'] > 12.3].tolist())
T243max = min(df.index[df['M28-x'] > 12.45].tolist())

Tambmin = min(df.index[df['M28-x'] > 0.2].tolist())
Tambmax = min(df.index[df['M28-x'] > 1.2].tolist())

# Tbackgroundmin = min(df.index[df['M28-x'] > 0.2].tolist())
# Tbackgroundmax = min(df.index[df['M28-x'] > 1.2].tolist())

# Averaging signal for CO2 --> CO2 conversion 

# Signal_T_background = df["M28-y-new"][Tbackgroundmin:Tbackgroundmax].mean()
Signal_T_M44 = pd.DataFrame()
Signal_T_M44['temperature'] = []
Signal_T_M44['mean'] = []

Signal_T_M44.loc[0] = {'temperature':30,'mean':df["M44-y-new"][Tambmin:Tambmax].mean()}
Signal_T_M44.loc[1] = {'temperature':92,'mean':df["M44-y-new"][T92min:T92max].mean()}
Signal_T_M44.loc[2] = {'temperature':107,'mean':df["M44-y-new"][T107min:T107max].mean()}
Signal_T_M44.loc[3] = {'temperature':122,'mean':df["M44-y-new"][T122min:T122max].mean()}
Signal_T_M44.loc[4] = {'temperature':137,'mean':df["M44-y-new"][T137min:T137max].mean()}
Signal_T_M44.loc[5] = {'temperature':152,'mean':df["M44-y-new"][T152min:T152max].mean()}
Signal_T_M44.loc[6] = {'temperature':167,'mean':df["M44-y-new"][T167min:T167max].mean()}
Signal_T_M44.loc[7] = {'temperature':182,'mean':df["M44-y-new"][T182min:T182max].mean()}
Signal_T_M44.loc[8] = {'temperature':197,'mean':df["M44-y-new"][T197min:T197max].mean()}
Signal_T_M44.loc[9] = {'temperature':213,'mean':df["M44-y-new"][T213min:T213max].mean()}
Signal_T_M44.loc[10] = {'temperature':228,'mean':df["M44-y-new"][T228min:T228max].mean()}
Signal_T_M44.loc[11] = {'temperature':243,'mean':df["M44-y-new"][T243min:T243max].mean()}

Signal_T_amb = df["M28-y-new"][Tambmin:Tambmax].mean() - 0.1*Signal_T_M44["mean"][0]
Signal_T_92 = df["M28-y-new"][T92min:T92max].mean() -  0.1*Signal_T_M44["mean"][1]
Signal_T_107 = df["M28-y-new"][T107min:T107max].mean() -  0.1*Signal_T_M44["mean"][2]
Signal_T_122 = df["M28-y-new"][T122min:T122max].mean() -  0.1*Signal_T_M44["mean"][3]
Signal_T_137 = df["M28-y-new"][T137min:T137max].mean() -  0.1*Signal_T_M44["mean"][4]
Signal_T_152 = df["M28-y-new"][T152min:T152max].mean() -  0.1*Signal_T_M44["mean"][5]
Signal_T_167 = df["M28-y-new"][T167min:T167max].mean() -  0.1*Signal_T_M44["mean"][6]
Signal_T_182 = df["M28-y-new"][T182min:T182max].mean() -  0.1*Signal_T_M44["mean"][7]
Signal_T_197 = df["M28-y-new"][T197min:T197max].mean() -  0.1*Signal_T_M44["mean"][8]
Signal_T_213 = df["M28-y-new"][T213min:T213max].mean() -  0.1*Signal_T_M44["mean"][9]
Signal_T_228 = df["M28-y-new"][T228min:T228max].mean() -  0.1*Signal_T_M44["mean"][10]
Signal_T_243 = df["M28-y-new"][T243min:T243max].mean() -  0.1*Signal_T_M44["mean"][11]


Conversion_T = pd.DataFrame()
Conversion_T['temperature'] = []
Conversion_T['conversion'] = []

Conversion_T.loc[0] = {'temperature':30,'conversion':((Signal_T_amb-Signal_T_amb)/Signal_T_amb) * 100}
Conversion_T.loc[1] = {'temperature':92,'conversion':((Signal_T_amb-Signal_T_92)/Signal_T_amb) * 100}
Conversion_T.loc[2] = {'temperature':107,'conversion':((Signal_T_amb-Signal_T_107)/Signal_T_amb) * 100}
Conversion_T.loc[3] = {'temperature':122,'conversion':((Signal_T_amb-Signal_T_122)/Signal_T_amb) * 100}
Conversion_T.loc[4] = {'temperature':137,'conversion':((Signal_T_amb-Signal_T_137)/Signal_T_amb) * 100}
Conversion_T.loc[5] = {'temperature':152,'conversion':((Signal_T_amb-Signal_T_152)/Signal_T_amb) * 100}
Conversion_T.loc[6] = {'temperature':167,'conversion':((Signal_T_amb-Signal_T_167)/Signal_T_amb) * 100}
Conversion_T.loc[7] = {'temperature':182,'conversion':((Signal_T_amb-Signal_T_182)/Signal_T_amb) * 100}
Conversion_T.loc[8] = {'temperature':197,'conversion':((Signal_T_amb-Signal_T_197)/Signal_T_amb) * 100}
Conversion_T.loc[9] = {'temperature':213,'conversion':((Signal_T_amb-Signal_T_213)/Signal_T_amb) * 100}
Conversion_T.loc[10] = {'temperature':228,'conversion':((Signal_T_amb-Signal_T_228)/Signal_T_amb) * 100}
Conversion_T.loc[11] = {'temperature':243,'conversion':((Signal_T_amb-Signal_T_243)/Signal_T_amb) * 100}

# Averaging signal for M31 (MeOH) ---> MeOH production vs T


Signal_T_M28 = pd.DataFrame()
Signal_T_M28['temperature'] = []
Signal_T_M28['mean'] = []

Signal_T_M28.loc[0] = {'temperature':30,'mean':Signal_T_amb}
Signal_T_M28.loc[1] = {'temperature':92,'mean':Signal_T_92}
Signal_T_M28.loc[2] = {'temperature':107,'mean':Signal_T_107}
Signal_T_M28.loc[3] = {'temperature':122,'mean':Signal_T_122}
Signal_T_M28.loc[4] = {'temperature':137,'mean':Signal_T_137}
Signal_T_M28.loc[5] = {'temperature':152,'mean':Signal_T_152}
Signal_T_M28.loc[6] = {'temperature':167,'mean':Signal_T_167}
Signal_T_M28.loc[7] = {'temperature':182,'mean':Signal_T_182}
Signal_T_M28.loc[8] = {'temperature':197,'mean':Signal_T_197}
Signal_T_M28.loc[9] = {'temperature':213,'mean':Signal_T_213}
Signal_T_M28.loc[10] = {'temperature':228,'mean':Signal_T_228}
Signal_T_M28.loc[11] = {'temperature':243,'mean':Signal_T_243}

#_______________________Averaging and trend - CO SIGNAL______________________#
#_________________________________Second Ramp________________________________#
# index of signals at different temperatures

T30min_2ndramp = min(df.index[df['M28-x'] > 28.5].tolist())
T30max_2ndramp = min(df.index[df['M28-x'] > 29].tolist())

T47min_2ndramp = min(df.index[df['M28-x'] > 30].tolist())
T47max_2ndramp = min(df.index[df['M28-x'] > 30.1].tolist())

T62min_2ndramp = min(df.index[df['M28-x'] > 31].tolist())
T62max_2ndramp = min(df.index[df['M28-x'] > 31.1].tolist())

T77min_2ndramp = min(df.index[df['M28-x'] > 32].tolist())
T77max_2ndramp = min(df.index[df['M28-x'] > 32.15].tolist())

T92min_2ndramp = min(df.index[df['M28-x'] > 33].tolist())
T92max_2ndramp = min(df.index[df['M28-x'] > 33.15].tolist())

T107min_2ndramp = min(df.index[df['M28-x'] > 34].tolist())
T107max_2ndramp = min(df.index[df['M28-x'] > 34.15].tolist())

T122min_2ndramp = min(df.index[df['M28-x'] > 35].tolist())
T122max_2ndramp = min(df.index[df['M28-x'] > 35.15].tolist())

T137min_2ndramp = min(df.index[df['M28-x'] > 36].tolist())
T137max_2ndramp = min(df.index[df['M28-x'] > 36.2].tolist())

T152min_2ndramp = min(df.index[df['M28-x'] > 37].tolist())
T152max_2ndramp = min(df.index[df['M28-x'] > 37.2].tolist())

T167min_2ndramp = min(df.index[df['M28-x'] > 38].tolist())
T167max_2ndramp = min(df.index[df['M28-x'] > 38.2].tolist())

T182min_2ndramp = min(df.index[df['M28-x'] > 39].tolist())
T182max_2ndramp = min(df.index[df['M28-x'] > 39.20].tolist())

T197min_2ndramp = min(df.index[df['M28-x'] > 40].tolist())
T197max_2ndramp = min(df.index[df['M28-x'] > 40.2].tolist())

T213min_2ndramp = min(df.index[df['M28-x'] > 41].tolist())
T213max_2ndramp = min(df.index[df['M28-x'] > 41.2].tolist())

T228min_2ndramp = min(df.index[df['M28-x'] > 42].tolist())
T228max_2ndramp = min(df.index[df['M28-x'] > 42.3].tolist())

T243min_2ndramp = min(df.index[df['M28-x'] > 43].tolist())
T243max_2ndramp = min(df.index[df['M28-x'] > 43.3].tolist())

# =========================================================================== #
#                               Figure section                                #
# =========================================================================== #

###############################################################################
#              FIGURE 1 - SIGNAL TREND WITH TEMPERATURE RAMPS                 #
###############################################################################

fig1, ax1 = plt.subplots()

#___________________________Left part of the graph____________________________#
#_____Uncomment the lines related to the molecules you are interested in______#

#______________________________RAW data plot__________________________________#

# ax1.plot(df["M28-x"], df["M28-y"], label ='CO, M28', linewidth = 0.7, c ='#C17E25')
# ax1.plot(df["M32-x"], df["M32-y"], label ='O$_2$, M32', linewidth = 0.7, c ='#527D8E')
# ax1.plot(df["M40-x"], df["M40-y"], label ='Ar, M40', linewidth = 0.7, c ='#9653B0')
# ax1.plot(df["M44-x"], df["M44-y"], label ='CO$_2$, M44', linewidth = 0.7, c ='#8EAE7F')

#___________________________BASELINED data plot_______________________________#

ax1.plot(df["M28-x"], df["M28-y-new"], label ='M28, CO', linewidth = 1.5, c ='#C17E25')
ax1.plot(df["M32-x"], df["M32-y-new"], label ='M32, O$_2$', linewidth = 1.5, c ='#527D8E')
# ax1.plot(df["M40-x"], df["M40-y-new"], label ='M40, Ar', linewidth = 1.5, c ='#9653B0')
ax1.plot(df["M44-x"], df["M44-y-new"], label ='M44, CO$_2$', linewidth = 1.5, c ='#8EAE7F')

#____________________Right part of the graph (Temperature)____________________#
#_________________comment/uncomment for RTD or TC controlled__________________#

ax2 = ax1.twinx() 

#_____________________________if RTD controlled_______________________________#

# ax2.plot(df["RTD temperature-x"], df["RTD temperature-y"], label ='T [°C]',
                # color = 'tab:red', linewidth = 0.5, alpha=1,
                # linestyle=linestyles_dict['loosely dotted'])

#______________________________if TC controlled_______________________________#

df["TC temperature-y"] = 1.5054704595186*df["TC temperature-y"]-13.306783369803
ax2.plot(df["TC temperature-x"], df["TC temperature-y"], label ='T [°C]', 
                color = 'tab:red', linewidth = 1, alpha=1, 
                linestyle=linestyles_dict['dotted'])

# Plotting settings

# plt.title('MRFR8 - 2,1 nm Au (cluster source), CO oxidation at 1 bar, 11-12-2022', 
#             fontsize=10, ha="center", fontweight="bold")

ax1.set_xlabel('Time [$\it{h}$]', fontsize = 10) 
ax1.set_ylabel('Flow $\it{[nmol/min]}$', fontsize = 10)
# ax1.set_ylabel('SEM current $\it{[A]}$', fontsize = 10)
ax1.tick_params(axis ='y', labelcolor = 'black') 
# ax1.legend(loc="best", fontsize = 'x-small', facecolor="white",
#            edgecolor="white", framealpha= 1)
# ax1.legend(fontsize = 'small', loc='center left', bbox_to_anchor=(1.15, 0.5),
#           fancybox=True, shadow=True, ncol=1)
# ax1.set_yscale('log')
ax1.set_ylim(2e-1, 15)
ax1.set_xlim(28,70)
# ax1.grid()
plt.yticks(fontsize = 10)

ax2.set_ylabel('Temperature $\it{(°C)}$', color = 'tab:red', fontsize=10) 
ax2.tick_params(axis ='y', labelcolor = 'tab:red')
# ax2.legend(loc="lower right", fontsize = 'small', facecolor="white", 
#             edgecolor="white", framealpha= 1)
# ax2.set_yscale('log')
# ax2.set_ylim(0, 210)
# ax2.grid()

###############################################################################
#         FIGURE 2 - CO2 CONVERSION VS TEMPERATURE AND M31 SIGNAL             #
###############################################################################

fig2, ax1 = plt.subplots()

ax1.plot(Conversion_T['temperature'], Conversion_T['conversion'], label = 'CO conversion', color = 'k')  
# ax1.plot(DataSub['T'], Subtracted, label = 'Subtracted Data')
# ax1.plot(DataSub['T'], ToSubtract, label = 'Function to subtract')

# ax1.set_ylim(0, 0.175e-11)

ax2 = ax1.twinx() 
ax2.plot(Signal_T_M28["temperature"], Signal_T_M28["mean"], label ='CO Signal (M28)', 
                color = 'tab:red', linewidth = 1.5, alpha=1, 
                linestyle=linestyles_dict['densely dotted'])
# ax2.legend(fontsize = 'x-small',  bbox_to_anchor=(0.263, 0.95),
#           fancybox=False, shadow=False, frameon=False, ncol=1)

# plt.title('CO Conversion', fontsize=10, ha="center", fontweight="bold")
# ax1.legend(fontsize = 'x-small', bbox_to_anchor=(0.28, 0.15),
#           fancybox=False, shadow=False, frameon=False, ncol=1)

ax1.set_xlabel('Temperature $\it{[°C]}$', fontsize = 10) 
ax1.set_ylabel('Conversion, $\chi$ [%]', fontsize = 10)
ax2.set_ylabel('Flow $\it{[nmol^{CO}/min]}$', color = 'tab:red', fontsize=10) 
ax2.tick_params(axis ='y', labelcolor = 'tab:red')
# ax2.set_ylim(0, 510)
# ax2.grid()
plt.show()

# =========================================================================== #
#                                 File saving                                 #
# =========================================================================== #

fig1.savefig('C:/Users/filro/OneDrive - Danmarks Tekniske Universitet/Skrivebord/Data/QMS/QMS_Images/MRFR8_Au_2nm_ClusterSource_CO_Oxidation_1bar_Stability.png',
            format='png', dpi=1200, bbox_inches="tight")

# fig2.savefig('C:/Users/filro/OneDrive - Danmarks Tekniske Universitet/Skrivebord/Data/QMS/QMS_Images/MRFR8_Au_2nm_ClusterSource_CO_Oxidation_1bar_TCcontrolled_2to1_O2toCO_Conversion.png',
#             format='png', dpi=1200, bbox_inches="tight")

# df.to_csv('C:/Users/filro/OneDrive - Danmarks Tekniske Universitet/Skrivebord/Data/QMS/QMS_ProcessedData/ProcessedData_MRFR4_NiGa3_2nmSputtering_CO2_reduction_Methanol_250C_1bar_Activation250.txt', header=True)
