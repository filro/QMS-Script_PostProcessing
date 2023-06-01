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
from collections import OrderedDict
import math
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
df = pd.read_csv("C:/Users/filro/OneDrive - Danmarks Tekniske Universitet/Skrivebord/Data/QMS/MRFR14_3,5nmAuTi_onTiO2_CO_Oxidation_1bar_ReproducibilityTest_Stability.txt",
                 sep="\t", skiprows=20)

print(df)

# =========================================================================== #
#                              Processing section                             #
# =========================================================================== #

#_____find plateau______#
Index_tmin = min(df.index[df['M40-x'] > 37].tolist())
Index_tmax = max(df.index[df['M40-x'] < 37.2].tolist())
Baseline_Ar = np.mean(df["M40-y"][Index_tmin:Index_tmax])

FM_CO2 = 13.50 #A*s/mol
FM_Ar = 12.79 #A*s/mol
FM_H2 = 8.87 #A*s/mol
FM_CO = 15.64 #A*s/mol
FM_O2 = 12.08 #A*s/mol

r = 0.5 #cm
Area_microreactor = math.pi*r**2 #cm^2
MW_Au = 196.97 #ng/nmol
LoadingAu_AuTi_W = 48 #ng/cm^2
LoadingAu_AuTi_mol = LoadingAu_AuTi_W/MW_Au #nmol/cm^2

Actual_LoadingAu_AuTi = Area_microreactor*LoadingAu_AuTi_mol #nmol

AreaSingleParticle = math.pi*(3.7/2)**2 #nm^2
TotalAreaCovered_Microreactor = 5.7/100*Area_microreactor*1e14 #5 % of the area in nm2
N_of_particles = TotalAreaCovered_Microreactor/AreaSingleParticle
AreaParticle = 4*math.pi*(3.7/2)**2 #nm^2
Ratio_Au_Ti = (582-284)/949 #atomic ratio in the Ti overlayer
ActiveArea_Au = N_of_particles*AreaParticle*Ratio_Au_Ti/1e14 # cm2 = 0.05623086703253555

#_____________________Ar baselined data and FM corrected_______________________________#
# df = df_old/df40*BaselineArgon*mol_to_nmol/SensitivityFactor/loading*second_to_min

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

T40min = min(df.index[df['M28-x'] > 1.2].tolist())
T40max = min(df.index[df['M28-x'] > 1.4].tolist())

T60min = min(df.index[df['M28-x'] > 2].tolist())
T60max = min(df.index[df['M28-x'] > 2.2].tolist())

T80min = min(df.index[df['M28-x'] > 27.9].tolist())
T80max = min(df.index[df['M28-x'] > 28].tolist())

T100min = min(df.index[df['M28-x'] > 27.5].tolist())
T100max = min(df.index[df['M28-x'] > 27.8].tolist())

T120min = min(df.index[df['M28-x'] > 26.6].tolist())
T120max = min(df.index[df['M28-x'] > 26.7].tolist())

T140min = min(df.index[df['M28-x'] > 25.5].tolist())
T140max = min(df.index[df['M28-x'] > 25.6].tolist())

T160min = min(df.index[df['M28-x'] > 24.4].tolist())
T160max = min(df.index[df['M28-x'] > 24.5].tolist())

T180min = min(df.index[df['M28-x'] > 23.4].tolist())
T180max = min(df.index[df['M28-x'] > 23.5].tolist())

T200min = min(df.index[df['M28-x'] > 22].tolist())
T200max = min(df.index[df['M28-x'] > 22.9].tolist())

T220min = min(df.index[df['M28-x'] > 150].tolist())
T220max = min(df.index[df['M28-x'] > 155].tolist())

T240min = min(df.index[df['M28-x'] > 187].tolist())
T240max = min(df.index[df['M28-x'] > 188].tolist())

T260min = min(df.index[df['M28-x'] > 202].tolist())
T260max = min(df.index[df['M28-x'] > 203].tolist())

T280min = min(df.index[df['M28-x'] > 18.1].tolist())
T280max = min(df.index[df['M28-x'] > 18.2].tolist())

T300min = min(df.index[df['M28-x'] > 17.1].tolist())
T300max = min(df.index[df['M28-x'] > 17.2].tolist())

T320min = min(df.index[df['M28-x'] > 15.3].tolist())
T320max = min(df.index[df['M28-x'] > 15.4].tolist())

Tambmin = min(df.index[df['M28-x'] > 38].tolist())
Tambmax = min(df.index[df['M28-x'] > 38.1].tolist())

# Tbackgroundmin = min(df.index[df['M28-x'] > 0.2].tolist())
# Tbackgroundmax = min(df.index[df['M28-x'] > 1.2].tolist())

# Averaging signal for CO2

# Signal_T_background = df["M28-y-new"][Tbackgroundmin:Tbackgroundmax].mean()
Signal_T_M44 = pd.DataFrame()
Signal_T_M44['temperature'] = []
Signal_T_M44['mean'] = []

Signal_T_M44.loc[0] = {'temperature':25,'mean':df["M44-y-TOF"][Tambmin:Tambmax].mean()}
Signal_T_M44.loc[1] = {'temperature':40,'mean':df["M44-y-TOF"][T40min:T40max].mean()}
Signal_T_M44.loc[2] = {'temperature':60,'mean':df["M44-y-TOF"][T60min:T60max].mean()}
Signal_T_M44.loc[3] = {'temperature':80,'mean':df["M44-y-TOF"][T80min:T80max].mean()}
Signal_T_M44.loc[4] = {'temperature':100,'mean':df["M44-y-TOF"][T100min:T100max].mean()}
Signal_T_M44.loc[5] = {'temperature':120,'mean':df["M44-y-TOF"][T120min:T120max].mean()}
Signal_T_M44.loc[6] = {'temperature':140,'mean':df["M44-y-TOF"][T140min:T140max].mean()}
Signal_T_M44.loc[7] = {'temperature':160,'mean':df["M44-y-TOF"][T160min:T160max].mean()}
Signal_T_M44.loc[8] = {'temperature':180,'mean':df["M44-y-TOF"][T180min:T180max].mean()}
Signal_T_M44.loc[9] = {'temperature':200,'mean':df["M44-y-TOF"][T200min:T200max].mean()}
Signal_T_M44.loc[10] = {'temperature':220,'mean':df["M44-y-TOF"][T220min:T220max].mean()}
Signal_T_M44.loc[11] = {'temperature':240,'mean':df["M44-y-TOF"][T240min:T240max].mean()}
Signal_T_M44.loc[12] = {'temperature':260,'mean':df["M44-y-TOF"][T260min:T260max].mean()}
Signal_T_M44.loc[13] = {'temperature':280,'mean':df["M44-y-TOF"][T280min:T280max].mean()}
Signal_T_M44.loc[14] = {'temperature':300,'mean':df["M44-y-TOF"][T300min:T300max].mean()}
Signal_T_M44.loc[15] = {'temperature':320,'mean':df["M44-y-TOF"][T320min:T320max].mean()}

Error_M44 = pd.DataFrame()
Error_M44['Error'] = []

Error_M44.loc[0] = np.std(df["M44-y-TOF"][Tambmin:Tambmax])
Error_M44.loc[1] = np.std(df["M44-y-TOF"][T40min:T40max])
Error_M44.loc[2] = np.std(df["M44-y-TOF"][T60min:T60max])
Error_M44.loc[3] = np.std(df["M44-y-TOF"][T80min:T80max])
Error_M44.loc[4] = np.std(df["M44-y-TOF"][T100min:T100max])
Error_M44.loc[5] = np.std(df["M44-y-TOF"][T120min:T120max])
Error_M44.loc[6] = np.std(df["M44-y-TOF"][T140min:T140max])
Error_M44.loc[7] = np.std(df["M44-y-TOF"][T160min:T160max])
Error_M44.loc[8] = np.std(df["M44-y-TOF"][T180min:T180max])
Error_M44.loc[9] = np.std(df["M44-y-TOF"][T200min:T200max])
Error_M44.loc[10] = np.std(df["M44-y-TOF"][T220min:T220max])
Error_M44.loc[11] = np.std(df["M44-y-TOF"][T240min:T240max])
Error_M44.loc[12] = np.std(df["M44-y-TOF"][T260min:T260max])
Error_M44.loc[13] = np.std(df["M44-y-TOF"][T280min:T280max])
Error_M44.loc[14] = np.std(df["M44-y-TOF"][T300min:T300max])
Error_M44.loc[15] = np.std(df["M44-y-TOF"][T320min:T320max])

## M28 averaging for conversion

Signal_T_amb = df["M28-y-new"][Tambmin:Tambmax].mean()
Signal_T_40 = df["M28-y-new"][T40min:T40max].mean()
Signal_T_60 = df["M28-y-new"][T60min:T60max].mean()
Signal_T_80 = df["M28-y-new"][T80min:T80max].mean()
Signal_T_100 = df["M28-y-new"][T100min:T100max].mean()
Signal_T_120 = df["M28-y-new"][T120min:T120max].mean()
Signal_T_140 = df["M28-y-new"][T140min:T140max].mean()
Signal_T_160 = df["M28-y-new"][T160min:T160max].mean()
Signal_T_180 = df["M28-y-new"][T180min:T180max].mean()
Signal_T_200 = df["M28-y-new"][T200min:T200max].mean()
Signal_T_220 = df["M28-y-new"][T220min:T220max].mean()
Signal_T_240 = df["M28-y-new"][T240min:T240max].mean() 
Signal_T_260 = df["M28-y-new"][T260min:T260max].mean()
Signal_T_280 = df["M28-y-new"][T280min:T280max].mean()
Signal_T_300 = df["M28-y-new"][T300min:T300max].mean()
Signal_T_320 = df["M28-y-new"][T320min:T320max].mean()


Conversion_T = pd.DataFrame()
Conversion_T['temperature'] = []
Conversion_T['conversion'] = []

Conversion_T.loc[0] = {'temperature':25,'conversion':((Signal_T_amb-Signal_T_amb)/Signal_T_amb) * 100}
Conversion_T.loc[1] = {'temperature':40,'conversion':((Signal_T_amb-Signal_T_40)/Signal_T_amb) * 100}
Conversion_T.loc[2] = {'temperature':60,'conversion':((Signal_T_amb-Signal_T_60)/Signal_T_amb) * 100}
Conversion_T.loc[3] = {'temperature':80,'conversion':((Signal_T_amb-Signal_T_80)/Signal_T_amb) * 100}
Conversion_T.loc[4] = {'temperature':100,'conversion':((Signal_T_amb-Signal_T_100)/Signal_T_amb) * 100}
Conversion_T.loc[5] = {'temperature':120,'conversion':((Signal_T_amb-Signal_T_120)/Signal_T_amb) * 100}
Conversion_T.loc[6] = {'temperature':140,'conversion':((Signal_T_amb-Signal_T_140)/Signal_T_amb) * 100}
Conversion_T.loc[7] = {'temperature':160,'conversion':((Signal_T_amb-Signal_T_160)/Signal_T_amb) * 100}
Conversion_T.loc[8] = {'temperature':180,'conversion':((Signal_T_amb-Signal_T_180)/Signal_T_amb) * 100}
Conversion_T.loc[9] = {'temperature':200,'conversion':((Signal_T_amb-Signal_T_200)/Signal_T_amb) * 100}
Conversion_T.loc[10] = {'temperature':220,'conversion':((Signal_T_amb-Signal_T_220)/Signal_T_amb) * 100}
Conversion_T.loc[11] = {'temperature':240,'conversion':((Signal_T_amb-Signal_T_240)/Signal_T_amb) * 100}
Conversion_T.loc[12] = {'temperature':260,'conversion':((Signal_T_amb-Signal_T_260)/Signal_T_amb) * 100}
Conversion_T.loc[13] = {'temperature':280,'conversion':((Signal_T_amb-Signal_T_280)/Signal_T_amb) * 100}
Conversion_T.loc[14] = {'temperature':300,'conversion':((Signal_T_amb-Signal_T_300)/Signal_T_amb) * 100}
Conversion_T.loc[15] = {'temperature':320,'conversion':((Signal_T_amb-Signal_T_320)/Signal_T_amb) * 100}

# Averaging signal for M31 (MeOH) ---> MeOH production vs T


Signal_T_M28 = pd.DataFrame()
Signal_T_M28['temperature'] = []
Signal_T_M28['mean'] = []

Signal_T_M28.loc[0] = {'temperature':25,'mean':Signal_T_amb}
Signal_T_M28.loc[1] = {'temperature':40,'mean':Signal_T_40}
Signal_T_M28.loc[2] = {'temperature':60,'mean':Signal_T_60}
Signal_T_M28.loc[3] = {'temperature':80,'mean':Signal_T_80}
Signal_T_M28.loc[4] = {'temperature':100,'mean':Signal_T_100}
Signal_T_M28.loc[5] = {'temperature':120,'mean':Signal_T_120}
Signal_T_M28.loc[6] = {'temperature':140,'mean':Signal_T_140}
Signal_T_M28.loc[7] = {'temperature':160,'mean':Signal_T_160}
Signal_T_M28.loc[8] = {'temperature':180,'mean':Signal_T_180}
Signal_T_M28.loc[9] = {'temperature':200,'mean':Signal_T_200}
Signal_T_M28.loc[10] = {'temperature':220,'mean':Signal_T_220}
Signal_T_M28.loc[11] = {'temperature':240,'mean':Signal_T_240}
Signal_T_M28.loc[12] = {'temperature':260,'mean':Signal_T_260}
Signal_T_M28.loc[13] = {'temperature':280,'mean':Signal_T_280}
Signal_T_M28.loc[14] = {'temperature':300,'mean':Signal_T_300}
Signal_T_M28.loc[15] = {'temperature':320,'mean':Signal_T_320}

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

# ax1.plot(df["M2-x"], df["M2-y"], label ='H$_2$, M2', linewidth = 0.7)
# ax1.plot(df["M4-x"], df["M4-y"], label ='He, M4', linewidth = 0.7)
# ax1.plot(df["M15-x"], df["M15-y"], label ='CH$_4$, M15', linewidth = 0.7)
# ax1.plot(df["M18-x"], df["M18-y"], label ='H$_2$O, M18', linewidth = 0.7)
# ax1.plot(df["M28-x"], df["M28-y"], label ='CO, M28', linewidth = 0.7)
# ax1.plot(df["M31-x"], df["M31-y"], label ='CH$_3$OH, M31', linewidth = 0.7)
# ax1.plot(df["M32-x"], df["M32-y"], label ='O$_2$, M32', linewidth = 0.7)
# ax1.plot(df["M40-x"], df["M40-y"], label ='Ar, M40', linewidth = 0.7)
# ax1.plot(df["M44-x"], df["M44-y"], label ='CO$_2$, M44', linewidth = 0.7)
# ax1.plot(df["M45-x"], df["M45-y"], label ='C$_2$H$_6$OH, M45', linewidth = 0.7)


#___________________________BASELINED data plot_______________________________#

ax1.plot(df["M28-x"], df["M28-y-new"], label ='M28, CO', linewidth = 1.5, c ='#C17E25')
ax1.plot(df["M32-x"], df["M32-y-new"], label ='M32, O$_2$', linewidth = 1.5, c ='#527D8E')
# ax1.plot(df["M40-x"], df["M40-y-new"], label ='M40, Ar', linewidth = 1, c ='#9653B0')
ax1.plot(df["M44-x"], df["M44-y-new"], label ='M44, CO$_2$', linewidth = 1.5, c ='#8EAE7F')
plt.gca().get_yaxis().set_ticklabels([])

#____________________Right part of the graph (Temperature)____________________#
#_________________comment/uncomment for RTD or TC controlled__________________#

ax2 = ax1.twinx() 

#_____________________________if RTD controlled_______________________________#

ax2.plot(df["RTD temperature-x"], df["RTD temperature-y"], label ='T [°C]',
                color = '#CB434F', linewidth = 1, alpha=1,
                linestyle=linestyles_dict['densely dotted'])

#______________________________if TC controlled_______________________________#

# df["TC temperature-y"] = 1.5054704595186*df["TC temperature-y"]-13.306783369803
# ax2.plot(df["TC temperature-x"], df["TC temperature-y"], label ='T [°C]', 
#                 color = 'tab:red', linewidth = 0.5, alpha=1, 
#                 linestyle=linestyles_dict['dotted'])

# Plotting settings

# plt.title('MRFR13 - 3,5 nm AuTi (cluster source), CO oxidation at 1 bar, 15-03-2023', 
#             fontsize=10, ha="center", fontweight="bold")

# ax1.set_xlabel('Time [$\it{h}$]', fontsize = 10) 
# ax1.set_ylabel('Flow $\it{[nmol/min]}$', fontsize = 10)
ax1.tick_params(axis ='y', labelcolor = 'black') 
# ax1.legend(loc="best", fontsize = 'x-small', facecolor="white",
#            edgecolor="white", framealpha= 1)
# ax1.legend(fontsize = 'small', loc='center left', bbox_to_anchor=(1.15, 0.5),
#           fancybox=True, shadow=True, ncol=1)
# ax1.set_yscale('log')
ax1.set_ylim(6e-2, 15)
ax1.set_xlim(0, 30)
# ax1.grid()
plt.yticks(fontsize = 10)
plt.gca().get_xaxis().set_ticklabels([])

# ax2.set_ylabel('Temperature $\it{(°C)}$', color = 'tab:red', fontsize=10) 
ax2.tick_params(axis ='y', labelcolor = 'tab:red')
# ax2.legend(loc="lower right", fontsize = 'small', facecolor="white", 
#             edgecolor="white", framealpha= 1)
# ax2.set_yscale('log')
ax2.set_ylim(0, 350)
# ax2.grid()
plt.gca().get_yaxis().set_ticklabels([])

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
ax2.legend(fontsize = 'x-small',  bbox_to_anchor=(0.283, 0.95),
          fancybox=False, shadow=False, frameon=False, ncol=1)

# plt.title('CO Conversion', fontsize=10, ha="center", fontweight="bold")
ax1.legend(fontsize = 'x-small', bbox_to_anchor=(0.28, 0.15),
          fancybox=False, shadow=False, frameon=False, ncol=1)

ax1.set_xlabel('Temperature $\it{[°C]}$', fontsize = 10) 
ax1.set_ylabel('Conversion, $\chi$ [%]', fontsize = 10)
ax2.set_ylabel('Flow $\it{[nmol/min]}$', color = 'tab:red', fontsize=10) 
ax2.tick_params(axis ='y', labelcolor = 'tab:red')
# ax2.set_ylim(0, 510)
# ax2.grid()
plt.show()

# =========================================================================== #
#                                 File saving                                 #
# =========================================================================== #

# fig1.savefig('C:/Users/filro/OneDrive - Danmarks Tekniske Universitet/Skrivebord/Data/QMS/QMS_Images/MRFR14_AuTi_3,5nm_ClusterSource_CO_Oxidation_1bar_RTDcontrolled_2to1_O2toCO_FullExp_15032023.tif',
#             format='tif', dpi=1200, bbox_inches="tight")

# fig2.savefig('C:/Users/filro/OneDrive - Danmarks Tekniske Universitet/Skrivebord/Data/QMS/QMS_Images/MRFR14_AuTi_3,5nm_ClusterSource_CO_Oxidation_1bar_RTDcontrolled_2to1_O2toCO_Conversion.png',
#             format='png', dpi=1200, bbox_inches="tight")

# df.to_csv('C:/Users/filro/OneDrive - Danmarks Tekniske Universitet/Skrivebord/Data/QMS/QMS_ProcessedData/ProcessedData_MRFR4_NiGa3_2nmSputtering_CO2_reduction_Methanol_250C_1bar_Activation250.txt', header=True)
