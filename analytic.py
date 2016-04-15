# Data and Seaborn Packages (always load)
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import scipy as spi


# Perhaps here you import a 'header' file that sets up how
# you want your plots to look
# from localfile import *

# Adjusts region plot takes up in figure
plt.subplots_adjust(left=0.15,bottom=0.15,wspace=0.0001,hspace=0.0001)

# Set the general style
# =====================
#sns.set_style("darkgrid") # gray background, white grid (no ticks)
sns.set_style("whitegrid") # white background, gray grid (no ticks)
#sns.set_style("dark") # gray background, no grid or ticks or border
#sns.set_style("white") # white background and black border, no grid or ticks
#sns.set_style("ticks") # white background, tick marks on black axes

# Axis styles
# ================
sns.set_style( {'axes.linewidth':2 } ) # thicken lines
sns.set_context( {'xtick.major.width': 2.0,'ytick.major.width': 2.0})

# Axis and Tick colors
sns.set_style( { 'axes.edgecolor': '0', 'xtick.color': '0', 'ytick.color': '0' } )
sns.set_style( {'grid.color': '.85' })

# Tick marks
# ==================
sns.set_style( {'xtick.major.size':8 , 'ytick.major.size':8})
sns.set_style( {'xtick.direction':u'in', 'ytick.direction':u'in'})
sns.set_context( {'xtick.major.pad': 11.0 , 'xtick.major.pad': 11.0 })
sns.set_context( {'xtick.labelsize': 15.0, 'ytick.labelsize': 15.0} )

# Color palette
# =================
sns.set_palette('deep') # inputs any seaborn or matplotlib color map
#sns.set_palette(sns.dark_palette("purple"))

# Labels and Such
plt.xlabel(r'Radius ',fontsize=20)
plt.ylabel('Density',fontsize=20)
#plt.title('Integration Measures',fontsize=20)
#plt.axis([10,1e4,0.001,0.00106])

# Custom major tick locations and labels?
# ==============================
#plt.gca().set_xticks([2,5,7])
#plt.gca().set_xticklabels(['a','b','c'])
#plt.gca().set_yticks([2,5,7])
#plt.gca().set_yticklabels(['a','b','c'])


# Do all your data plotting here
# ==================

# Read in data from a file?
#data = np.genfromtxt('Errors.out',unpack=True)

# Assign data to appropriate variables
#x = data[1,:] + 1  # + 1 if you start at 0, which means there's only been one solve
#v = data[2,:]
#a = data[3,:]

## Number of N_z = 21, R/Z ratio = 1.0
RhoTop = np.array([1.164459085671002,1.164343565289169,1.16399855270063,1.163423172606706,1.162618854507716,1.161585658421739,1.160324890436914,1.158837536628536,1.157124783152768,1.155188519839381,1.153029810265238,1.150651428863217,1.148054326774537,1.145241996120269,1.142215642040211,1.138978959879492,1.135533887185871,1.131883952108906,1.128031889879728,1.12398105726374,1.119734950291426,1.115296767373707,1.110670711043054,1.105859822148996,1.100868933662268,1.095700943760192,1.090361073780069,1.084852669634858,1.079180713801885,1.07334913301802,1.067362672033316,1.061225785702046,1.054943001701583,1.048519230419118,1.041958791809162,1.035266967177879,1.028447885217112,1.021507076282644,1.01444869766812,1.007278200755853,0.9999999997830338])
Radius = np.array([0,0.006610643947204104,0.01322128789440821,0.01983193184161231,0.02644257578881641,0.03305321973602052,0.03966386368322462,0.04627450763042872,0.05288515157763283,0.05949579552483693,0.06610643947204103,0.07271708341924514,0.07932772736644925,0.08593837131365335,0.09254901526085745,0.09915965920806155,0.1057703031552657,0.1123809471024698,0.1189915910496739,0.125602234996878,0.1322128789440821,0.1388235228912862,0.1454341668384903,0.1520448107856944,0.1586554547328985,0.1652660986801026,0.1718767426273067,0.1784873865745108,0.1850980305217149,0.191708674468919,0.1983193184161231,0.2049299623633272,0.2115406063105313,0.2181512502577354,0.2247618942049395,0.2313725381521436,0.2379831820993477,0.2445938260465518,0.251204469993756,0.25781511394096,0.2644257578881641])
plt.plot(Radius,RhoTop,label=r'Code')


plt.plot(Radius,1.0+0*Radius,'k--')
# Analytic Solution

pi = 3.14159265359
beta = 1.0
Dcst = (pi*RhoTop[0]/2.0)/(1.0+1.0/beta)

print Dcst

#debug
#Radius = np.array(np.linspace(0,0.75,100))

Rho_anal = RhoTop[0]/(1.0+Dcst*(Radius)**2)**2
plt.plot(0.9*Radius,Rho_anal,label=r'Scaled Analytic')
plt.plot(Radius,Rho_anal,label=r'Analytic')



#debug
RhoTest = np.array([1.007821337150717,1.007814111613822,1.007792449295447,1.007756351595061,1.007705820817345,1.007640860232731,1.007561474050048,1.007467667388549,1.007359446277891,1.007236817657905,1.007099789379787,1.006948359883281,1.006782546022668,1.006602361857024,1.006407818988196,1.006198929923675,1.005975708076599,1.005738167765751,1.00548632421556,1.005220193555945,1.00493979278336,1.004645132929195,1.004336233314915,1.004013120449186,1.003675815055961,1.003324338735363,1.002958713963695,1.002578964093428,1.002185113353212,1.001777186847869,1.001355210558394,1.000919207753025,1.000469199977104,1.000005225673208,1,1,1,1,1,1])
#plt.plot(Radius,RhoTest,label=r'Test')



plt.legend(loc=1,numpoints=1,ncol=2,labelspacing=0.1,columnspacing=0.9,fontsize=15)

#plt.text(500,0.06,r"$N_z=21$",fontsize=15)

# Post plotting alterations
# =============================================================
sns.despine(offset=0,trim=False,left=False,right=True,top=True,bottom=False)
            # for white and ticks styles, removes axes marked as True.
            # offset = pushes axes away from data
            # trim = limits range of surviving spines
            # Ticks and their labels still shown

# Production
# ==================
#plt.show()
plt.savefig("Analytic.png")
