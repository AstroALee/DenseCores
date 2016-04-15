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
plt.xlabel(r'Num of grid cells     $N_r$ ',fontsize=20)
plt.ylabel('Mass Error',fontsize=20)
plt.title('Integration Measures',fontsize=20)
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
N = np.array([21,41,81,161,321,641,1281])

errC21 = np.array([0.002889857612242387,0.00282143447696593,0.00278721614076634,0.0027701051750394,0.002761549269718954,0.002757271204657558 ,0.002755132145731898  ])
errT21 = np.array([0.002621862932491298,0.002621889248002543,0.002621895959509248, 0.002621897603236747, 0.002621898022760518, 0.002621898125473687 , 0.002621898151692177  ])
errS21 = np.array([0.002621898454463879,0.002621898130959228,0.002621898179989353, 0.002621898158354433 , 0.002621898161542243, 0.002621898160168936 , 0.002621898160366057  ])


## Number of N_z = 81, R/Z ratio = 1.0
errC81 = np.array([0.002786648411805127,0.00272066895993146,0.002687672707167593,0.002671172847359429,0.002662922510086135,0.00265879723306273 ,0.002656734569098628  ])
errT81 = np.array([0.002621862932491267, 0.002621889248002568, 0.002621895959509276, 0.00262189760323678, 0.002621898022760515, 0.002621898125473749 , 0.002621898151692182   ])
errS81 = np.array([0.002621898454463866, 0.002621898130959247, 0.002621898179989344 , 0.002621898158354458, 0.002621898161542252, 0.002621898160168988, 0.002621898160366034 ])

errC321 = np.array([0.00276084611169579,0.002695477580672843,0.002662786848767783,0.002646439765439329,0.002638265820177967,0.002634178740164309,0.002632135174940175 ])
errT321 = np.array([ 0.002621862932491307,0.002621889248002618, 0.002621895959509125, 0.002621897603236693,0.002621898022760564, 0.002621898125474131,  0.00262189815169206  ])
errS321 = np.array([0.002621898454463839,  0.0026218981309592,0.002621898179989235,0.00262189815835433, 0.002621898161542187, 0.002621898160169221, 0.002621898160366392   ])






#Mathematica value
analval = 0.0026246309302537037

errC21 = np.fabs(errC21-analval)/analval
errT21 = np.fabs(errT21-analval)/analval
errS21 = np.fabs(errS21-analval)/analval

errC81 = np.fabs(errC81-analval)/analval
errT81 = np.fabs(errT81-analval)/analval
errS81 = np.fabs(errS81-analval)/analval

errC321 = np.fabs(errC321-analval)/analval
errT321 = np.fabs(errT321-analval)/analval
errS321 = np.fabs(errS321-analval)/analval

# Mathematica solution
#plt.semilogx(N,0.0026246309302537037 + 0*N,label=r'Mathematica')

if(1):
    plt.loglog(N,errC81,'-',label=r'Riemann')
    plt.loglog(N,errT81,'-',label=r'Trapezoid')
    plt.loglog(N,errS81,'-',label=r'Simpson')
    plt.loglog(N,errC21,'--')
    plt.loglog(N,errT21,'--')
    plt.loglog(N,errS21,'--')
    plt.loglog(N,errC321,'-.')
    plt.loglog(N,errT321,'-.')
    plt.loglog(N,errS321,'-.')
else:
    plt.semilogx(N,errC81,'-',label=r'Riemann')
    plt.semilogx(N,errT81,'-',label=r'Trapezoid')
    plt.semilogx(N,errS81,'-',label=r'Simpson')
    plt.semilogx(N,errC21,'--')
    plt.semilogx(N,errT21,'--')
    plt.semilogx(N,errS21,'--')
    plt.semilogx(N,errC321,'-.')
    plt.semilogx(N,errT321,'-.')
    plt.semilogx(N,errS321,'-.')



plt.legend(loc=3,numpoints=1,ncol=2,labelspacing=0.1,columnspacing=0.9,fontsize=15)

plt.text(500,0.06,r"$N_z=21$",fontsize=15)
plt.text(500,0.017,r"$N_z=81$",fontsize=15)
plt.text(500,0.005,r"$N_z=321$",fontsize=15)

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
plt.savefig("IntegralMeasures.png")
