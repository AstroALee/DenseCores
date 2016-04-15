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
plt.xlabel('Number of Finite Difference Solves',fontsize=20)
plt.ylabel('Error Measure',fontsize=20)
plt.title('Convergence',fontsize=20)
#plt.axis([0,15,-8,10])

# Custom major tick locations and labels?
# ==============================
#plt.gca().set_xticks([2,5,7])
#plt.gca().set_xticklabels(['a','b','c'])
#plt.gca().set_yticks([2,5,7])
#plt.gca().set_yticklabels(['a','b','c'])


# Do all your data plotting here
# ==================

# Read in data from a file?
data = np.genfromtxt('Errors.out',unpack=True)

# Assign data to appropriate variables
x = data[1,:] + 1  # + 1 if you start at 0, which means there's only been one solve
v = data[2,:]
a = data[3,:]


plt.semilogy(x,v,'--',label=r'Gravitational Potential Gradient')
plt.semilogy(x,a,label=r'Magnetic Potential Gradient')

plt.legend(loc=1,numpoints=1,ncol=2,labelspacing=0.1,columnspacing=0.9)



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
plt.savefig("ErrorMeasure.png")
