# Imports
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *
import scipy as spi
import sys

# Colors
dred = [0.6,0,0]

# Reads in Data as column arrays
data = genfromtxt('test.txt',delimiter=",",unpack=True)

#doesn't include last end point in array. If you want the first five
#items, do 0:6.
r = data[0,:]
z = data[1,:]
PlotMe = data[2,:]

pLength = float(sys.argv[1])

DeltaR = pLength/(max(r)+1.0)
DeltaZ = pLength/(max(z)+1.0)


RGrid = arange(max(r)+2)*DeltaR ## Grid is the edges, not the center of cells
ZGrid = arange(max(z)+2)*DeltaZ

PlotPlot = PlotMe.reshape(max(r)+1,max(z)+1).T # Need the transpose (.T)

plotmin = min(PlotMe)
plotmax = max(PlotMe)
plot2min = min(n for n in PlotMe if n!=plotmin)


#plt.setp(ax,xticklabels=[]) #turns off labels

# ---===---===---===---===---===---===---===---===---===
plt.subplot(1,1,1)
plt.subplots_adjust(left=0.18,bottom=0.15,wspace=0.0001,hspace=0.0001)

#plt.axes([0.16,0.13,0.79,0.82]) #left,bot,wid,height
plt.axis([0, pLength, 0, 1.0*pLength])

# Axis labels
plt.xlabel(r'Distance    $R$') #TeX requires the r
plt.ylabel(r'Distance    $Z$')
plt.title(r'Gravitational Potential')
#Moves the x-axis down a little
plt.gca().xaxis.labelpad=9 #5 is default

#Makes ticks thicker/longer
ax = plt.gca()
ax.tick_params(which='both',width=2)
ax.tick_params(which='major',size=12)
ax.tick_params(which='minor',size=7)
	 
#moves numbers out a little
ax.tick_params(axis='y',pad=9)
ax.tick_params(axis='x',pad=9)

# Adds dashed lines
#xline = [-5,5]
#yline = [2,2]
#plt.plot(xline,yline,'k:',linewidth=2.0)




plt.pcolor(RGrid,ZGrid,log10(PlotPlot),cmap='Paired',vmin=log10(plot2min),vmax=log10(plotmax))
plt.colorbar()


Cont = array([0.484375,0.483296,0.481134,0.477885,0.473559,0.468187,0.461843,0.454656,0.447168,0.439292,0.431268,0.423517,0.416623,0.410899,0.40677,0.404603])

area = pi*2.0


plt.scatter(Cont,ZGrid[:-1],s=area,c='black')


#Text
#plt.text(0.05,0.06,r'${\cal M} = 4.47$',fontsize=26.0)

#plt.legend(loc=4,numpoints=1,ncol=2,labelspacing=0.1,columnspacing=0.9)


# Creates Plot
plt.show()
#plt.savefig('Plot.pdf')
#plt.savefig('figure3.eps')

#Restores defaults
mpl.rcdefaults() 