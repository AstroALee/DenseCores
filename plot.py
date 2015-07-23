# Imports
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *
import scipy as spi

# Colors
dred = [0.6,0,0]

# Reads in Data as column arrays
data = genfromtxt('test.txt',delimiter=",",unpack=True)

#doesn't include last end point in array. If you want the first five
#items, do 0:6.
r = data[0,:]
z = data[1,:]
PlotMe = data[2,:]


RGrid = arange(max(r)+2) ## Grid is the edges, not the center of cells
ZGrid = arange(max(z)+2)
PlotPlot = PlotMe.reshape(max(r)+1,max(z)+1).T

plotmin = min(PlotMe)
plotmax = max(PlotMe)
plot2min = min(n for n in PlotMe if n!=plotmin)


#plt.setp(ax,xticklabels=[]) #turns off labels

# ---===---===---===---===---===---===---===---===---===
plt.subplot(1,1,1)
plt.subplots_adjust(left=0.18,bottom=0.15,wspace=0.0001,hspace=0.0001)

#plt.axes([0.16,0.13,0.79,0.82]) #left,bot,wid,height
plt.axis([min(r), max(r)+1, min(z), max(z)+1])

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


Cont = array([50.5,50.4886,50.4658,50.4316,50.386,50.329,50.2605,50.1807,50.0894,49.9867,49.8726,49.7472,49.6105,49.4637,49.309,49.1435,48.9673,48.7806,48.5836,48.3796,48.1679,47.9469,47.7169,47.4789,47.2375,46.9887,46.733,46.4719,46.2095,45.9426,45.6721,45.4006,45.13,44.8589,44.5887,44.3227,44.0609,43.804,43.5533,43.3119,43.08,42.8587,42.6496,42.4542,42.2745,42.1109,41.9645,41.8365,41.7279,41.6392,41.5713])

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