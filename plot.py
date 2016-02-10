# Imports
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *
import scipy as spi
import sys

# Colors
dred = [0.6,0,0]

# Reads in raw data as column arrays
data = genfromtxt('Output.out',delimiter=",",unpack=True,skip_header=2)

# Read in contour array
flen = len(data[0,:])
VCont = genfromtxt('Output.out',delimiter=",",unpack=True,skip_header=1,skip_footer=flen)

# M, N, zL (pc), rRatio (entire box), mExcess (sol), beta, n
Params = genfromtxt('Output.out',delimiter=",",unpack=True,skip_footer=flen+1)

r = data[0,:]
z = data[1,:]

idx = int(sys.argv[1])
PlotMe = data[idx,:]

plotmin = min(PlotMe)
plotmax = max(PlotMe)
plot2min = min(n for n in PlotMe if n!=plotmin)
plot2max = min(n for n in PlotMe if n!=plotmax)

DeltaR = data[2,-1]/(max(r)+1)
DeltaZ = data[3,-1]/(max(z)+1)

RGrid = arange(max(r)+2)*DeltaR ## Grid is the edges, not the center of cells
ZGrid = arange(max(z)+2)*DeltaZ

PrepPlot = PlotMe.reshape(max(r)+1,max(z)+1).T  # Need the transpose (.T)


#plt.setp(ax,xticklabels=[]) #turns off labels

# ---===---===---===---===---===---===---===---===---===
plt.subplot(1,1,1)
plt.subplots_adjust(left=0.18,bottom=0.15,wspace=0.0001,hspace=0.0001)

#plt.axes([0.16,0.13,0.79,0.82]) #left,bot,wid,height
plt.axis([0, data[2,-1], 0, data[3,-1]])

# Axis labels
plt.xlabel(r'Distance   $R$  (pc)') #TeX requires the r
plt.ylabel(r'Distance   $Z$  (pc)')

if(idx==0): plt.title(r'R Grid id')
if(idx==1): plt.title(r'Z Grid id')
if(idx==2): plt.title(r'R Distance (pc)')
if(idx==3): plt.title(r'Z Distance (pc)')
if(idx==4): plt.title(r'Gravitational Potential')
if(idx==5): plt.title(r'A')
if(idx==6): plt.title(r'Q')
if(idx==7): plt.title(r'Density')


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


# Overwrite min max plotted
#plotmin = 0
#plotmax = 20


#plt.pcolor(RGrid,ZGrid,log10(PlotPlot),cmap='Paired',vmin=log10(plot2min),vmax=log10(plotmax))
plt.pcolor(RGrid,ZGrid,PrepPlot,cmap='Paired',vmin=(plotmin),vmax=(plotmax))

plt.colorbar()


# Filament Boundary
area = pi*2.0
plt.scatter(VCont,ZGrid[:-1]+0.5*DeltaZ,s=1.5*area,c='black')
plt.plot(VCont,ZGrid[:-1]+0.5*DeltaZ,c='black',linewidth=1,alpha=0.3)


#Text
#plt.text(0.05,0.06,r'${\cal M} = 4.47$',fontsize=26.0)

#plt.legend(loc=4,numpoints=1,ncol=2,labelspacing=0.1,columnspacing=0.9)


# Creates Plot
#plt.show()
plt.savefig('Plot.pdf')
#plt.savefig('figure3.eps')

#Restores defaults
mpl.rcdefaults()
