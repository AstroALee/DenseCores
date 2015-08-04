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
data2 = genfromtxt('test2.txt',delimiter=",",unpack=True)

dataC = genfromtxt('VContour.txt',delimiter=",",unpack=True)
Cont = dataC[:-1]

#array([0.97619,0.974944,0.972449,0.968704,0.963713,0.957488,0.950057,0.941474,0.931819,0.921583,0.910675,0.899113,0.887149,0.875269,0.863747,0.852763,0.842738,0.834093,0.827274,0.822515,0.82006])


idx = int(sys.argv[1])

#doesn't include last end point in array. If you want the first five
#items, do 0:6.
r = data[0,:]
z = data[1,:]


if(idx==5):
    Gpot = data[2,:]
    Q = data[4,:]
    ePot = array([exp(-g) for g in Gpot])
    PlotMe = multiply(Q,ePot)
    Gpot2 = data2[2,:]
    Q2 = data2[4,:]
    ePot2 = array([exp(-g) for g in Gpot2])
    PlotMe2 = multiply(Q2,ePot2)
else:
    PlotMe = data[idx,:]
    PlotMe2 = data2[idx,:]


plotmin = min(PlotMe)
plotmax = max(PlotMe)
plotmin2 = min(PlotMe2)
plotmax2 = max(PlotMe2)
print plotmin
print plotmax
print plotmin2
print plotmax2


# If potential or Q, adjust so the min value is 0.
if(idx==2 or idx==4): PlotMe = array([y - plotmin for y in PlotMe])
if(idx==2 or idx==4): PlotMe2 = array([y - plotmin2 for y in PlotMe2])



pLength = float(sys.argv[2])

DeltaR = pLength/(max(r)+1.0)
DeltaZ = pLength/(max(z)+1.0)


RGrid = arange(max(r)+2)*DeltaR ## Grid is the edges, not the center of cells
ZGrid = arange(max(z)+2)*DeltaZ


diff = 100*divide( subtract(PlotMe2,PlotMe) , PlotMe2 )
where_are_nans = isnan(diff)
diff[where_are_nans] = 0

Plotdiff = diff.reshape(max(r)+1,max(z)+1).T # Need the transpose (.T)

diffmin = min(diff)
diffmax = max(diff)


#plt.setp(ax,xticklabels=[]) #turns off labels

# ---===---===---===---===---===---===---===---===---===
plt.subplot(1,1,1)
plt.subplots_adjust(left=0.18,bottom=0.15,wspace=0.0001,hspace=0.0001)

#plt.axes([0.16,0.13,0.79,0.82]) #left,bot,wid,height
plt.axis([0, pLength, 0, 1.0*pLength])

# Axis labels
plt.xlabel(r'Distance    $R$') #TeX requires the r
plt.ylabel(r'Distance    $Z$')


if(idx==2): plt.title(r'Gravitational Potential Diff %')
if(idx==3): plt.title(r'A Diff %')
if(idx==4): plt.title(r'Q Diff %')
if(idx==5): plt.title(r'Density Diff %')


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





#plt.pcolor(RGrid,ZGrid,log10(PlotPlot),cmap='Paired',vmin=log10(plot2min),vmax=log10(plotmax))
plt.pcolor(RGrid,ZGrid,Plotdiff,cmap='Paired',vmin=min(diff),vmax=max(diff))

plt.colorbar()


area = pi*2.0


plt.scatter(Cont,ZGrid[:-1]+0.5*DeltaZ,s=1.5*area,c='black')
plt.plot(Cont,ZGrid[:-1]+0.5*DeltaZ,c='black',linewidth=1,alpha=0.3)


#Text
#plt.text(0.05,0.06,r'${\cal M} = 4.47$',fontsize=26.0)

#plt.legend(loc=4,numpoints=1,ncol=2,labelspacing=0.1,columnspacing=0.9)


# Creates Plot
plt.show()
#plt.savefig('Plot.pdf')
#plt.savefig('figure3.eps')

#Restores defaults
mpl.rcdefaults() 