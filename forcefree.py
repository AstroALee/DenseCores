# Imports
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *
import scipy as spi
import scipy.special as special
import sys

# Colors
dred = [0.6,0,0]


Num = int(sys.argv[1])

Len = float(sys.argv[2])

BesselZero = 2.4048255
SolnK = BesselZero/Len
SolnA = float(sys.argv[3])

DeltaX = Len/(1.0+float(Num))

r = zeros(Num*Num)
z = zeros(Num*Num)
plot = zeros(Num*Num)

for i in range(0,Num*Num):
    r[i] = i/Num
    z[i] = i%Num
    plot[i] = SolnA*special.jv(0,SolnK*DeltaX*r[i])*cosh(SolnK*DeltaX*z[i])


plotmin = min(plot)
plotmax = max(plot)
print plotmin
print plotmax

RGrid = arange(max(r)+2)*DeltaX ## Grid is the edges, not the center of cells
ZGrid = arange(max(z)+2)*DeltaX

PlotPlot = plot.reshape(max(r)+1,max(z)+1).T # Need the transpose (.T)



#plt.setp(ax,xticklabels=[]) #turns off labels

# ---===---===---===---===---===---===---===---===---===
plt.subplot(1,1,1)
plt.subplots_adjust(left=0.18,bottom=0.15,wspace=0.0001,hspace=0.0001)

#plt.axes([0.16,0.13,0.79,0.82]) #left,bot,wid,height
plt.axis([0, Len-DeltaX, 0, 1.0*Len-DeltaX])

# Axis labels
plt.xlabel(r'Distance    $R$') #TeX requires the r
plt.ylabel(r'Distance    $Z$')

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
plt.pcolor(RGrid,ZGrid,(PlotPlot),cmap='Paired',vmin=(plotmin),vmax=(plotmax))

plt.colorbar()

plt.title(r'$A(r,z) = A_0 J_0(kr) \cosh(kz)$')
# Creates Plot
plt.show()
#plt.savefig('Plot.pdf')
#plt.savefig('figure3.eps')

#Restores defaults
mpl.rcdefaults()
