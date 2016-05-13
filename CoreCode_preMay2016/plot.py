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


Pc2Code = Params[-2]
Sol2Code = Params[-1]


r = data[0,:]
z = data[1,:]

Rpos = data[2,:]
Zpos = data[3,:]

idx = int(sys.argv[1])

Bfield = 'N'
NBCont = 2
if( len(sys.argv) > 2 ):
    if(sys.argv[2].upper()=='Y'):
        Bfield = 'Y'
        if( len(sys.argv) > 3 ): NBCont = int(sys.argv[3])

if(idx==8): PlotMe = multiply(Rpos,data[5,:])
else: PlotMe = data[idx,:]

plotmin = min(PlotMe)
plotmax = max(PlotMe)
plot2min = min(n for n in PlotMe if n!=plotmin)
plot2max = min(n for n in PlotMe if n!=plotmax)

DeltaR = data[2,-1]/(max(r)+1)
DeltaZ = data[3,-1]/(max(z)+1)

RGrid = arange(max(r)+2)*DeltaR ## Grid is the edges, not the center of cells
ZGrid = arange(max(z)+2)*DeltaZ

PrepPlot = PlotMe.reshape(max(r)+1,max(z)+1).T  # Need the transpose (.T)


# IF plotting B fields, get some contours
if(Bfield=='Y'):
    PhiData = multiply(Rpos,data[5,:])
    BContours = []
    DeltaPhi = (max(PhiData)-min(PhiData))/1.0/NBCont
    #print DeltaPhi
    for i in range(NBCont):
        curPhi = DeltaPhi*(i+1.0)
        print curPhi
        curCont = []
        for j in range(int(Params[1])): #for each Z row
            curRow = [ PhiData[k] for k in range(int(Params[0]*Params[1])) if (k%(int(Params[1]))==j) ] #PhiData[ Params[1]*i : Params[1]*(i+1) ]
            print curRow
            leftidx = 0
            for k in range(len(curRow)-1):
                if(( curPhi <= curRow[k+1] and curPhi >= curRow[k] ) or ( curPhi >= curRow[k+1] and curPhi <= curRow[k] )):
                    leftidx = k
                    break
            lPhi = curRow[leftidx]
            rPhi = curRow[leftidx+1]
            lRos = DeltaR*(leftidx)
            rRos = DeltaR*(leftidx+1)
            m = (rPhi-lPhi)/DeltaR
            curCont.append( (curPhi-lPhi)/m + lRos )
        BContours.append(curCont)




#plt.setp(ax,xticklabels=[]) #turns off labels

# ---===---===---===---===---===---===---===---===---===
plt.subplot(1,1,1)
plt.subplots_adjust(left=0.18,bottom=0.15,wspace=0.0001,hspace=0.0001)

#plt.axes([0.16,0.13,0.79,0.82]) #left,bot,wid,height
plt.axis([0, data[2,-1]/Pc2Code, 0, data[3,-1]/Pc2Code])

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
if(idx==8): plt.title(r'Phi')


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

cmapmin = plotmin
cmapmax = plotmax
#cmapmin = 0.80
#cmapmax = 1.2

cmap = 'Paired'
cmap = 'RdBu_r'
#cmap = 'Spectral'
#cmap = 'Blues_r'
plt.pcolor(RGrid/Pc2Code,ZGrid/Pc2Code,PrepPlot,cmap=cmap,vmin=cmapmin,vmax=cmapmax)

plt.colorbar()


# Filament Boundary
area = pi*2.0
plt.scatter(VCont/Pc2Code,(ZGrid[:-1]+0.5*DeltaZ)/Pc2Code,s=1.5*area,c='black')
plt.plot(VCont/Pc2Code,(ZGrid[:-1]+0.5*DeltaZ)/Pc2Code,c='black',linewidth=1,alpha=0.3)


# B fields
if(Bfield=='Y'):
    for i in range(NBCont):
        plt.plot(BContours[i]/Pc2Code,(ZGrid[:-1]+0.5*DeltaZ)/Pc2Code,c='black',linewidth=2,alpha=0.7)
        #plt.scatter(BContours[i],ZGrid[:-1]+0.5*DeltaZ,s=2*1.5*area,c='blue',marker=(4,1))

#Text
#plt.text(0.05,0.06,r'${\cal M} = 4.47$',fontsize=26.0)

#plt.legend(loc=4,numpoints=1,ncol=2,labelspacing=0.1,columnspacing=0.9)


# Creates Plot
#plt.show()
plt.savefig('Plot.pdf')
#plt.savefig('figure3.eps')

#Restores defaults
mpl.rcdefaults()
