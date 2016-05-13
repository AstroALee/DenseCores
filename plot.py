# Imports
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *
import scipy as spi
import sys

# INPUTS
# ----
# python plot.py A Y/N Y/N #
# A = numerical index for what to plot
# Y/N = unit convert?
# Y/N = plot B-field lines?
# # = number of B-field lines (must be >=2)


# Colors
dred = [0.6,0,0]

# READ IN THE DATA

# Reads in raw data as column arrays
# R cell, Z Cell, R dis, Z dis, V, A, Q, Rho
data = genfromtxt('Output.out',delimiter=",",unpack=True,skip_header=2)
# Read in contour array (line 2)
flen = len(data[0,:])
VCont = genfromtxt('Output.out',delimiter=",",unpack=True,skip_header=1,skip_footer=flen)
# M, N, zL, rRatio (entire box), contR, beta, n, VContour value, lambda, Pc2Code, Sol2Code
# All the physical quantities are non-dimensional
Params = genfromtxt('Output.out',delimiter=",",unpack=True,skip_footer=flen+1)


# DEFINE USEFUL VALUES
Pc2Code = Params[-2]
Sol2Code = Params[-1]
M = Params[0]
N = Params[1]
zL = Params[2]
rL = Params[2]*Params[3]
DeltaR = rL/M
DeltaZ = zL/N


# Index for plotting
idxS = str(sys.argv[1])
idx = -1
if(idxS=='0'): idx = 0
if(idxS=='1'): idx = 1
if(idxS=='2'): idx = 2
if(idxS=='3'): idx = 3
if(idxS=='4' or idxS.lower()=='v'): idx = 4
if(idxS=='5' or idxS.lower()=='a'): idx = 5
if(idxS=='6' or idxS.lower()=='q'): idx = 6
if(idxS=='7' or idxS.lower()=='rho'): idx = 7
if(idxS=='8' or idxS.lower()=='phi'): idx = 8

if(idx==-1):
    print "ERROR with idx"
    exit(1)

Units = 0
if( str(sys.argv[2]).lower()=='y'): Units = 1

# USEFUL ARRAYS
# Position arrays
Rpos = data[2,:]
Zpos = data[3,:]

# Grid arrays, THIS IS WHAT IS USED FOR PLOTTING
RGrid = arange(M+1)*DeltaR # +1 because matplotlib grid is the edges, not the center of cells
ZGrid = arange(N+1)*DeltaZ

# THIS IS WHAT IS PLOTTED
if(idx==8): PlotMe = multiply(Rpos,data[5,:]) # Phi = r*A
else: PlotMe = data[idx,:]
PreparedPlot = PlotMe.reshape(M,N).T  # Need the transpose (.T)


# PlotMe may be a derived quantity, here are some values
plotmin = min(PlotMe)
plot2min = min(n for n in PlotMe if n!=plotmin)
plotmax = max(PlotMe)
plot2max = min(n for n in PlotMe if n!=plotmax)

# POtentially we have some B fields to draw
Bfield = 'N'
NBCont = 1
if( len(sys.argv) > 3 ):
    if(sys.argv[3].upper()=='Y'):
        Bfield = 'Y'
        if( len(sys.argv) > 4 and int(sys.argv[4]) > 1): NBCont = int(sys.argv[4]) # will plot min, at least

# If B fields, plot field lines with equal
BContours = []
if(Bfield=='Y'):
    curR = rL/10.0 #first anchor point
    lastR = 0
    PhiArray = multiply(Rpos,data[5,:]) # we draw contours of constant phi
    botRowR = [ Rpos[k] for k in range(int(M*N)) if (k%(int(N))==0) ] # bottom row R values
    botRowPhi = [ PhiArray[k] for k in range(int(M*N)) if (k%(int(N))==0) ] # bottom row Phi values
    botPhi = interp(curR,botRowR,botRowPhi) # phi value at curR
    lastPhi = botPhi
    goalDelta = botPhi/curR # assumes curPhi(r=0) = 0, goal is cst DeltaPhi/r
    #goalDelta = botPhi      # assumes curPhi(r=0) = 0, goal is cst DeltaPhi
    # now finds the list of r anchor points
    rAnchor=[]
    rAnchor.append(curR)
    for i in range(1,NBCont):
        locDelta = [(x - lastPhi)/(y+0.000000001) for x,y in zip(botRowPhi,botRowR)] # first goal
        #print locDelta
        #locDelta = [(x - lastPhi) for x in botRowPhi] # second goal
        #print locDelta
        # now find new r
        curR = interp(goalDelta,locDelta,botRowR) # if Phi is monotonically increasing, this will always work
        locPhi = interp(curR,botRowR,botRowPhi)
        rAnchor.append(curR)
        print curR
        lastPhi = locPhi
    #With the set of anchor points, now find contours
    for i in range(NBCont):
        curCont = []
        curR = rAnchor[i]
        #using the new value of r (curR), find the phi value on the bottom
        curPhi = interp(curR,botRowR,botRowPhi)
        #now for each row above, find the r value where curPhi occurs (B fields are contours of Phi)
        curCont.append(curR)
        for j in range(1,int(N)): #for each Z row
            locRowPhi = [ PhiArray[k] for k in range(int(M*N)) if (k%(int(N))==j) ] # local row Phi values
            locR = interp(curPhi,locRowPhi,botRowR) # each row as the same set of R values, hence botRowR
            curCont.append(locR)
        BContours.append(curCont)

# UNIT CONVERSIONS
Code2Pc = 1.0/Pc2Code

PlotUnit = 1.0

#plt.setp(ax,xticklabels=[]) #turns off labels

# ---===---===---===---===---===---===---===---===---===
plt.subplot(1,1,1)
plt.subplots_adjust(left=0.18,bottom=0.15,wspace=0.0001,hspace=0.0001)

#plt.axes([0.16,0.13,0.79,0.82]) #left,bot,wid,height
rPlotmax = rL
zPlotmax = zL
if(Units):
    rPlotmax = rPlotmax*Code2Pc
    zPlotmax = zPlotmax*Code2Pc

plt.axis([0, rPlotmax, 0, zPlotmax])

# Axis labels
if(Units):
    plt.xlabel(r'Distance   $R$  (pc)') #TeX requires the r
    plt.ylabel(r'Distance   $Z$  (pc)')
else:
    plt.xlabel(r'Distance   $R$   (nd)') #TeX requires the r
    plt.ylabel(r'Distance   $Z$   (nd)')

# Plot Titles
if(idx==0): plt.title(r'R Grid id')
if(idx==1): plt.title(r'Z Grid id')
if(idx==2): plt.title(r'R Distance')
if(idx==3): plt.title(r'Z Distance')
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

# Max and Min for ColorMap
# Overwrite min max plotted
#plotmin = 0
#plotmax = 20
cmapmin = plotmin
cmapmax = plotmax

# ColorMap choice
cmap = 'Paired'
cmap = 'RdBu_r'
#cmap = 'Spectral'
#cmap = 'Blues_r'

if(Units):
    RGrid = RGrid*Code2Pc
    ZGrid = ZGrid*Code2Pc
    DeltaZ = DeltaZ*Code2Pc
    DeltaR = DeltaR*Code2Pc
    VCont = VCont*Code2Pc
    PreparedPlot = PreparedPlot*PlotUnit
    BContours=array(BContours)*Code2Pc

plt.pcolor(RGrid,ZGrid,PreparedPlot,cmap=cmap,vmin=cmapmin,vmax=cmapmax)
plt.colorbar()

# Filament Boundary
area = pi*2.0
plt.scatter(VCont,(ZGrid[:-1]+0.5*DeltaZ),s=1.5*area,c='black')
plt.plot(VCont,(ZGrid[:-1]+0.5*DeltaZ),c='black',linewidth=1,alpha=0.3)

# B fields
if(Bfield=='Y'):
    for i in range(NBCont):
        plt.plot(BContours[i],(ZGrid[:-1]+0.5*DeltaZ),c='black',linewidth=2,alpha=0.7)

#Text
#plt.text(0.05,0.06,r'${\cal M} = 4.47$',fontsize=26.0)


# Creates Plot
#plt.show()
plt.savefig('Plot.pdf')

#Restores defaults
mpl.rcdefaults()
