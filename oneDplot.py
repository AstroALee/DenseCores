# Imports
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *
import scipy as spi
import sys

# Colors
dmagenta = '#DB1672'
dgreen = '#259B41'
dpurp  = '#BB46C9'

# Comparing 1D data
comp = 1
CompX = [0.0, 0.00265789, 0.0053157899999999999, 0.0079736800000000004, 0.0106316, 0.013289499999999999, 0.0159474, 0.018605300000000002, 0.021263199999999999, 0.023921100000000001, 0.026578899999999999, 0.0292368, 0.031894699999999998, 0.034552600000000003, 0.037210500000000001, 0.039868399999999998, 0.042526300000000003, 0.045184200000000001, 0.047842099999999999, 0.050500000000000003]
CompY = [0.0, 5.1524699999999998e-05, 0.000201258, 0.00044955800000000002, 0.00080049299999999999, 0.00124957, 0.0017974899999999999, 0.0024475399999999998, 0.0031955799999999999, 0.00404266, 0.0049912899999999998, 0.0060376800000000001, 0.0071832099999999998, 0.0084296300000000005, 0.0097734899999999993, 0.011216500000000001, 0.012759700000000001, 0.0143999, 0.016139299999999999, 0.017972700000000001]



# Reads in raw data as column arrays
data = genfromtxt('Output.out',delimiter=",",unpack=True,skip_header=2)

# Read in contour array
flen = len(data[0,:])
VCont = genfromtxt('Output.out',delimiter=",",unpack=True,skip_header=1,skip_footer=flen)

# M, N, zL (pc), rRatio (entire box), mExcess (sol), beta, n
Params = genfromtxt('Output.out',delimiter=",",unpack=True,skip_footer=flen+1)


# What to plot?
idx = int(sys.argv[1])

# Slice in R or Z?
if(sys.argv[2].upper() == 'Z'): IDXseek = data[1,:]
else: IDXseek = data[0,:]

if(sys.argv[2].upper() == 'Z'): plotXall = data[2,:]
else: plotXall = data[3,:]

DoTwo = 0
if( len(sys.argv) > 4 ): DoTwo = 1

# What coordinate to plot?
grid = int(sys.argv[3])
grid2 = grid
if(DoTwo): grid2 = int(sys.argv[4])
if(grid2>Params[1]-1 and sys.argv[2].upper() == 'Z'): grid2 = Params[1]-1
elif(grid2>Params[0]-1): grid2 = Params[0]-1

plotYall = data[idx,:]

RightIDX = [idxX for idxX in range(len(IDXseek)) if IDXseek[idxX]==grid]

plotX = [ plotXall[idxX] for idxX in RightIDX ]
plotY = [ plotYall[idxX] for idxX in RightIDX ]

VvalR = -1
VvalR2= -1
if(sys.argv[2].upper() == 'Z'): VvalR = VCont[grid]
if(sys.argv[2].upper() == 'Z'): VvalR2 = VCont[grid2]



if(DoTwo):
    RightIDX2 = [idxX for idxX in range(len(IDXseek)) if IDXseek[idxX]==grid2]
    plotX2 = [ plotXall[idxX] for idxX in RightIDX2 ]
    plotY2 = [ plotYall[idxX] for idxX in RightIDX2 ]


# if Vpot, get the contour value
if(idx==4): Vbdy = Params[-1]



# ---===---===---===---===---===---===---===---===---===
plt.subplot(1,1,1)
plt.subplots_adjust(left=0.2,bottom=0.15,wspace=0.0001,hspace=0.0001)

#plt.axes([0.16,0.13,0.79,0.82]) #left,bot,wid,height

maxY = max(plotY)
minY = min(plotY)
if(maxY==minY):
    maxY = 1.01*maxY
    minY = 0.99*minY

minY2  = minY
maxY2 = maxY

if(DoTwo):
    maxY2 = max(plotY2)
    minY2 = min(plotY2)
    if(maxY2==minY2):
        maxY2 = 1.01*maxY2
        minY2 = 0.99*minY2

deltaRange = max(maxY,maxY2) - min(minY,minY2)

minY = min( minY,minY2 ) - deltaRange/20.0
maxY = max( maxY,maxY2 ) + deltaRange/20.0

print plotX
print plotY
if(DoTwo): print plotY2

# Axis limits
plt.axis([0, max(plotX), minY,maxY])
#plt.axis([0.0, max(plotX), 0.8,1.1*maxY])

# Axis labels
plt.xlabel(r'Distance  (pc)') #TeX requires the r

if(idx==0): ytit = 'R Grid id'
if(idx==1): ytit = 'Z Grid id'
if(idx==2): ytit = 'R Distance (pc)'
if(idx==3): ytit = 'Z Distance (pc)'
if(idx==4): ytit = 'Gravitational Potential'
if(idx==5): ytit = 'A'
if(idx==6): ytit = 'Q'
if(idx==7): ytit = 'Density'
plt.title(ytit)

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
plt.plot(plotX,plotY,'k')
if(DoTwo): plt.plot(plotX2,plotY2,'g')


plt.plot([VvalR,VvalR],[minY,maxY],'k--',linewidth=2)
if(DoTwo): plt.plot([VvalR2,VvalR2],[minY,maxY],'g--',linewidth=2)

if(idx==4): plt.plot(plotX,len(plotX)*[Vbdy],'k--',linewidth=2)

if(comp): plt.plot(CompX,CompY,'--',color=dpurp)
#Text
#plt.text(0.05,0.06,r'${\cal M} = 4.47$',fontsize=26.0)

#plt.legend(loc=4,numpoints=1,ncol=2,labelspacing=0.1,columnspacing=0.9)


# Creates Plot
#plt.show()
plt.savefig('OneDPlot.pdf')
#plt.savefig('figure3.eps')

#Restores defaults
mpl.rcdefaults()
