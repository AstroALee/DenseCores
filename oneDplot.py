# Imports
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import scipy as spi
import sys

import AnalyticSolutions as AS

# Colors
dmagenta = '#DB1672'
dgreen = '#259B41'
dpurp  = '#BB46C9'

# INPUTS
# -----------
#            0           1    2    3     4--
# python newoneDplot.py idx units dir [ rows ]
# idx = parameter to plot
# units = convert to dimensional?
# dir = slice direction
# [rows] list of rows to plot in slice (requires at least one)


# Checks
if(len(sys.argv)<5):
    print("Need more arguments")
    print("python newoneDplot.py idx units(y or n) dir(Z or R) [ rows ]")
    exit(1)

# Are we plotting the analytic solutions for comparison?
PlotAS = 1
PlotEQ = 1 # from equations vs code?
# Are we plotting the contour locations?
PlotCT = 1



# Reads in raw data as column arrays
data = genfromtxt('Output.out',delimiter=",",unpack=True,skip_header=2)
# Read in contour array
flen = len(data[0,:])
VCont = genfromtxt('Output.out',delimiter=",",unpack=True,skip_header=1,skip_footer=flen)
# M, N, zL (pc), rRatio (entire box), mExcess (sol), beta, n, Rbdy, Lambda, Pc2Code, Sol2Code
Params = genfromtxt('Output.out',delimiter=",",unpack=True,skip_footer=flen+1)

M = int(Params[0])
N = int(Params[1])
zL = float(Params[2])
rL = float(Params[3])*zL
DeltaR = rL/float(M)
DeltaZ = zL/float(N)
beta = float(Params[5])
nCyl = float(Params[6])
Rbdy = float(Params[7]) # Rho at edge of filament
dLam = float(Params[8])

Pc2Code = float(Params[-2])
Sol2Code = float(Params[-1])

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
if(idxS=='9' or idxS.lower()=='j'): idx = 9

if(idx==-1):
    print "ERROR with idx"
    exit(1)

# Are we converting to dimensional units?
Units = 0
if( str(sys.argv[2]).lower()=='y'): Units = 1

# Get the distance array depending on what slice we are doing
Slice = sys.argv[3].upper()
if(Slice=='Z'): XCord = sort(list(set(data[2,:]))) # Slicing in Z, list of R values
else: XCord = sort(list(set(data[3,:]))) # Slicing in R, list of Z values

# Which rows or columns are we plotting
PlotThese = sort(list(set(sys.argv[4:]))) # gets rid of duplicates
NumLines = len(PlotThese)
PlotThese = [ int(PlotThese[i]) for i in range(NumLines) ]
for i in range(0,NumLines):
    locIdx = PlotThese[i]
    if(Slice=='Z' and locIdx>N-1): locIdx = N-1
    elif(locIdx>M-1): locIdx = M-1
    PlotThese[i] = locIdx
PlotThese = sort(list(set(PlotThese))) # gets rid of duplicates
NumLines = len(PlotThese)
PlotThese = [ str(PlotThese[i]) for i in range(NumLines) ]


# Get the indices we'll be searching through to find the values in 'PlotThese'
if(Slice=='Z'): IDXseek = data[1,:] # Slice in Z, list of Z indices
else: IDXseek = data[0,:] # Slice in R, list of R indices

# The indices we will need for each line we plot
Indices = []
for i in range(0,NumLines):
    locIdx = [idxX for idxX in range(len(IDXseek)) if IDXseek[idxX]==int(PlotThese[i])]
    Indices.append(locIdx)

# Create the data we will plot
YCordAll = []
if(idx < 8):   # data from output file (non-derived quantities)
    YCordAll = data[idx,:]
elif(idx==8):  # Phi
    YCordAll = multiply( data[2,:] , data[5,:] )
elif(idx==9): # computes r*exp(-V)*dQ/dPhi (propto j)
    Qval = data[6,:]
    Phival = multiply(data[2,:],data[5,:])
    dQdPhiMat = zeros(len(data[2,:]))
    for i in range(0,len(data[2,:])): # dQ= Q(i+1) + Q(j+1) - Q(i-1) - Q(j-1)
        dQval = 0
        dPval = 0
        if( i%N > 0  and  i%N < N-1 ):  #not on left or right side where there's no r derivatives by symmetry
            dQval = dQval + Qval[i+1] - Qval[i-1]
            dPval = dPval + Phival[i+1] - Phival[i-1]
        if( i/N > 0 and i/N < M-1 ): # not on top or bottom, where no z derivatives
            dQval = dQval + Qval[i+N] - Qval[i-N]
            dPval = dPval + Phival[i+N] - Phival[i-N]
        if(dPval != 0):
            dQdPhiMat[i] = dQval/dPval
    #YCordAll = multiply( multiply(data[2,:],dQdPhiMat) , exp(-data[4,:]) )
    YCordAll = multiply( multiply(data[2,:],data[8,:]) , exp(-data[4,:]) )
    #YCordAll = data[8,:]


# Now pull out the data we will plot
YCord = []
for i in range(0,NumLines):
    locAry = [ YCordAll[j] for j in Indices[i] ]
    YCord.append(locAry)
print YCord


# If plotting analytic solutions, construct those before unit conversions happen
# data[7,N-1] = rho(r=0) at top
if(PlotEQ):
    print(beta)
    if(idx==4):
        Yanal = [ AS.Veqn(i,VCont[0],beta,data[7,N-1]) for i in XCord]
    elif(idx==5):
        Yanal = [ AS.Aeqn(i,VCont[0],beta,data[7,N-1],Rbdy,nCyl) for i in XCord]
    elif(idx==6):
        Yanal = [ AS.Qeqn(i,VCont[0],beta,data[7,N-1]) for i in XCord]
    elif(idx==7):
        Yanal = [ AS.RHOeqn(i,VCont[0],beta,data[7,N-1]) for i in XCord]
    elif(idx==8):
        Yanal = [ AS.PHIeqn(i,VCont[0],beta,data[7,N-1],Rbdy,nCyl) for i in XCord]
    elif(idx==9):
        Yanal = [ AS.DQDZeqn(i,VCont[0],beta,data[7,N-1],Rbdy,nCyl) for i in XCord]



# ---===---===---===---===---===---===---===---===---===
plt.subplot(1,1,1)
plt.subplots_adjust(left=0.2,bottom=0.15,wspace=0.0001,hspace=0.0001)

XCord = array(XCord)
YCord = array(YCord)

# Units
plotConvert = 1.0
disConvert = 1.0

if(Units==1):
    disConvert = 1.0/Pc2Code

# Converts units if needed
VCont = VCont*disConvert
XCord = XCord*disConvert
YCord = YCord*plotConvert

# Max and Min values
minY =  YCord.min()
maxY =  YCord.max()

if(PlotEQ and PlotAS):
    minY = min(minY,min(Yanal))
    maxY = max(maxY,max(Yanal))

if(maxY==minY):
    maxY = 1.01*maxY
    minY = 0.99*minY
deltaRange = maxY-minY

# Give a little padding above and below
if(minY==0): minY = 0
else: minY = minY - deltaRange/20.0
maxY = maxY + deltaRange/20.0

# Axis limits
plt.axis([0, XCord.max(), minY,maxY])
#plt.axis([0, XCord.max(), minY,1])

# Axis labels
if(Units): plt.xlabel(r'Distance  (pc)') #TeX requires the r
else: plt.xlabel(r'Distance  (ndim)') #TeX requires the r

if(idx==0): ytit = 'R Grid id'
if(idx==1): ytit = 'Z Grid id'
if(idx==2): ytit = 'R Distance (pc)'
if(idx==3): ytit = 'Z Distance (pc)'
if(idx==4): ytit = 'Gravitational Potential'
if(idx==5): ytit = 'A'
if(idx==6): ytit = 'Q'
if(idx==7): ytit = 'Density'
if(idx==8): ytit = 'Phi'
if(idx==9): ytit = 'J (Current)'

if(Units): ytit = ytit + '  (dim)'
else: ytit = ytit + '  (ndim)'
plt.title(ytit)

#Moves the x-axis down a little
plt.gca().xaxis.labelpad=9 #5 is default

#Makes ticks thicker/longer
ax = plt.gca()
ax.tick_params(which='both',width=2)
ax.tick_params(which='major',size=12)
ax.tick_params(which='minor',size=7)

# moves numbers out a little
ax.tick_params(axis='y',pad=9)
ax.tick_params(axis='x',pad=9)

for i in range(NumLines):
    alp = 1.0 - float(i)/(1.0+float(NumLines)) # darker lines are closer to midplane
    plt.plot(XCord,YCord[i],'k',alpha=alp)
    #plt.plot(log10(XCord),log10(YCord[i]),'k')

# If we're slicing in Z, plot the filament boundary locations
if(Slice=='Z' and PlotCT):
    for i in range(NumLines):
        alp = 1.0 - float(i)/(1.0+float(NumLines))
        Vloc = VCont[int(PlotThese[i])]
        plt.plot([Vloc,Vloc],[minY,maxY],'k--',alpha=alp,linewidth=2)

#Anayltic solution
if(PlotAS and Slice=='Z'):
    linewidth=3
    clr = 'red'
    if(PlotEQ):
        if(idx>=4 and idx<=9): plt.plot(XCord,Yanal,'--',color=clr,linewidth=linewidth)
    else:
        if(idx==4): plt.plot(AS.R[beta],AS.V[beta],'--',color=clr,linewidth=linewidth)
        if(idx==5): plt.plot(AS.R[beta],AS.A[beta],'--',color=clr,linewidth=linewidth)
        if(idx==6): plt.plot(AS.R[beta],AS.Q[beta],'--',color=clr,linewidth=linewidth)

#plt.axis([0, 0.5,0,1])

# Creates Plot
#plt.show()
plt.savefig('OneDPlot.pdf')

#Restores defaults
mpl.rcdefaults()
