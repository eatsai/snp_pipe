# Purpose:    Make the QC plots that can be found from the controls dashboard
## NB: May need a new module.
##     Instead of loading all of GTC.py, we only need the controls probes (which is fast))
# from optparse import OptionParser
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from BPM import *
from GTC import *
import numpy as np
import os


### Parse Inputs from Command Line

parser = ArgumentParser()

parser.add_argument("-B","--bpm",
                    type=str,
                    dest="bpm",
                    default='/Users/et85/testing/gtc/MEGA_Consortium_15063755_B1.csv',
                    help="csv-version of bpm file")

parser.add_argument("-O","--out",
                    type=str,
                    dest="out",
                    default='/Users/et85/workspace/zCall/custom/outdir/qcplot',
                    help="prefix of output")

parser.add_argument("-L","--gtclist",       # --gtclist a.gtc b.gtc c.gtc 
                    type=str,
                    dest='list',
                    nargs='+')

parser.add_argument("-D","--gtcdir",        # --gtcdir /path/to/dir
                    type=str,
                    dest='dir',
                    nargs=1)

parser.add_argument("-F","--gtcfilelist",       # --gtcfilelist samp_x100.filepath 
                    type=str,
                    dest='filelist')

args = parser.parse_args()

bpmPath = args.bpm
# print bpmPath

gtcPathList = []

# initialize variables for each plot
X_staining = []
Y_staining = []
X_extension = []
Y_extension = []
X_targetRemoval = []
Y_targetRemoval = []
X_hybridization = []
Y_hybridization = []
X_stringency = []
Y_stringency = []
X_nonSpecificBinding = []
Y_nonSpecificBinding = []
X_nonPolymorphic = []
Y_nonPolymorphic = []
X_restoration = []
Y_restoration = []
callRate = []
lrrSD = []

if args.list:
    print args.list
    gtcPathList = args.list
elif args.dir:
    allFiles = os.listdir(args.dir[0])
    for thisFile in allFiles:
        gtcPathList.append(os.path.join(args.dir[0], thisFile))
elif args.filelist:
    thisFH = open(args.filelist,"r")
    for thisLine in thisFH.readlines():
        thisLine = thisLine.rstrip()
        gtcPathList.append(thisLine)
else:
    print "ERROR! GTC file [list], [dir], or [filelist] must be specified."

# NB: list always gets used first, if someone specifies both, need to output warning
if sum(bool(x) for x in (args.list, args.dir, args.filelist)) != 1:
    print "WARNING: More than one input is specified. Will use the first specified in this order: GTC list, GTC dir, GTC filepath."

if args.out:
    outDirPath = args.out
else:
    print "Something went wrong, no args.path"

# Initialize flag for data QC check
flagSamp = np.zeros(len(gtcPathList),dtype=int)
flagStaining = np.zeros(len(gtcPathList),dtype=int)
flagExtension = np.zeros(len(gtcPathList),dtype=int)
flagTargetRemoval = np.zeros(len(gtcPathList),dtype=int)
flagHybridization = np.zeros(len(gtcPathList),dtype=int)
flagStringency = np.zeros(len(gtcPathList),dtype=int)
flagNonSpecificBinding = np.zeros(len(gtcPathList),dtype=int)
flagNonPolymorphic = np.zeros(len(gtcPathList),dtype=int)

for i in range(0,len(gtcPathList)):
    thisGTCPath = gtcPathList[i]
    thisGTC = GTC(thisGTCPath,'qcplot')     # TODO: the following could potentially be moved into GTC.py as a method
#     print thisGTC.rawXcontrol
#     print thisGTC.rawYcontrol
    X_staining.append(thisGTC.rawXcontrol[0:4*4])
    Y_staining.append(thisGTC.rawYcontrol[0:4*4])
    X_extension.append(thisGTC.rawXcontrol[4*4:4*8])
    Y_extension.append(thisGTC.rawYcontrol[4*4:4*8])
    X_targetRemoval.append(thisGTC.rawXcontrol[4*8:4*9])
    Y_targetRemoval.append(thisGTC.rawYcontrol[4*8:4*9])
    X_hybridization.append(thisGTC.rawXcontrol[4*9:4*12])
    Y_hybridization.append(thisGTC.rawYcontrol[4*9:4*12])
    X_stringency.append(thisGTC.rawXcontrol[4*12:4*14])
    Y_stringency.append(thisGTC.rawYcontrol[4*12:4*14])
    X_nonSpecificBinding.append(thisGTC.rawXcontrol[4*14:4*18]) # note: g/b are flipped on bpm.csv
    Y_nonSpecificBinding.append(thisGTC.rawYcontrol[4*14:4*18])
    X_nonPolymorphic.append(thisGTC.rawXcontrol[4*18:4*22])
    Y_nonPolymorphic.append(thisGTC.rawYcontrol[4*18:4*22])
    X_restoration.append(thisGTC.rawXcontrol[4*22:4*23])
    Y_restoration.append(thisGTC.rawYcontrol[4*22:4*23])
    callRate.append(thisGTC.callRate)
    lrrSD.append(thisGTC.lrrDev)
    
    if(thisGTC.callRate < 0.97):        # Automatically FAIL if Call Rate < 0.97
        print thisGTCPath + " " + str(thisGTC.callRate)


# with PdfPages('qcplots.pdf') as pdf:

# print outDirPath
# sys.exit()

with PdfPages(outDirPath + ".pdf") as pdf:
    
#     set_markeredgewidth(0)
    
    f, (red, grn) = plt.subplots(1,2,sharey=True)
    plt.suptitle('Staining',fontsize=16)
    thisMax = max([max(max(X_staining)),max(max(Y_staining))])  # shared axis, find max
    red.axis([0.5,len(X_staining)+0.5,0,thisMax+0.05*thisMax])
    grn.axis([0.5,len(X_staining)+0.5,0,thisMax+0.05*thisMax])
    for i in range(0,len(gtcPathList)):
        red.plot([i+1]*4, X_staining[i][0:4],'ro',[i+1]*4, X_staining[i][4:4*2],'mo',
                 [i+1]*4, X_staining[i][4*2:4*3],'go', [i+1]*4, X_staining[i][4*3:4*4],'bo',
                 markeredgewidth=0, markersize=1.5)
        grn.plot([i+1]*4, Y_staining[i][0:4],'ro',[i+1]*4, Y_staining[i][4:4*2],'mo',
                 [i+1]*4, Y_staining[i][4*2:4*3],'go', [i+1]*4, Y_staining[i][4*3:4*4],'bo',
                 markeredgewidth=0, markersize=2)
        
        # TODO: draw thresholds
        red.plot([0,len(callRate)+1], [1000, 1000], 'k--', lw=0.5, color='0.75')
        grn.plot([0,len(callRate)+1], [1000, 1000], 'k--', lw=0.5, color='0.75')
        
        
        
        # QC check: Staining
#         if(any(j > 1000 for j in X_staining[i][4:4*4])):
#             flagStaining[i] = 1
#         
#         if(any(j > 1000 for j in Y_staining[i][0:4*2]) | any(j > 1000 for j in Y_staining[i][4*3:4*4])):
#             flagStaining[i] = 1

        # Check if X_stain > Y_stain (red), X_stain > Y_stain (green)
        if(any(j > min(X_staining[i][0:4]) for j in X_staining[i][4:4*4])):
            flagSamp[i] = 1
        
        if(any(j > min(Y_staining[i][4*2:4*3]) for j in Y_staining[i][4:4*2]) | 
           any(j > min(Y_staining[i][4*2:4*3]) for j in Y_staining[i][4*3:4*4])):
            flagSamp[i] = 1
        
    
    red.set_title('Red')
    grn.set_title('Green')
    pdf.savefig()   # saves the current figure into a pdf page
    plt.close()
        
    
    f, (red, grn) = plt.subplots(1,2,sharey=True)
    plt.suptitle('Extension',fontsize=16)
    thisMax = max([max(max(X_extension)),max(max(Y_extension))])  # shared axis, find max
    red.axis([0.5,len(X_extension)+0.5,0,thisMax+0.05*thisMax])
    grn.axis([0.5,len(X_extension)+0.5,0,thisMax+0.05*thisMax])
    for i in range(0,len(gtcPathList)):
        red.plot([i+1]*4, X_extension[i][0:4],'ro',[i+1]*4, X_extension[i][4:4*2],'mo',
                 [i+1]*4, X_extension[i][4*2:4*3],'go', [i+1]*4, X_extension[i][4*3:4*4],'bo',
                 markeredgewidth=0, markersize=2)
        grn.plot([i+1]*4, Y_extension[i][0:4],'ro',[i+1]*4, Y_extension[i][4:4*2],'mo',
                 [i+1]*4, Y_extension[i][4*2:4*3],'go', [i+1]*4, Y_extension[i][4*3:4*4],'bo',
                 markeredgewidth=0, markersize=2)
        
        # QC check: Extension
#         if(any(j > 4000 for j in X_extension[i][4*2:4*3]) | any(j > 4000 for j in X_extension[i][4*3:4*4])):
#             flagExtension[i] = 1
#         
#         if(any(j > 1000 for j in Y_extension[i][0:4]) | any(j > 1000 for j in Y_extension[i][4:4*2])):
#             flagExtension[i] = 1
        
        # Check to make sure that signal > noise
        if(any(j > max(X_extension[i][0:4*2]) for j in X_extension[i][4*2:4*4])):
            flagSamp[i] = 1
        
        if(any(j > max(Y_extension[i][4*2:4*4]) for j in Y_extension[i][0:4*2])):
            flagSamp[i] = 1
        
        # Draw thresholds
        red.plot([0,len(callRate)+1], [4000, 4000], 'k--', lw=0.5, color='0.75')
        grn.plot([0,len(callRate)+1], [1000, 1000], 'k--', lw=0.5, color='0.75')
        
       
        
    red.set_title('Red')
    grn.set_title('Green')
    pdf.savefig()   # saves the current figure into a pdf page
    plt.close()
    
    f, (red, grn) = plt.subplots(1,2,sharey=True)
    plt.suptitle('Target Removal',fontsize=16)
    thisMax = max([max(max(X_targetRemoval)),max(max(Y_targetRemoval))])  # shared axis, find max
    red.axis([0.5,len(X_targetRemoval)+0.5,0,thisMax+0.05*thisMax])
    grn.axis([0.5,len(X_targetRemoval)+0.5,0,thisMax+0.05*thisMax])
    for i in range(0,len(gtcPathList)):
        red.plot([i+1]*4, X_targetRemoval[i][0:4],'bo',
                 markeredgewidth=0, markersize=2)
        grn.plot([i+1]*4, Y_targetRemoval[i][0:4],'bo',
                 markeredgewidth=0, markersize=2)

        # QC check: TargetRemoval
        if(any(j > 5000 for j in X_targetRemoval[i][0:4])):
            flagTargetRemoval[i] = 1
        
        if(any(j > 1000 for j in Y_targetRemoval[i][0:4])):
            flagTargetRemoval[i] = 1
        
        # Draw thresholds
        red.plot([0,len(callRate)+1], [5000, 5000], 'k--', lw=0.5, color='0.75')
        grn.plot([0,len(callRate)+1], [1000, 1000], 'k--', lw=0.5, color='0.75')
        


    red.set_title('Red')
    grn.set_title('Green')
    pdf.savefig()   # saves the current figure into a pdf page
    plt.close()
    
    f, (red, grn) = plt.subplots(1,2,sharey=True)
    plt.suptitle('Hybridization',fontsize=16)
    thisMax = max([max(max(X_hybridization)),max(max(Y_hybridization))])  # shared axis, find max
    red.axis([0.5,len(X_hybridization)+0.5,0,thisMax+0.05*thisMax])
    grn.axis([0.5,len(X_hybridization)+0.5,0,thisMax+0.05*thisMax])
    for i in range(0,len(gtcPathList)):
        red.plot([i+1]*4, X_hybridization[i][0:4],'go',[i+1]*4, X_hybridization[i][4:4*2],'bo',
                 [i+1]*4, X_hybridization[i][4*2:4*3],'ko',
                 markeredgewidth=0, markersize=2)
        grn.plot([i+1]*4, Y_hybridization[i][0:4],'go',[i+1]*4, Y_hybridization[i][4:4*2],'bo',
                 [i+1]*4, Y_hybridization[i][4*2:4*3],'ko',
                 markeredgewidth=0, markersize=2)
        
        # QC check: Hybridization
#         if(any(j > 5000 for j in X_hybridization[i][0:4*3])):
#             flagHybridization[i] = 1

        # Make sure high > medium > low for both red/green channels         
        if(max(X_hybridization[i][4:4*2]) > min(X_hybridization[i][0:4]) | 
           max(X_hybridization[i][4*2:4*3]) > min(X_hybridization[i][4:4*2])):
            flagSamp[i] = 1
       
        if(max(Y_hybridization[i][4:4*2]) > min(Y_hybridization[i][0:4]) | 
           max(Y_hybridization[i][4*2:4*3]) > min(Y_hybridization[i][4:4*2])):
            flagSamp[i] = 1
        
        # Draw thresholds
        red.plot([0,len(callRate)+1], [5000, 5000], 'k--', lw=0.5, color='0.75')     
        
        
    red.set_title('Red')
    grn.set_title('Green')
    pdf.savefig()   # saves the current figure into a pdf page
    plt.close()
    
    f, (red, grn) = plt.subplots(1,2,sharey=True)
    plt.suptitle('Stringency',fontsize=16)
    thisMax = max([max(max(X_stringency)),max(max(Y_stringency))])  # shared axis, find max
    red.axis([0.5,len(X_stringency)+0.5,0,thisMax+0.05*thisMax])
    grn.axis([0.5,len(X_stringency)+0.5,0,thisMax+0.05*thisMax])
    for i in range(0,len(gtcPathList)):
        red.plot([i+1]*4, X_stringency[i][0:4],'ro',[i+1]*4, X_stringency[i][4:4*2],'mo',
                 markeredgewidth=0, markersize=2)
        grn.plot([i+1]*4, Y_stringency[i][0:4],'ro',[i+1]*4, Y_stringency[i][4:4*2],'mo',
                 markeredgewidth=0, markersize=2)
        
        # Check to make sure signal > noise
        if(any(j > max(X_stringency[i][0:4]) for j in X_stringency[i][4:4*2])):
            flagStringency[i] = 1
        
        # Control for background
#         if(any(j > 1000 for j in Y_stringency[i][4:4*2])):
#             flagSamp[i] = 1
        
        # Draw thresholds
        grn.plot([0,len(callRate)+1], [1000, 1000], 'k--', lw=0.5, color='0.75')     
        
    red.set_title('Red')
    grn.set_title('Green')
    pdf.savefig()   # saves the current figure into a pdf page
    plt.close()
    
    f, (red, grn) = plt.subplots(1,2,sharey=True)
    plt.suptitle('Non-specific Binding',fontsize=16)
    thisMax = max([max(max(X_nonSpecificBinding)),max(max(Y_nonSpecificBinding))])  # shared axis, find max
    red.axis([0.5,len(X_nonSpecificBinding)+0.5,0,thisMax+0.05*thisMax])
    grn.axis([0.5,len(X_nonSpecificBinding)+0.5,0,thisMax+0.05*thisMax])
    for i in range(0,len(gtcPathList)):
        red.plot([i+1]*4, X_nonSpecificBinding[i][0:4],'ro',[i+1]*4, X_nonSpecificBinding[i][4:4*2],'mo',
                 [i+1]*4, X_nonSpecificBinding[i][4*3:4*4],'go', [i+1]*4, X_nonSpecificBinding[i][4*2:4*3],'bo',
                 markeredgewidth=0, markersize=2)
        grn.plot([i+1]*4, Y_nonSpecificBinding[i][0:4],'ro',[i+1]*4, Y_nonSpecificBinding[i][4:4*2],'mo',
                 [i+1]*4, Y_nonSpecificBinding[i][4*3:4*4],'go', [i+1]*4, Y_nonSpecificBinding[i][4*2:4*3],'bo',
                 markeredgewidth=0, markersize=2)
        
        # Control for background noise
        if(any(j > 1500 for j in X_nonSpecificBinding[i][0:4*4])):
            flagNonSpecificBinding[i] = 1
        
        if(any(j > 1500 for j in Y_nonSpecificBinding[i][0:4*4])):
            flagNonSpecificBinding[i] = 1
        
        # TODO: draw thresholds
        red.plot([0,len(callRate)+1], [1000, 1000], 'k--', lw=0.5, color='0.75')
        grn.plot([0,len(callRate)+1], [1000, 1000], 'k--', lw=0.5, color='0.75')
        
    red.set_title('Red')
    grn.set_title('Green')
    pdf.savefig()   # saves the current figure into a pdf page
    plt.close()
    
    
    #########
    ### Boxplot
    #########
    f, (red, grn) = plt.subplots(1,2,sharey=True)
    plt.suptitle('Non-specific Binding Boxplot',fontsize=16)
    mean_x = np.mean(X_nonSpecificBinding)
    sd_x = np.sqrt(np.var(X_nonSpecificBinding)*len(X_nonSpecificBinding)/(len(X_nonSpecificBinding)-1))
    mean_y = np.mean(Y_nonSpecificBinding)
    sd_y = np.sqrt(np.var(Y_nonSpecificBinding)*len(Y_nonSpecificBinding)/(len(Y_nonSpecificBinding)-1))
    red.boxplot(X_nonSpecificBinding)
    grn.boxplot(Y_nonSpecificBinding)
    
    # Draw thresholds
    red.plot([0,len(callRate)+1], [1500, 1500], 'k--', lw=0.5, color='0.75')
    grn.plot([0,len(callRate)+1], [1500, 1500], 'k--', lw=0.5, color='0.75')
    
    red.set_title('Red')
    grn.set_title('Green')
    pdf.savefig()   # saves the current figure into a pdf page
    plt.close()
    
    f, (red, grn) = plt.subplots(1,2,sharey=True)
    plt.suptitle('Non-polymorphic',fontsize=16)
    thisMax = max([max(max(X_nonPolymorphic)),max(max(Y_nonPolymorphic))])  # shared axis, find max
    red.axis([0.5,len(X_nonPolymorphic)+0.5,0,thisMax+0.05*thisMax])
    grn.axis([0.5,len(X_nonPolymorphic)+0.5,0,thisMax+0.05*thisMax])
    for i in range(0,len(gtcPathList)):
        red.plot([i+1]*4, X_nonPolymorphic[i][0:4],'ro',[i+1]*4, X_nonPolymorphic[i][4:4*2],'mo',
                 [i+1]*4, X_nonPolymorphic[i][4*2:4*3],'go', [i+1]*4, X_nonPolymorphic[i][4*3:4*4],'bo',
                 markeredgewidth=0, markersize=2)
        grn.plot([i+1]*4, Y_nonPolymorphic[i][0:4],'ro',[i+1]*4, Y_nonPolymorphic[i][4:4*2],'mo',
                 [i+1]*4, Y_nonPolymorphic[i][4*2:4*3],'go', [i+1]*4, Y_nonPolymorphic[i][4*3:4*4],'bo',
                 markeredgewidth=0, markersize=2)
        
        # Check to make sure signal > noise
        if(any(j > max(X_nonPolymorphic[i][0:4*2]) for j in X_nonPolymorphic[i][4*2:4*4])):
            flagNonPolymorphic[i] = 1
        
        # Check to make sure signal > noise
        if(any(j > max(Y_nonPolymorphic[i][4*2:4*4]) for j in Y_nonPolymorphic[i][0:4*2])):
            flagNonPolymorphic[i] = 1
        
        
        # control for background
#         if(any(j > 1000 for j in X_nonPolymorphic[i][4*2:4*4])):
#             flagNonPolymorphic[i] = 1
#         
#         if(any(j > 1000 for j in Y_nonPolymorphic[i][0:4*2])):
#             flagNonPolymorphic[i] = 1
        
        # Draw thresholds
        red.plot([0,len(callRate)+1], [1000, 1000], 'k--', lw=0.5, color='0.75')
        grn.plot([0,len(callRate)+1], [1000, 1000], 'k--', lw=0.5, color='0.75')
        
        
    red.set_title('Red')
    grn.set_title('Green')
    pdf.savefig()   # saves the current figure into a pdf page
    plt.close()
    
    f, (red, grn) = plt.subplots(1,2,sharey=True)
    plt.suptitle('Restoration (FFPE)',fontsize=16)
    thisMax = max([max(max(X_restoration)),max(max(Y_restoration))])  # shared axis, find max
    red.axis([0.5,len(X_restoration)+0.5,0,thisMax+0.05*thisMax])
    grn.axis([0.5,len(X_restoration)+0.5,0,thisMax+0.05*thisMax])
    for i in range(0,len(gtcPathList)):
        red.plot([i+1]*4, X_restoration[i][0:4],'go',
                 markeredgewidth=0, markersize=2)
        grn.plot([i+1]*4, Y_restoration[i][0:4],'go',
                 markeredgewidth=0, markersize=2)
    red.set_title('Red')
    grn.set_title('Green')
    pdf.savefig()   # saves the current figure into a pdf page
    plt.close()
    
    plt.figure()
    plt.suptitle('Call Rate',fontsize=16)
    callRate_U995 = [x for x in callRate if (x >= 0.995)]
    indCallRate_U995 = [i+1 for i,x in enumerate(callRate) if (x >= 0.995)]
    callRate_L995 = [x for x in callRate if (0.99 <= x < 0.995)]
    indCallRate_L995 = [i+1 for i,x in enumerate(callRate) if (0.99 <= x < 0.995)]
    callRate_L99 = [x for x in callRate if (x < 0.99)]
    indCallRate_L99 = [i+1 for i,x in enumerate(callRate) if (x < 0.99)]
    thisMin = min(callRate)
    plt.axis([0,len(callRate)+1,thisMin - thisMin*0.05,1])
    plt.plot(indCallRate_U995, callRate_U995,'ko',
         indCallRate_L995, callRate_L995,'yo',
         indCallRate_L99, callRate_L99,'ro', markeredgewidth=0, markersize=2)
#     # Ideal thresholds
    plt.plot([0,len(callRate)+1], [0.995, 0.995], 'y--', lw=0.5)
    plt.plot([0,len(callRate)+1], [0.99, 0.99], 'r--', lw=0.5)
    # Thresholds with our dataset
#     plt.plot([0,len(callRate)+1], [0.99, 0.99], 'y--', lw=0.5)
#     plt.plot([0,len(callRate)+1], [0.97, 0.97], 'r--', lw=0.5)
    pdf.savefig()   # saves the current figure into a pdf page
    plt.close()
    
    plt.figure()
    plt.suptitle('SD(LogR Ratio)',fontsize=16)
    bins = [0.05,.1,0.15,0.175,0.2,0.225,0.25,0.3,0.325,0.35,0.4,0.45,0.5,1]
    # the histogram of the data with histtype='step'
#     n, bins, patches = plt.hist(lrrSD, bins, normed=1, histtype='bar', rwidth=0.8)
    n, bins, patches = plt.hist(lrrSD, bins, normed=0, histtype='bar', rwidth=0.8)
    pdf.savefig()   # saves the current figure into a pdf page
    plt.close()
    
    # Check to make sure sample passes set flags
    thisInd = [index for index,value in enumerate(callRate) if value < 0.99]
#     print callRate
    for i in thisInd:
        flagSamp[i] = 1
    

outFH = open(outDirPath + ".txt", "w")

outFH.write('\t'.join(['SentrixID','QC_Flag', 'CallRate', 'SD(LRR)', 
                       'StainingFlag', 'ExtensionFlag', 'TargetRemovalFlag', 'HybridizationFlag', 'StringencyFlag', 'NSBindingFlag', 'NonPolymorphicFlag',
                       'X_Staining-R', 'X_Staining-P', 'X_Staining-G', 'X_Staining-B', 'Y_Staining-R', 'Y_Staining-P', 'Y_Staining-G', 'Y_Staining-B', 'X_Extension-R', 'X_Extension-P', 'X_Extension-G', 'X_Extension-B', 'Y_Extension-R', 'Y_Extension-P', 'Y_Extension-G', 'Y_Extension-B', 'X_TargetRemoval', 'Y_TargetRemoval', 'X_Hybridization_High', 'X_Hybridization_Medium', 'X_Hybridization_Low', 'Y_Hybridization_High', 'Y_Hybridization_Medium', 'Y_Hybridization_Low', 'X_Stringency_PM', 'X_Stringency_MM', 'Y_Stringency_PM', 'Y_Stringency_MM', 'X_NonSpecificBinding-R', 'X_NonSpecificBinding-P', 'X_NonSpecificBinding-B', 'X_NonSpecificBinding-G', 'Y_NonSpecificBinding-R', 'Y_NonSpecificBinding-P', 'Y_NonSpecificBinding-B', 'Y_NonSpecificBinding-G', 'X_NonPolymorphic-R', 'X_NonPolymorphic-P', 'X_NonPolymorphic-G', 'X_NonPolymorphic-B', 'Y_NonPolymorphic-R', 'Y_NonPolymorphic-P', 'Y_NonPolymorphic-G', 'Y_NonPolymorphic-B', 'X_Restoration', 'Y_Restoration']) + '\n')

for i in range(0,len(gtcPathList)):
    thisSentrixID = os.path.basename(gtcPathList[i]).split(".")[0]
    outFH.write(thisSentrixID + '\t')
    outFH.write(str(flagSamp[i]) + '\t')
    outFH.write(str(callRate[i]) + '\t')
    outFH.write(str(lrrSD[i]) + '\t')
    outFH.write(str(flagStaining[i]) + '\t')
    outFH.write(str(flagExtension[i]) + '\t')
    outFH.write(str(flagTargetRemoval[i]) + '\t')
    outFH.write(str(flagHybridization[i]) + '\t')
    outFH.write(str(flagStringency[i]) + '\t')
    outFH.write(str(flagNonSpecificBinding[i]) + '\t')
    outFH.write(str(flagNonPolymorphic[i]) + '\t')
    outFH.write('\t'.join(str(x) for x in X_staining[i][0::4]) + '\t') # taking every fourth value
    outFH.write('\t'.join(str(x) for x in Y_staining[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in X_extension[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in Y_extension[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in X_targetRemoval[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in Y_targetRemoval[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in X_hybridization[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in Y_hybridization[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in X_stringency[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in Y_stringency[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in X_nonSpecificBinding[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in Y_nonSpecificBinding[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in X_nonPolymorphic[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in Y_nonPolymorphic[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in X_restoration[i][0::4]) + '\t')
    outFH.write('\t'.join(str(x) for x in Y_restoration[i][0::4]))
    outFH.write('\n')
    
    


    ##TODO: trying to figure out how to place legend
#     plt.axis([0, 4, 0, 30000])
#     plt.plot([1]*4, X_staining[0][0:4],'ro',[1]*4, X_staining[0][4:4*2],'mo',
#              [1]*4, X_staining[0][4*2:4*3],'go', [1]*4, X_staining[0][4*3:4*4],'bo')
#     plt.title('Staining - Red')
#     line_up = plt.plot([1,2,3], label='Line 2')
#     line_down, = plt.plot([3,2,1], label='Line 1')
#     plt.legend(handles=[line_up, line_down])
#     plt.legend([line_up,line_down],['Line Up','Line Down'])
#     plt.legend(bbog_to_anchor=(1.05,1),loc=2,borderaxespad=0.)
#     pdf.savefig()   # saves the current figure into a pdf page
#     plt.close()
    

