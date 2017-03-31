#!/usr/bin/env python
import gzip
import os
import os.path
import argparse
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

import numpy as np
import pandas as pd

#1. GET READS
def getReads(fileDescriptor=None):
    '''
    Reads a single line from the file descriptor
    Returns an empty string when the end of the files has been reached
    fileDescriptor - File Descriptor to get a line
    '''
    return fileDescriptor.readline()

#2. GET WINDOW
def getWindow(lineRead=None):
    '''
    Split lineRead in fields, chromosome Name, start, end, methylation value,windowMatched
    return list of fields
    lineRead - array of window fields
    '''
    tmpList = lineRead.rstrip().split('\t')
    if len(tmpList) == 4:
        return [tmpList[0],int(tmpList[1]),int(tmpList[2]),float(tmpList[3]),False]
    return tmpList

#3. GET METHYLATION RANGE
def getMethRange(methValue=0):
    '''
    Get methylation ranges:   |--range 1--|--range 2--|--range 3--|
                              0          0.3         0.7          1
    methValue - Value to extract the range
    return range
    '''
    if methValue <= 0.3:
        return 1
    elif methValue > 0.3 and methValue <= 0.7:
        return 2
    else:
        return 3

#4.COMPARE VALUES
def compareMethValues(methValueFirst=0,methValueSecond=0):
    '''
    Compare methylation values.
    Returns True if values are considered as equal, otherwise returns false
    methValueFirst - First Methylation Value
    methValueSecond - Second Methylation Value
    '''
    diff = methValueFirst - methValueSecond
    if abs(diff) <= 0.1:
        return True
    else:
        if getMethRange(methValueFirst) == getMethRange(methValueSecond):
            return True

    return False

#5. COMPARE WINDOWS
def compareWindows(arrayWindowFirst=None,arrayWindowSecond=None,previousContig="",relaxedCompare=False):
    '''
    Compare two given windows by each location. Return 0 if both are equal, 1 if the firstone is smaller otherwise 2.
    If one of the contigs are different from the previous one then the one which is like the previous contig name is considered as smaller.
    arrayWindowFirst - Array of fields for the first window.
    arrayWindowSecond - Array of fields for the second window.
    previousContig - Previous contig name
    relaxedCompare - Compares if a position in a second windows is within first window ranges
    return 0 windows have the same location, 1 First window is smaller than secondone, 2 Second window is smaller than the firstone, 11 First includes second, 12 Second Includes First
    '''
    if arrayWindowFirst[0] == arrayWindowSecond[0]:
        if arrayWindowFirst[1] == arrayWindowSecond[1] and arrayWindowFirst[2] == arrayWindowSecond[2]:
            arrayWindowFirst[4] = True
            arrayWindowSecond[4] = True
            return 0
        elif relaxedCompare == True:
            if arrayWindowFirst[1] == (arrayWindowSecond[1] + 1) and arrayWindowFirst[2] == (arrayWindowSecond[2] + 1):
                arrayWindowFirst[4] = True
                arrayWindowSecond[4] = True
                return 0 
            elif ( (arrayWindowFirst[1] <= arrayWindowSecond[1] and arrayWindowSecond[2] <= arrayWindowFirst[2]) or
                  (arrayWindowFirst[1] <= (arrayWindowSecond[1] + 1) and (arrayWindowSecond[2] + 1) <= arrayWindowFirst[2]) ):
                arrayWindowFirst[4] = True 
                return 11
            elif ( (arrayWindowFirst[1] >= arrayWindowSecond[1] and arrayWindowSecond[2]  >= arrayWindowFirst[2]) or
                 (arrayWindowFirst[1] >= (arrayWindowSecond[1] + 1) and (arrayWindowSecond[2] + 1) >= arrayWindowFirst[2]) ): 
                arrayWindowSecond[4] = True
                return 12
        
        if arrayWindowFirst[1] < arrayWindowSecond[1] and arrayWindowFirst[2] < arrayWindowSecond[2]:
            return 1
        elif arrayWindowFirst[1] > arrayWindowSecond[1] and arrayWindowFirst[2] > arrayWindowSecond[2]:
            return 2
    elif arrayWindowFirst[0] == previousContig:
        return 1
    elif arrayWindowSecond[0] == previousContig:
        return 2

#6. COMPUTE COUNTERS
def computeCounters(methValue=0.0,counterUnmeth=0,counterUndef=0,counterMeth=0):
    '''
    Compute counters according to window methylation value
    methValue - Window methylation value
    counterUnmeth - Counter Unmethylation value
    counterUndef - Counter Undefined methylation value
    counterMeth - Counter Methylation value
    '''
    methRange = getMethRange(methValue)
    if methRange == 1:
        counterUnmeth = counterUnmeth + 1
    elif methRange == 2:
        counterUndef = counterUndef + 1
    elif methRange == 3:
        counterMeth = counterMeth + 1

    return counterUnmeth,counterUndef,counterMeth

#7.PLOT STATS
def plotsStats(statsUniqueFirst,statsUniqueSecond,statsShared,different,label_first,label_second,outpng):
    '''
    Plots a bar chart with the comarison stats
    statsUniqueFirst - Tuple of unique stats for file one
    statsUniqueSecond - Tuple of unique stats for file two
    statsShared - Tuple of shared stats 
    different - number of different stats
    label_first - Label given for the first file
    label_second - Label given for the second file
    outpng - png output file
    '''
    #Plot based on: http://chrisalbon.com/python/matplotlib_grouped_bar_plot.html
    

    raw_data = {'concept': ['Unique %s' %(label_first), 'Unique %s' %(label_second), 'Equal Meth', 'Diff Meth'],
                'unMeth': [statsUniqueFirst[0],statsUniqueSecond[0],statsShared[0],0],
                'unDef': [statsUniqueFirst[1],statsUniqueSecond[1],statsShared[1],0],
                'meth': [statsUniqueFirst[2],statsUniqueSecond[2],statsShared[2],0],
                'different': [0,0,0,different]}

    df = pd.DataFrame(raw_data, columns = ['concept', 'unMeth', 'unDef', 'meth','different'])
    
    # Data Frame
    #  ===+=========+========+========+=======+=======+
    #   0 |concept  |unMeth  |unDef   |meth   |diff   |
    #  ===+=========+========+========+=======+=======+
    #   1 |uniq 1st |unMeth  |unDef   |meth   |diff   |
    #  ---+---------+--------+--------+-------+-------+
    #   2 |uniq 2nd |unMeth  |unDef   |meth   |diff   |
    #  ---+---------+--------+--------+-------+-------+
    #   3 |eq meth  |unMeth  |unDef   |meth   |diff   |
    #  ---+---------+--------+--------+-------+-------+
    #   4 |diff meth|unMeth  |unDef   |meth   |diff   |
    #  ---+---------+--------+--------+-------+-------+

    #Setting the positions and width for the bars
    pos = list(range(len(df['unMeth'])))
    width = 0.25

    # Plotting the bars
    fig, ax = plt.subplots(figsize=(10,5))

    #Create a bar with unMeth data in position pos
    plt.bar(pos,df['unMeth'],width,alpha=0.5,color='#EE3224',label=df['concept'][0])

    #Create a bar with unDef data in position pos + some width buffer
    plt.bar([p + width for p in pos],df['unDef'],width,alpha=0.5,color='#F78F1E',label=df['concept'][1])

    #Create a bar with meth data in pos + some width buffer
    plt.bar([p + width*2 for p in pos],df['meth'],width,alpha=0.5,color='#FFC222',label=df['concept'][2])

    #Create a bar with different data in pos + some width buffer
    plt.bar([p + width*3 for p in pos],df['different'],width,alpha=0.5,color='#b29d87',label=df['concept'][2])

    #Set the y axis label
    ax.set_ylabel('#CG dinucleotides',fontsize=8)

    #Set the chart's title
    ax.set_title('Comparison methylation values and CG reported between %s and %s' %(label_first,label_second), fontsize=9)

    #Set the position of the x ticks
    ax.set_xticks([p + 1.5 * width for p in pos])

    #Set the labels for the x ticks
    ax.set_xticklabels(df['concept'],fontsize=8)

    #Setting the x-axis and y-axis limits
    plt.xlim(min(pos)-width, max(pos)+width*4)
    plt.ylim([0,max(df['unMeth'] + df['unDef'] + df['meth'] + df ['different'])])

    # Adding the legend and swowinf the plot
    plt.legend(['UnMethylated [0-0.3]','UnDefined [0.3-0.7]','Methylated [0.7-1]','Different'],fontsize=10)
    plt.grid()
    pylab.savefig(outpng)


#8. Add windows Unique for a given bed file
def addUniqueBedFile(bedFile,window):
    with open(bedFile , 'a') as bedFile:
        bedFile.write("%s\t%i\t%i\t%.2f\n" %( window[0],window[1],window[2],window[3]))


###################################################
#              ARGUMENTS PARSING                  #
###################################################

#1.Create object for argument parsinng
parser = argparse.ArgumentParser(prog="compareMethylationValues.py",description="Compares two bedgraph files outputting a set of comparison stats.")     

#1.1 Definition of input parameters
input_group = parser.add_argument_group('Inputs')
input_group.add_argument('--first', dest="first", metavar="FILE", help='First Bedgraph file.')
input_group.add_argument('--second', dest="second", metavar="FILE", help='Second Bedgraph file.')
input_group.add_argument('--first-label', dest="first_label", metavar="NAME", help='First label.')
input_group.add_argument('--second-label', dest="second_label", metavar="NAME", help='Second label.')
input_group.add_argument('-r',dest="relaxed", action='store_true', help='Relaxed Location Comparison.To check windowd positions in second file located within windows in first file.', default=False)

#1.2 Ouput
output_group = parser.add_argument_group('Output')
output_group.add_argument('--png', dest="png", metavar="FILE", help='Ouput PNG file.')
output_group.add_argument('--csv', dest="csv", metavar="FILE", help='Ouput CSV file.')

#Optional
output_group.add_argument('--outFirstUnique', dest="out_first_unique", metavar="FILE", help='Ouput first unique file.')
output_group.add_argument('--outSecondUnique', dest="out_second_unique", metavar="FILE", help='Ouput second unique file.')


#1.3. Argument parsing
args = parser.parse_args()

#1.4. Reset Unique bed files
if args.out_first_unique:
    open(args.out_first_unique, 'w').close()

if args.out_second_unique:
    open(args.out_second_unique, 'w').close()

###################################################
#              RUN COMPARISON                     #
###################################################

#2.1 Open file descriptors
fileFirst = open(args.first, 'r')
fileSecond = open(args.second, 'r')

#2.2 Define stats
#2.2.1 Unique to first
uniqueCGfirstUnmethylated, uniqueCGfirstUndefined,uniqueCGfirstMethylated = 0,0,0
#2.2.2 Unique to second
uniqueCGsecondUnmethylated, uniqueCGsecondUndefined, uniqueCGsecondMethylated = 0,0,0
#2.2.3 Shared First Second
sharedCGunmethylated, sharedCGundefined, sharedCGmethylated = 0,0,0
#2.2.4 Different Methylation Values First Second
differentCGmethValues = 0


#2.3 Control variables
flagReadFirst = True
flagReadSecond = True
previousChromosome = ""
windowFirst = []
windowSecond = []

#3. PARSING OF INPUT FILES
windowFirst = getWindow(getReads(fileFirst))
windowSecond = getWindow(getReads(fileSecond))

while len(windowFirst) == 5 or len(windowSecond) == 5:
    #3.1 Check if some of files has reached the end
    if len(windowFirst) < 5 and len(windowSecond) == 5:
        flagReadFirst = False
        flagReadSecond = True
        uniqueCGsecondUnmethylated, uniqueCGsecondUndefined, uniqueCGsecondMethylated = computeCounters(windowSecond[3],uniqueCGsecondUnmethylated,uniqueCGsecondUndefined,uniqueCGsecondMethylated)
    elif len(windowFirst) == 5 and len(windowSecond) < 5:
        flagReadFirst = True
        flagReadSecond = False
        uniqueCGfirstUnmethylated, uniqueCGfirstUndefined,uniqueCGfirstMethylated = computeCounters(windowFirst[3],uniqueCGfirstUnmethylated,uniqueCGfirstUndefined,uniqueCGfirstMethylated)
    #3.2 Both files have a window to be checked
    elif len(windowFirst) == 5 and len(windowSecond) == 5:
        #3.2.1 Update previous chromosome
        if previousChromosome == "":
            if windowFirst[0] < windowSecond[0]:
                previousChromosome = windowFirst[0]
            else:
                previousChromosome = windowSecond[0]
        elif windowFirst[0] == windowSecond[0]:
            previousChromosome = windowFirst[0]
        else:
            if windowFirst[0] < windowSecond[0]:
                previousChromosome = windowFirst[0]
            else:
                previousChromosome = windowSecond[0]

        #3.2.2 Comparing windows
        compare = compareWindows(arrayWindowFirst=windowFirst,arrayWindowSecond=windowSecond,previousContig=previousChromosome,relaxedCompare=args.relaxed)
  
        if compare == 0 or compare == 11 or compare == 12:
            #3.2.2.1 Windows are equal location
            if (compareMethValues(methValueFirst=windowFirst[3],methValueSecond=windowSecond[3])):
                #3.2.2.1.1 Windows have equal methylation value
                sharedCGunmethylated, sharedCGundefined, sharedCGmethylated = computeCounters(windowFirst[3],sharedCGunmethylated,sharedCGundefined,sharedCGmethylated)
            else:
                #3.2.2.1.2 Windows have different methylation value
                differentCGmethValues = differentCGmethValues + 1
       
            if compare == 0:
                flagReadFirst = True
                flagReadSecond = True
            elif compare == 11:
                #First Window Includes Second than just get new second window
                flagReadFirst = False
                flagReadSecond = True
            elif compare == 12:
                #Second Window Includes First than just get new first window
                flagReadFirst = True
                flagReadSecond = False
        elif compare == 1:
            #3.2.2.2 First window is smaller than second
            if not windowFirst[4]:
                uniqueCGfirstUnmethylated, uniqueCGfirstUndefined,uniqueCGfirstMethylated = computeCounters(windowFirst[3],uniqueCGfirstUnmethylated,uniqueCGfirstUndefined,uniqueCGfirstMethylated)

            flagReadFirst = True
            flagReadSecond = False

            if args.out_first_unique:
                if not windowFirst[4]:
                    addUniqueBedFile(args.out_first_unique,windowFirst)

        elif compare == 2:
            #3.2.2.3 Second window is smaller than first
            if not windowSecond[4]:
                uniqueCGsecondUnmethylated, uniqueCGsecondUndefined, uniqueCGsecondMethylated = computeCounters(windowSecond[3],uniqueCGsecondUnmethylated,uniqueCGsecondUndefined,uniqueCGsecondMethylated)

            flagReadFirst = False
            flagReadSecond = True

            if args.out_second_unique:
                if not windowSecond[4]:
                    addUniqueBedFile(args.out_second_unique,windowSecond)

    #3.3 Perform new read calls
    if flagReadFirst:
        windowFirst = getWindow(getReads(fileFirst))
    if flagReadSecond:
        windowSecond = getWindow(getReads(fileSecond))

#4. Print results
print("Unique CG to %s:\n" %(args.first_label))  
print("    Total: %i\n" %(uniqueCGfirstUnmethylated+uniqueCGfirstUndefined+uniqueCGfirstMethylated))
print("        Unmethylated (0-0.3)  : %i\n" %(uniqueCGfirstUnmethylated))
print("        Undefined    (0.3-0.7): %i\n" %(uniqueCGfirstUndefined))
print("        Methylated   (0.7-1)  : %i\n" %(uniqueCGfirstMethylated))

print("\n")
print("Unique CG to %s:\n" %(args.second_label))  
print("    Total: %i\n" %(uniqueCGsecondUnmethylated+uniqueCGsecondUndefined+uniqueCGsecondMethylated))
print("        Unmethylated (0-0.3)  : %i\n" %(uniqueCGsecondUnmethylated))
print("        Undefined    (0.3-0.7): %i\n" %(uniqueCGsecondUndefined))
print("        Methylated   (0.7-1)  : %i\n" %(uniqueCGsecondMethylated))

print("\n")
print("Shared CG between %s and %s:\n" %(args.first_label,args.second_label)) 
print("    Total Equal Methylated Value: %i\n" %(sharedCGunmethylated+sharedCGundefined+sharedCGmethylated))
print("        Unmethylated (0-0.3)  : %i\n" %(sharedCGunmethylated))
print("        Undefined    (0.3-0.7): %i\n" %(sharedCGundefined))
print("        Methylated   (0.7-1)  : %i\n" %(sharedCGmethylated))
print("    Total Different Methylated Value: %i\n" %(differentCGmethValues))

#4.1 Run Plotting
plotsStats((uniqueCGfirstUnmethylated,uniqueCGfirstUndefined,uniqueCGfirstMethylated),\
            (uniqueCGsecondUnmethylated,uniqueCGsecondUndefined,uniqueCGsecondMethylated),\
            (sharedCGunmethylated,sharedCGundefined,sharedCGmethylated),\
            differentCGmethValues,args.first_label,args.second_label,args.png)


#4.2 Save Results as CSV file
with open(args.csv , 'w') as csvFile:
    csvFile.write("Uniq CG to %s,%i,%i,%i,%i\n" \
                  %(args.first_label,(uniqueCGfirstUnmethylated+uniqueCGfirstUndefined+uniqueCGfirstMethylated),uniqueCGfirstUnmethylated,uniqueCGfirstUndefined,uniqueCGfirstMethylated))
    csvFile.write("Uniq CG to %s,%i,%i,%i,%i\n" \
                  %(args.second_label,(uniqueCGsecondUnmethylated+uniqueCGsecondUndefined+uniqueCGsecondMethylated),uniqueCGsecondUnmethylated,uniqueCGsecondUndefined,uniqueCGsecondMethylated))
    csvFile.write("Shared between %s and %s,%i,%i,%i,%i\n" \
                  %(args.first_label,args.second_label,(sharedCGunmethylated+sharedCGundefined+sharedCGmethylated),sharedCGunmethylated,sharedCGundefined,sharedCGmethylated))
    csvFile.write("Shared between %s and %s different meth.,%i\n"\
                  %(args.first_label,args.second_label,differentCGmethValues))


#5. Close files descriptors
fileFirst.close()
fileSecond.close()

























