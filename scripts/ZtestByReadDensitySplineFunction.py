#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Mon Jan 30 13:27:14 2017

an example test call from terminal:
python ./ZtestByReadDensitySplineFunction.py --InputFile /Users/benfair/Documents/School/GradSchool/Rotation2/LabNotebook_20140619-Present/20161206_ScreenSeqAnalysis-Srp2Knockout/ZtestByReadDensitySplineFunctionTestData.txt --OutputFile /dev/null --NumeratorColumnName SumUnspliced --DenominatorColumnName SumSpliced --PlotMA myma.pdf

See main() docstring for more information on this script

@author: benfair
"""

import sys
import numpy as np
import scipy.stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse
import math
from scipy.interpolate import UnivariateSpline
import matplotlib.lines as mlines
import warnings
import qvalue

def Log10Ticks(MyMin, MyMax):
    '''
    Helper function for adding "log10 scale ticks" on a graph whose values are already log10 transformed. For example, on MA plots I often plot log10(read count) on the x-axis. If I want to add tick marks at the 100, 200, 300, ..., 1000, 2000, etc, this function will help you do that
    
    Returns a tuple of two objects. First object is a ist of tick positions. Second is a list of tick labels where only powers of 10 are labelled (no tick labels at 200, 300, etc.; only at 10, 100, 1000, etc.)
    
    To replace tick marks with log10 scale ticks to a graph in log10 space, do this:
    matplotlib.pyplot.xticks(*Log10Ticks(*axes.get_xlim)))
    '''
    np.set_printoptions(precision=1)
    mylist = range(int(math.floor(MyMin)), int(math.ceil(MyMax)))
    listout =[]
    ticklabelsout=[]
    for i in mylist:
        # listout += range(10**i, 10**(i+1), 10**i)
        listout += np.linspace(10**i, 10**(i+1), num=9, endpoint=False).tolist()
    listout += np.around([10**math.ceil(MyMax)]).tolist()
    for i in listout:
        if np.log10(i).is_integer():
            if i >= 1: ticklabelsout.append(str(int(i)))
            else: ticklabelsout.append(str(i))
        else:
            ticklabelsout.append('')
    listoutlogspace = np.log10(listout)
    return listoutlogspace, ticklabelsout

def HeatScatter(ax, x, y, **kwargs):
    """
     A helper function to make a scatter plot where the points are colored according to the point density. Similar to the 'Heatscatter' function from the R package 'LSD'.

    Parameters
    ----------
    ax : Axes
        The axes to draw to

    data1 : array
       The x data

    data2 : array
       The y data

    kwargs
       Will accept all the same arguments as matplotlib.pyplot.scatter()
       
    Returns
    -------
    out :
        Adds all the same artists as calling ax.scatter()
    """
    if 'c' not in kwargs: #unless color ('c' parameter) is specified, use rainbow color according to point density
        #Before calculating point density we must make sure x and y are masked in the same points

        # Calculate the point density
        xy = np.vstack([x,y])
        z = scipy.stats.gaussian_kde(xy)(xy)
        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        kwargs['c'] = z
        kwargs['edgecolor'] = ''
    out = ax.scatter(x, y, **kwargs)
    return out

def main(InputFile,
    OutputFile = None,
    MinimumReadThreshold = 1000,
    MeanSplineFunctionOrder = 2,
    StdSplineFunctionOrder = 1,
    NumberBins = 10,
    BinSpacing = 'EvenNumberSamplesPerBin',
    TestType = 'TwoSidedTest',
    MultipleHypothesisCorrection = 'BenjaminiHochberg',
    NoSplineFitting = False,
    FDR = 0.05,
    PlotMA = None,
    PlotMAPercent = None,
    PlotQQ = None,
    LabelSignificantPointsInPlot = False,
    ShowSplineFitConfidenceInterval = False,
    CopyInputFileContentsToOutput = False,
    IdentifierColumnName = None,
    NumeratorColumnName = [],
    DenominatorColumnName = [],
    Verbose = False,
    **PlotParameterDict):
    """
    The main function of this script is used to test for significance of a ratio of two measured numbers value when the expected mean ratio and error (standard deviation) changes as a function of the the total signal (or more precisely, the sum of the measured numbers). In practice, we use to this to find samples which have a significantly changed ratio of unspliced/spliced reads as a function of total read count. This is the method used to identify significant changes in splicing among S. pombe knockout strains in Larson et al (2016) G3. The procedure is to first make a scatter plot of log(Numerator/Denominator) vs log(Numerator+Denominator). This is basically an MA-plot similar to those often used in microarray gene-expression or RNA-seq data. Noise in the measured ratio [log(Numerator/Denominator)] usually decreases as total signal [log(Numerator+Denominator)] increases. The points within this scatter plot (excluding those points below a total signal minimum threshold) are binned by total signal. Within each bin, the mean of the ratio and the standard deviation of the ratio are calculated and used as points to which a smoothed function is fit to using splines. The resulting spline functions (one function for measured ratio as a function of total signal, and another function for standard deviation as a function of total signal) are used to form a null hypothesis standard distribution at any given the total signal level. The Z-score (number of standard deviations from the mean given total signal) is then calculated for each point, which is then converted to a P-value, which is among the values reported for each sample in the output. The P-values are adjusted for multiple hypothesis correction by the Benjamini-Hochberg method. This script does not take into account replicate measurements. What we have done Larson et al (2016) to identifiy significant splicing changes, is simply combine biological replicates before using this statistical test by summing spliced read counts amongst replicates, and summing read counts for unspliced amongst replicates.
The required input is a tab-delimited text file where each row corresponds to a different sample and the columns correspond to the sample name, the numerator (UnsplicedReadCount) and a denominator (SplicedReadCount), respectively. Other columns may be present as well but they will not be used. The first row must be a header row which gives the titles of the columns. Optionally, this main function can also output a pdf of the MA-plot which helps visualize how the dataset looks and how the spline functions estimate a null hypothesis.
     
    DEPENDENCIES:
    This python script requires numpy, scipy, and matplotlib
    """
    warnings.simplefilter(action = "ignore", category = FutureWarning)
    np.seterr(divide='ignore')
    if InputFile == 'stdin':
        InputFile = sys.stdin
    TableIn = np.genfromtxt(InputFile, dtype=None, delimiter='\t', missing_values="", names=True, usemask=True)
    if not IdentifierColumnName: IdentifierColumnName = TableIn.dtype.names[0]
    if not NumeratorColumnName: NumeratorColumnName.append(TableIn.dtype.names[1])
    if not DenominatorColumnName: DenominatorColumnName.append(TableIn.dtype.names[2])
    if not set([IdentifierColumnName] + NumeratorColumnName + DenominatorColumnName) <= set(TableIn.dtype.names):
        print('ERROR: Could not identify column names supplied. If you provided column name(s) in the IdentifierColumnName, NumeratorColumnName, or DenominatorColumnName or arguments, make sure that they exactly matches following column names:')
        print('\n'.join(TableIn.dtype.names))
        return
    
    if Verbose: print('Succesfully parsed input arguments and read InputFile into memory. Doing math now...')
    NumeratorArray = sum([TableIn[i] for i in NumeratorColumnName])
    DenominatorArray = sum([TableIn[i] for i in DenominatorColumnName])

    
    x = np.ma.log10(NumeratorArray + DenominatorArray)
    y = np.log2(np.ma.true_divide(NumeratorArray,DenominatorArray))
    x.mask = y.mask
    
    yOverReadThreshold = np.ma.masked_where((x<np.log10(MinimumReadThreshold)), y)
    xOverReadThreshold = np.ma.masked_where((x<np.log10(MinimumReadThreshold)), x)
#     yOverReadThreshold = np.ma.masked_where((x<np.log10(MinimumReadThreshold))&(NumeratorArray<10), y)
#     xOverReadThreshold = np.ma.masked_where((x<np.log10(MinimumReadThreshold))&(NumeratorArray<10), x)
    
    if BinSpacing == 'EvenNumberSamplesPerBin':
        bin_edges = scipy.stats.mstats.mquantiles(xOverReadThreshold.compressed(), prob=np.linspace(0,1,num=NumberBins+1, endpoint=True))
    elif BinSpacing == 'EvenlySpacedBins':
        bin_edges = np.linspace(xOverReadThreshold.compressed().min(), xOverReadThreshold.compressed().max(), num=NumberBins+1)
    bin_std = scipy.stats.binned_statistic(xOverReadThreshold.compressed(), yOverReadThreshold.compressed(), statistic='std', bins=bin_edges)[0]
    bin_mean = scipy.stats.binned_statistic(xOverReadThreshold.compressed(), yOverReadThreshold.compressed(), statistic='mean', bins=bin_edges)[0]
    bin_centers = [(bin_edges[i] + bin_edges[i+1])/2 for i in range(0,len(bin_edges)-1)] #Find the bin_centers by averaging the (i)th and the (i+1)th bin edge
    bin_centers[0] = bin_edges[0] #Replace the first and last bins_centers (which will be used as points for spline-fitting) to the bin edges.
    bin_centers[-1] = bin_edges[-1]
    std_as_function_of_counts = UnivariateSpline(bin_centers, bin_std, k=StdSplineFunctionOrder)
    mean_as_function_of_counts = UnivariateSpline(bin_centers, bin_mean, k=MeanSplineFunctionOrder)
    yTransformedToPercent = 100*(2**(y-mean_as_function_of_counts(x)+np.ma.median(yOverReadThreshold)))/(1+2**(y-mean_as_function_of_counts(x)+np.ma.median(yOverReadThreshold)))
    Zscores = (yOverReadThreshold-mean_as_function_of_counts(xOverReadThreshold))/std_as_function_of_counts(xOverReadThreshold)
    if NoSplineFitting:
        Zscores = yOverReadThreshold - np.ma.mean(yOverReadThreshold)/yOverReadThreshold.std()
    if TestType == 'TwoSidedTest':
        pvalues = np.ma.masked_array(scipy.stats.norm.sf(abs(Zscores)), mask=xOverReadThreshold.mask) * 2
    elif TestType == 'RightSidedTest':
        pvalues = np.ma.masked_array(scipy.stats.norm.sf(Zscores), mask=xOverReadThreshold.mask)
    elif TestType == 'LeftSidedTest':
        pvalues = np.ma.masked_array(scipy.stats.norm.cdf(Zscores), mask=xOverReadThreshold.mask)
    if MultipleHypothesisCorrection == 'BenjaminiHochberg':
        CorrectedP = pvalues*pvalues.count()/scipy.stats.mstats.rankdata(pvalues)
    elif MultipleHypothesisCorrection == 'Qvalue':
        CorrectedP = qvalue.estimate(pvalues)
    SignificantCount = CorrectedP.compressed()[CorrectedP.compressed()<FDR].size
    if Verbose: print('{} significant after multiple hypothesis correction...)\nWriting output file to now...').format(SignificantCount)
    
    if OutputFile:
        fhOut = open(OutputFile, 'w')
        fhIn = open(InputFile, 'rU')
        for i, line in enumerate(fhIn):
            if i == 0:
                if CopyInputFileContentsToOutput:
                    fhOut.write(line.strip('\n') + '\t')
                fhOut.write('Sample\tReadCount\tRatio\tZscore\tPvalue\tCorrectedPvalue\tRatioAsNormalizedPercent')
            else:
                indexPos = i-1
                fhOut.write('\n')
                if CopyInputFileContentsToOutput:
                    fhOut.write(line.strip('\n') + '\t')
                fhOut.write('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(TableIn[IdentifierColumnName][indexPos],x[indexPos],y[indexPos],Zscores[indexPos],pvalues[indexPos],CorrectedP[indexPos],yTransformedToPercent[indexPos]))
        fhIn.close()
        fhOut.close()

    if PlotMA:
        print 'plotmay', y
        if Verbose: print('Generating MA-plot...')
        fig, ax = plt.subplots(1, 1)
        mpl.rcParams['mathtext.fontset'] = 'custom'
        mpl.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
        mpl.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
        mpl.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
        HeatScatter(ax, x, y, **PlotParameterDict)
        plt.xticks(*Log10Ticks(*ax.get_xlim()))
        print ax.get_ylim()
        yticks, yticklabels = Log10Ticks(*[i/(np.log(10)/np.log(2)) for i in ax.get_ylim()])
        yticks = np.array(yticks) * (np.log(10)/np.log(2))
        plt.yticks(yticks, yticklabels)
        if ShowSplineFitConfidenceInterval:
            xnew = np.linspace(bin_centers[0], bin_centers[-1], num=100)
            plt.fill_between(xnew, mean_as_function_of_counts(xnew)+std_as_function_of_counts(xnew)*2,mean_as_function_of_counts(xnew)-std_as_function_of_counts(xnew)*2, color='gray', alpha=0.5) #Shade two standard deviations above and below the mean as a function of counts
            plt.plot(xnew,mean_as_function_of_counts(xnew), c='black')
            NullHypothesisLine = mlines.Line2D([], [], color='black', marker='o', label=r'$\rmH_0\ with\ points\ used\ for\ spline\ interpolation$')
            plt.legend(handles=[NullHypothesisLine], loc=3, fontsize=10)
            plt.scatter(bin_centers,bin_mean, marker='o', c='black', edgecolor='black')
        plt.xlabel('Read Count')
        plt.ylabel(r'$\frac{%s}{%s}$' %('+'.join(NumeratorColumnName).replace('_', ''), '+'.join(DenominatorColumnName).replace('_', '')))
        if LabelSignificantPointsInPlot:
            SignificantPoints = []
            for i, CorrectedPvalue in enumerate(np.nditer(CorrectedP)):
                if CorrectedPvalue < FDR and x[i]>=np.log10(MinimumReadThreshold):
                    SignificantPoints.append((x[i],y[i],TableIn[IdentifierColumnName][i]))
            for x, y, name in SignificantPoints:
                plt.text(x,y,name,size=6)
        plt.savefig(PlotMA)

    if PlotMAPercent:
        if Verbose: print('Generating MA-plot...')
        fig, ax = plt.subplots(1, 1)
        mpl.rcParams['mathtext.fontset'] = 'custom'
        mpl.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
        mpl.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
        mpl.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
        print 'hello', y
#         newy = np.ma.log2(100*(2**(y-mean_as_function_of_counts(x)+np.ma.median(yOverReadThreshold)))/(1+2**(y-mean_as_function_of_counts(x)+np.ma.median(yOverReadThreshold))))
        newy = 100*(2**(y-mean_as_function_of_counts(x)+np.ma.median(yOverReadThreshold)))/(1+2**(y-mean_as_function_of_counts(x)+np.ma.median(yOverReadThreshold)))

        print len(newy), len(newy.compressed())
        fig, ax = plt.subplots(1,1)
        ax.scatter(x, newy, c='gray', alpha=0.1, edgecolors='none', s=12)
        plt.xticks(*Log10Ticks(*ax.get_xlim()))
        plt.xlim(3,5.5)
        plt.ylim(0,16)
#         yticks, yticklabels = Log10Ticks(*[i/(np.log(10)/np.log(2)) for i in ax.get_ylim()])
#         yticks = np.array(yticks) * (np.log(10)/np.log(2))
#         plt.yticks(yticks, yticklabels)
        ax.tick_params(direction='in', top='on', right='on')
        plt.xlabel('Read Count')
        plt.ylabel("Percent Unspliced")
#         plt.ylabel(r'$\frac{%s}{%s}$' %('+'.join(NumeratorColumnName).replace('_', ''), '+'.join(DenominatorColumnName).replace('_', '')))
        if LabelSignificantPointsInPlot:
            SignificantPoints = []
            for i, CorrectedPvalue in enumerate(np.nditer(CorrectedP)):
                if (CorrectedPvalue < FDR and x[i]>=np.log10(MinimumReadThreshold)) or TableIn[IdentifierColumnName][i]=='102B12':
                    SignificantPoints.append((x[i],newy[i],TableIn[IdentifierColumnName][i]))
            for myx, myy, name in SignificantPoints:
                ax.scatter(myx,myy,c='white', edgecolor='gray', linewidths=1.5, s=12)
#                 plt.text(myx,myy,name,size=1)
        plt.savefig(PlotMAPercent)

    if PlotQQ:
        fig, ax = plt.subplots(1, 1)
        (osm, osr), (slope, intercept, r) = scipy.stats.probplot(Zscores.compressed(),plot=ax,dist=scipy.stats.norm,rvalue=True, fit=True)
#         print r**2
        x = np.linspace(*ax.get_xlim())
        line, = ax.plot(x, x)
        line.set_label('y=x line')
        ax.legend()
        plt.savefig(PlotQQ)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage = main.__doc__)
    parser.add_argument('--InputFile', metavar = "[InputFilepath]", help="REQUIRED. The input file. Should contain a header line and at minmum a column with the sample name, a column for the numerator and a column for the denominator. If piping from stdin, using the filename 'stdin'.", required=True)
    parser.add_argument('--OutputFile', metavar = "[OutputFilepath]", help="OPTIONAL. File that describes barcodes and filepaths to separate into. If none is provided, won't write out data.", required=False)
    parser.add_argument("--Verbose", help="OPTIONAL. Output print statements of progress to stdout", action="store_true")
    parser.add_argument("--MinimumReadThreshold", metavar = "[integer]", help="OPTIONAL. Minimum number of reads in a sample to consider for spline fitting and statistical testing, Default = 1000", type=int, default=1000)
    parser.add_argument("--MeanSplineFunctionOrder", metavar = "[integer]", help="OPTIONAL. Polynomial order of the spline fit to the mean of bins used to estimate the  mean under the null hypothesis at a given read depth. Default = 2", type=int, default=2)
    parser.add_argument("--StdSplineFunctionOrder", metavar = "[integer]", help="OPTIONAL. Polynomial order of the spline fit to the standard deviation of bins used to estimate the standard deviation under the null hypothesis at a given read depth. Default = 1", type=int, default=1)
    parser.add_argument("--NumberBins", metavar = "[integer]", help="OPTIONAL. Number of ReadCount bins to divide the dataset into to create points for spline fitting. Samples below the MinimumReadThreshold are not considered part of the dataset. Default = 10", type=int, default=10)
    parser.add_argument("--BinSpacing", help="OPTIONAL. How to determine bin spacing along ReadCount axis, Default = EvenNumberSamplesPerBin", choices=['EvenNumberSamplesPerBin','EvenlySpacedBins'], default='EvenNumberSamplesPerBin')
    parser.add_argument("--TestType", help="OPTIONAL. TwoSidedTest: And increase or decrease in log(numerator/denominator) can be called significant. RightSidedTest: Only an increase in log(numerator/denominator) can be called significant). LeftSidedTest: Only an decrease in log(numerator/denominator) can be called significant). Default = TwoSidedTest", choices=['TwoSidedTest','RightSidedTest', 'LeftSidedTest'], default='TwoSidedTest')
    parser.add_argument("--MultipleHypothesisCorrection", help="OPTIONAL. BenjaminiHochberg: Implemented as on the wiki-page. Qvalue: Implemented using the module 'qvalue' which must be installed. Default = BenjaminiHochberg", choices=['BenjaminiHochberg','Qvalue'], default='BenjaminiHochberg')
    parser.add_argument("--NoSplineFitting", help="OPTIONAL. No spline fitting, just get Z statistic from log(numerator/denominator)", action="store_true")
    parser.add_argument("--FDR", metavar = "[float]", help="OPTIONAL. False discovery rate. Default = 0.05", type=float, default=0.05)
    parser.add_argument("--PlotMA", metavar = "[filepath.pdf]", help="OPTIONAL. If a output filepath is supplied, create a MA-plot of the data and save it as the provided filepath. The type of image depends on the file extension (.svg, .png, .pdf) given. Default = None (No plot created)", default=None)
    parser.add_argument("--PlotMAPercent", metavar = "[filepath.pdf]", help="OPTIONAL. If a output filepath is supplied, create a MA-plot of the data and save it as the provided filepath. The type of image depends on the file extension (.svg, .png, .pdf) given. Default = None (No plot created)", default=None)
    parser.add_argument("--PlotQQ", metavar = "[filepath.pdf]", help="OPTIONAL. If a output filepath is supplied, create a normal QQ-plot of the Z-scores and save it as the provided filepath. The type of image depends on the file extension (.svg, .png, .pdf) given. Default = None (No plot created)", default=None)
    parser.add_argument("--PlotParameter", metavar = "[key]=[value]", help="OPTIONAL. Plotting parameters to pass to MA-plotting function. The same parameters that can be used with the matplotlib.pyplot.scatter() function can be used. Each key=value must be preceded by --PlotParameterDictionary. For example, if you wanted to the plot points to be blue and transparent, you would include --PlotParameter c=blue --PlotParameter alpha=0.5", action='append', type=lambda kv: kv.split("="), dest='PlotParameterList', default = None)
    parser.add_argument("--CopyInputFileContentsToOutput", help="OPTIONAL. Copy the file contents of the input file to the output file, while appending the same columns that would otherwise be in the output file.", action="store_true")
    parser.add_argument("--LabelSignificantPointsInPlot", help="OPTIONAL. Label samples which are below the acceptable FDR with a sample name (samples names are provided by the value of IdentifierColumn of the InputFile specified by --IdentifierColumnName)", action="store_true")
    parser.add_argument("--ShowSplineFitConfidenceInterval", help="OPTIONAL. Show a null hypothesis line with a shaded region that represent the spline-fit std and mean within each bin", action="store_true")
    parser.add_argument("--IdentifierColumnName", metavar = "[Column_name_string]", help="OPTIONAL. Column name in the header of InputFile that corresponds to sample names. If odd characters (+, space, etc) are used in the header, they must be replaced with '_'. Default = If none is provided, this script defaults to using first column")
    parser.add_argument("--NumeratorColumnName", metavar = "[Column_name_string]", help="OPTIONAL. Column name in the header of InputFile that corresponds to numerator to use when calculating a ratio. Normally the numerator corresponds to a count of unspliced reads. If odd characters (+, space, etc) are used in the header, they must be replaced with '_'. If more than one argument is provided, will use the sum of read counts in those columns. Default = If none is provided, this script defaults to using second column", nargs='+', default=[])
    parser.add_argument("--DenominatorColumnName", metavar = "[Column_name_string]", help="OPTIONAL. Column name in the header of InputFile that corresponds to denominator to use when calculating a ratio. Normally the numerator corresponds to a count of spliced reads. If odd characters (+, space, etc) are used in the header, they must be replaced with '_'. If more than one argument is provided, will use the sum of read counts in those columns. Default = If none is provided, this script defaults to using third column", nargs='+', default=[])
    args = parser.parse_args()
    PlotParameterDict={}
    if args.PlotParameterList:
        for k,v in dict(args.PlotParameterList).items():
            try:
                PlotParameterDict[k] = float(v)
            except ValueError:
                PlotParameterDict[k] = v
    ArgumentsDict = vars(args)
    del ArgumentsDict['PlotParameterList']
    AllMainArguments = dict(list(PlotParameterDict.items()) + list(ArgumentsDict.items()))
    print AllMainArguments
    main(**AllMainArguments)
