#!/usr/bin/env python

# Purpose:       Get individual genotypes and return it in VCF format
# 
# Usage:         python GetIndividualGenotypesVCF.py -c config.yml [OPTIONS]
# 
# Author:        et85, etsai@bwh.harvard.edu



from argparse import ArgumentParser
import Bio.bgzf # import biopython
import drmaa
import os, re, shutil, subprocess, sys, time, yaml

####################################################################################
#############################   Function Definitions   #############################
####################################################################################
def AssignVariableFromConfig(thisYaml, I_list, var1, var2, desc):
    try:
        thisVariable = thisYaml[var1][var2]
        
        # if this is a data file, please make sure it exists
        if(var1 == 'Files'):
            if(I_list):
                I_exists = True
                for thisFile in thisVariable:
                    if(not os.path.isfile(thisFile)):
                        print "Invalid " + desc + " path: " + thisFile
                        sys.exit()
                return thisVariable
            elif(os.path.isfile(thisVariable)):
                    return thisVariable
            else:
                print "Invalid " + desc + " path: " + thisVariable
                sys.exit()
        elif(desc.startswith("I_")):
            return thisVariable
        else:
            return thisVariable
    except:
        print "Please assign " + desc + " by adding this parameter or editing the config file"
        raise

def CheckSampleArray(thisFile,thisIdDictList):
    thisBatchList = [[] for i in range(len(thisIdDictList))]
    
    thisFH = open(thisFile,'r')
    # Decode Byte Order Mark
    s = thisFH.read()
    u = s.decode("utf-16")
    s = u.encode("utf-8")
    s = s.rstrip()
    thisIdList = s.split('\r\n')[2:]     # skip the first two header rows
    for thisBiobankID in thisIdList:
        count = 0
        # iterate through the different IdDictList until we find one
        for i in range(len(thisIdDictList)):
            if(thisIdDictList[i].has_key(thisBiobankID)):
                thisBatchList[i].append(thisBiobankID)
                count += 1
        if(count == 0):
            print thisBiobankID+': Cannot find sentrix id for this subject id'
        if(count > 1):
            print thisBiobankID+': Multiple biobank samples found for this subject id'
    
    thisNsampList = []
    for i in range(len(thisIdDictList)):
        thisNsampList.append(len(thisBatchList[i]))
    
    return thisBatchList, thisNsampList

def FormatUsername(thisBaseName):
    thisUsername = thisBaseName.split("_")[1]
    return thisUsername

def MakeConcordanceTable(idtablepath):
    thisDict = {}
    thisFH = open(idtablepath, 'r')
     
    for thisLine in thisFH.readlines()[1:]:      # skip first line (header)
        thisLine = thisLine.rstrip()
        thisRead = thisLine.split('\t')
        thisDict[thisRead[0]] = thisRead[1]
    return thisDict

def SetDefaultOutDir(thisBaseDir,thisBaseName):
    if (thisBaseDir == ''):
        thisBaseDir = '/data/ppm/biobank/release'
    
    thisUsername = FormatUsername(thisBaseName)

    #Check to make sure you're on panasas
    if(os.path.exists(thisBaseDir)):
        subprocess.call(['mkdir', '-p', thisBaseDir + '/' + thisUsername])
    else:
        sys.exit('Cannot find path: '+thisBaseDir)
    return thisBaseDir + '/' + thisUsername + '/' + thisBaseName

# thisJobIdList = ExtractPedDataset(s, thisSampleToBatchList, thisNSampCovarDict, thisIdDictList, pedDirList, outDir, tempDir, batchNamesList, batchVariableList, extractPedScript)

def ExtractVcfDataset(s, thisIdList, thisNSampCovarDict, thisIdDictList, thisOutDir, thisTempDir, thisBatchNameList, thisBatchVariableList, thisExtractVcfScript):
    thisOutfilePathList = []
    thisSnpinfoPathList = []
    thisJobIdList = []
    
    # Make output dir if it doesn't exist.
    if os.path.exists(thisOutDir):
        # temporarily disable check for empty outdir
        pass
#         if(os.listdir(thisOutDir) != []):
#             sys.exit("This out dir is not empty. Please delete the contents in the out dir before trying again:\n" + thisOutDir)
    else:
        os.makedirs(thisOutDir)
    
    # Make drmaa template
    print('creating job template for extract_vcf')
    jt = s.createJobTemplate()
     
    thisBsubOut = thisTempDir + '/extract_vcf_%J_%I.out'
    thisBsubErr = thisTempDir + '/extract_vcf_%J_%I.err'
    jt.nativeSpecification = '-R \'rusage[mem=8192]\' ' + \
                             '-M 12288 ' + \
                             '-q pcpgmwgs ' + \
                             '-J extract_vcf ' + \
                             '-o ' + thisBsubOut + ' ' \
                             '-e ' + thisBsubErr
    jt.remoteCommand = thisExtractVcfScript
    
    # Initialize covariate file if there are samples
    for thisKey in sorted(thisNSampCovarDict.keys()):
        if(thisNSampCovarDict[thisKey][0] > 0):
            thisCovarFH = open(thisOutDir+'/'+thisNSampCovarDict[thisKey][1]+'.covar.txt','w')
            thisCovarFH.write('SentrixID\tBiobankID\tBatch\n')
            thisMergeFH = open(thisNSampCovarDict[thisKey][3],'w')
        thisCovarFH.close()
        thisMergeFH.close()
    
    for i in range(len(thisBatchVariableList)):
        thisBatchKey = thisBatchVariableList[i][1].split('_')[0]
        thisBatchBase = thisNSampCovarDict[thisBatchKey][1]
        
        # skip this iteration if there are no samples for it
        if(len(thisIdList[i]) == 0):
            print "No requests in this list"
            continue
        
        # make filelist in outdir
        thisCovarFH = open(thisOutDir+'/'+thisBatchBase+'.covar.txt','a')
        thisMergeFH = open(thisNSampCovarDict[thisBatchKey][3],'a')
        thisSampListFile = thisTempDir+'/'+thisBatchVariableList[i][0]+'_' + str(thisBatchVariableList[i][2]) + '.samplist.txt'
        thisSampFH = open(thisSampListFile,'w')
        
        for thisBiobankID in thisIdList[i]:
            if(thisIdDictList[i].has_key(thisBiobankID)):
                thisSentrixID = thisIdDictList[i][thisBiobankID]
                thisSampleFile = thisSentrixID + '-' + thisBiobankID
                thisCovarFH.write(thisSentrixID+'\t'+thisBiobankID+'\t'+str(thisBatchVariableList[i][2])+'\n')
                thisSampFH.write(thisSentrixID+'-'+thisBiobankID+'\n')
                print thisSentrixID + '-' + thisBiobankID
            else:
                print "Cannot find key."
        
        thisOutfilePath = ''
        thisOutfilePath = thisTempDir + '/' + thisBatchVariableList[i][0] + '_' + str(thisBatchVariableList[i][2])
        thisMergeFH.write(thisOutfilePath + '.vcf.gz\n')
        
        thisCovarFH.close()
        thisSampFH.close()
        thisMergeFH.close()
        
        if(thisNSampCovarDict[thisBatchKey][0] > 0):
            thisBatchVcfFile = thisBatchVariableList[i][3]
            
            # run drmaa jobs
            jt.args = [thisSampListFile,thisOutfilePath,thisBatchVcfFile]
            jobid = s.runJob(jt)
            thisJobIdList.append(jobid)
        
    return thisJobIdList


def MergeVcfFiles(s, thisJobIdList, thisNSampCovarDict, batchVariableList, thisTempDir, thisMergeVcfScript):
    thisMergeVcfJobIdList = []
    
    print('creating job template for merge_vcf')
    jt = s.createJobTemplate()
      
    thisBsubOut = thisTempDir + '/merge_vcf_%J_%I.out'
    thisBsubErr = thisTempDir + '/merge_vcf_%J_%I.err'
    thisBsub = '-R \'rusage[mem=8192]\' ' + \
                             '-M 12288 ' + \
                             '-q pcpgmwgs ' + \
                             '-J merge_vcf ' + \
                             '-o ' + thisBsubOut + ' ' \
                             '-e ' + thisBsubErr + ' '
    
    if(len(thisJobIdList) > 0):
        thisDependency = '-w \''
        for i in range(len(thisJobIdList)):
            if(i == 0):
                thisDependency += 'done('+str(thisJobIdList[i])+')'
            else:
                thisDependency += ' && done('+str(thisJobIdList[i])+')'
        thisDependency += '\''
        thisBsub += thisDependency
    
    jt.nativeSpecification = thisBsub
    jt.remoteCommand = thisMergeVcfScript
    
    for thisKey in thisNSampCovarDict.keys():
        thisOutfilePath = thisNSampCovarDict[thisKey][2]
        thisSampMergeListPath = thisNSampCovarDict[thisKey][3]
        thisSnpinfoFile = thisNSampCovarDict[thisKey][4]
        
        # run drmaa jobs
        jt.args = [thisSampMergeListPath,thisOutfilePath,thisSnpinfoFile]
        jobid = s.runJob(jt)
        thisMergeVcfJobIdList.append(jobid)
    
    return thisMergeVcfJobIdList

# def MakeMergedReadmeFile(thisReadmeFile, thisOutDir, thisBaseNameList, thisNsampList):
def MakeMergedReadmeFile(thisReadmeFile, thisOutDir, thisNSampCovarDict):
    if not os.path.isfile(thisReadmeFile):
        sys.exit("README file does not exist. Please check file path: " + thisReadmeFile)
    
    thisFormDate = time.strftime('%Y-%m-%d')
    thisOutFH = open(thisOutDir + '/README.txt',"w")         
    
    thisBaseNameList = []
    for thisKey in thisNSampCovarDict.keys():
        thisBaseNameList.append(thisNSampCovarDict[thisKey][1])
    
    with open(thisReadmeFile, "r") as thisReadmeFH:
        for thisLine in thisReadmeFH.readlines():
            thisLine = thisLine.replace('INSERTFILENAME_VCFGZ', '.vcf.gz, '.join(thisBaseNameList)+'.vcf.gz')
            thisLine = thisLine.replace('INSERTFILENAME_VCFTBI', '.vcf.gz.tbi, '.join(thisBaseNameList)+'.vcf.gz.tbi')
            thisLine = thisLine.replace('INSERTFILENAME_SNPINFO', '.snpinfo.txt, '.join(thisBaseNameList)+'.snpinfo.txt')
            thisLine = thisLine.replace('INSERTFILENAME_COVAR', '.covar.txt, '.join(thisBaseNameList)+'.covar.txt')
            thisLine = thisLine.replace("INSERTDATE", thisFormDate)
            thisOutFH.write(thisLine)
    
    thisOutFH.close()

def BundleToDropbox(s, thisOutPath, thisTempDir, I_readonly, I_compress, I_md5sum, I_debug, thisJobIdList, thisBundleDataScript):
    thisRead = thisOutPath.split("/")
    thisOutDirName = thisRead[-1]
    
    print('creating job template for bundle_data')
    jt = s.createJobTemplate()
      
    thisBsubOut = thisTempDir + '/bundle_data_%J_%I.out'
    thisBsubErr = thisTempDir + '/bundle_data_%J_%I.err'
    thisBsub = '-R \'rusage[mem=8192]\' ' + \
                             '-M 12288 ' + \
                             '-q pcpgmwgs ' + \
                             '-J bundle_data ' + \
                             '-o ' + thisBsubOut + ' ' \
                             '-e ' + thisBsubErr + ' '
    
    if(len(thisJobIdList) > 0):
        thisDependency = '-w \''
        for i in range(len(thisJobIdList)):
            if(i == 0):
                thisDependency += 'done('+str(thisJobIdList[i])+')'
            else:
                thisDependency += ' && done('+str(thisJobIdList[i])+')'
        thisDependency += '\''
        thisBsub += thisDependency
        
    jt.nativeSpecification = thisBsub
    jt.remoteCommand = thisBundleDataScript
        
    jt.args = [thisOutPath, thisOutDirName, thisTempDir, str(I_readonly), str(I_compress), str(I_md5sum), str(I_debug)]
    jobid = s.runJob(jt)

def PrepareSampleProcessing(thisInFile,thisIdTableList,thisBatchNamesList,thisBatchCovarList,thisBatchVcfFileList, thisSnpinfoFileList, thisOutDir, thisTempDir):
    thisIdDictList = []
    for idTable in thisIdTableList:
        thisDict = MakeConcordanceTable(idTable)
        thisIdDictList.append(thisDict)
    
    # Check list to see which batch the biobank id belongs to
    (thisSampleToBatchList, thisNsampList) = CheckSampleArray(thisInFile,thisIdDictList)
        
    # Configure basenames for each batch (include # samples per batch)
    batchVariableList = []
    thisNSampCovarDict = {}
    for i in range(len(thisBatchNamesList)):
        thisBatchKey = thisBatchNamesList[i].split('_')[0]
        # batchVariableList consists of (snparray, covar)
        batchVariableList.append([thisBatchNamesList[i], thisBatchKey, thisBatchCovarList[i], thisBatchVcfFileList[i]])        
        
        if thisNSampCovarDict.has_key(thisBatchKey):
            thisNSampCovarDict[thisBatchKey][0] += thisNsampList[i]
        else:
            thisNSampCovarDict[thisBatchKey] = [thisNsampList[i]] # nSamp
            thisNSampCovarDict[thisBatchKey].append('') # baseName (e.g. Biobank_USERID_DATE_MEGA_xNN)
            thisNSampCovarDict[thisBatchKey].append('') # outPath
            thisNSampCovarDict[thisBatchKey].append('') # merge file list
            thisNSampCovarDict[thisBatchKey].append(thisSnpinfoFileList[i]) # snpinfo ref file
    
    #Now, create  entries in thisNSampCovarDict that require nSamp
    for thisKey in thisNSampCovarDict.keys():
    #     thisBatch = thisNSampCovarDict[thisKey][0]
        thisPrefix = '_'.join(baseName.split('_')[:-2]) + '_' + thisKey # remove _xNN_VCF
        thisNSamp = thisNSampCovarDict[thisKey][0]
        thisNSampCovarDict[thisKey][1] = thisPrefix+'_x' + str(thisNSamp)
        thisNSampCovarDict[thisKey][2] = thisOutDir+'/'+thisNSampCovarDict[thisKey][1]
        thisNSampCovarDict[thisKey][3] = thisTempDir+'/'+thisNSampCovarDict[thisKey][1]+'.mergelist.txt'
        
    return thisIdDictList, thisSampleToBatchList, thisNSampCovarDict, batchVariableList



####################################################################################
################################   Argument Parser   ###############################
####################################################################################
parser = ArgumentParser()

parser.add_argument("-i","--in",
                    type=str,
                    dest="infile",
                    help="csv-version of requested samples")

parser.add_argument("-o","--out",
                    type=str,
                    dest="outdir",
                    help="output directory for data")

parser.add_argument("--idtable",
                    nargs='+',
                    type=str,
                    dest="idtablelist",
                    help="List of concordance tables between biobank ID and Illumina MEGA array sentrix ID")

parser.add_argument("--batchvcffile",
                    nargs='+',
                    type=str,
                    dest="batchvcffilelist",
                    help="List of batch VCF files to extract subset datasets from")

parser.add_argument("--snpinfo",
                    nargs='+',
                    type=str,
                    dest="snpinfofilelist",
                    help="corresponding snpinfo files for each array")

parser.add_argument("--batchnames",
                    nargs='+',
                    type=str,
                    dest="batchnameslist",
                    help="Array of Illumina biobank array batch names")

parser.add_argument("--batchcovar",
                    nargs='+',
                    type=str,
                    dest="batchcovarlist",
                    help="Array of covariates matching the above batch names")

parser.add_argument("--readme",
                    type=str,
                    dest="readmefile",
                    help="corresponding README file for the final vcf files (MEGA and MEGAEX arrays)")

parser.add_argument("-y","--yaml",
                    type=str,
                    dest="yamlfile",
                    required=True,
                    help="corresponding YAML file for the final ped files")

parser.add_argument("-5","--md5sum",
                    action="store_true",
                    dest="I_md5sum",
                    default=False,
                    help="flag for md5sum creation")

parser.add_argument("--readonly",
                    action="store_true",
                    dest="I_readonly",
                    default=False,
                    help="flag for marking files read-only")

parser.add_argument("-c","--compress",
                    action="store_true",
                    dest="I_compress",
                    default=False,
                    help="flag for compression")

parser.add_argument("-b","--bin",
                    type=str,
                    dest="bindir",
                    help="biobank bin directory")

parser.add_argument("-d","--debug",
                    action="store_true",
                    dest="I_debug",
                    default=False,
                    help="Flag if you want to keep temp directory for debugging purposes")
 
args = parser.parse_args()


####################################################################################
################################    Main Function    ###############################
####################################################################################
# Check between input from config file and commandline parameter. Use the parameter from command line to overwrite 
# any variables assigned through the config file

yamlFile = args.yamlfile

with open(yamlFile, 'r') as stream:
    yml = yaml.load(stream)

##################
### Files
##################
# Assign Infile
if(args.infile != None):
    inFile = args.infile
else:
    inFile = AssignVariableFromConfig(yml, False, 'Files', 'Infile', "infile")
    if(inFile == None):
        sys.exit("Need to specify an inFile")

## define baseName
baseName = os.path.basename(inFile).replace('.csv','')

# Assign Readme
if(args.readmefile != None):
    readmeFile = args.readmefile
else:
    readmeFile = AssignVariableFromConfig(yml, False, 'Files', 'MergedReadme', "readme")

# Assign idTable
if(args.idtablelist != None):
    idTableList = args.idtablelist
else:
    idTableList = AssignVariableFromConfig(yml, True, 'Files', 'ConcordanceTable', "concordance tables")

# Assign batchVCFFile
if(args.batchvcffilelist != None):
    batchVcfFileList = args.batchvcffilelist
else:
    batchVcfFileList = AssignVariableFromConfig(yml, True, 'Files', 'BatchVcfFiles', "batch vcf files")

# Assign snpinfofile
if(args.snpinfofilelist != None):
    snpinfoFileList = args.snpinfofilelist
else:
    snpinfoFileList = AssignVariableFromConfig(yml, True, 'Files', 'SnpInfo', "snpinfo files")


##################
### Directories
##################
# Assign basedir
baseDir = AssignVariableFromConfig(yml, False, 'Directories', 'BaseDir', "basedir")

# Assign outdir
if(args.outdir != None):
    outDir = args.outdir
else:
    outDir = AssignVariableFromConfig(yml, False, 'Directories', 'OutDir', "outdir")
    if(outDir == None):
        if(baseName == ''):
            print "something wrong, no basename"
        elif(baseDir != None):
            outDir = SetDefaultOutDir(baseDir,baseName)
        else:
            outDir = SetDefaultOutDir('',baseName)

# Assign biobank bin dir
if(args.bindir != None):
    binDir = args.bindir
else:
    binDir = AssignVariableFromConfig(yml, False, 'Directories', 'BinDir', "bindir")
    
extractVcfScript = binDir+'/extract_vcf.sh'
mergeVcfScript = binDir+'/merge_vcf.sh'
bundleDataScript = binDir+'/bundle_data.sh'

# Assign tempDir/tempdir
try:
    tempDir = yml['Directories']['TempDir']
    if(tempDir == None):
        tempDir = outDir
except:
    tempDir = outDir
# create tempdir using basename
tempDir = tempDir+'/'+baseName
if not os.path.isdir(tempDir):
    os.makedirs(tempDir)
else:   # currently we are just going to continue saving to this directory. We can decide to throw an error if required alter
#     sys.exit('tempDir may not be empty. Please delete this directory and try again: '+tempDir)
    pass

##################
### Other Settings
##################
# Assign md5sum flag
if(args.I_md5sum == True):
    I_md5sum = args.I_md5sum
else:
    I_md5sum = AssignVariableFromConfig(yml, False, 'Settings', 'Md5sum', "I_md5sum")

# Assign read only flag
if(args.I_readonly == True):
    I_readonly = args.I_readonly
else:
    I_readonly = AssignVariableFromConfig(yml, False, 'Settings', 'ReadOnlyPermissions', "I_readonly")

# Assign compression flag
if(args.I_compress == True):
    I_compress = args.I_compress
else:
    I_compress = AssignVariableFromConfig(yml, False, 'Settings', 'CompressFiles', "I_compress")

# Assign batch names list
if(args.batchnameslist != None):
    batchNamesList = args.batchnameslist
else:
    batchNamesList = AssignVariableFromConfig(yml, True, 'Settings', 'ArrayBatchNames', "batch names")

# Assign batch covariate list
if(args.batchcovarlist != None):
    batchCovarList = args.batchcovarlist
else:
    batchCovarList = AssignVariableFromConfig(yml, True, 'Settings', 'ArrayBatchCovariates', "batch covariate names")

# Assign debug flag
I_debug = args.I_debug

######################### End of setting config variables, start of official script
# Prepare for PED processing
(thisIdDictList, thisSampleToBatchList, thisNSampCovarDict, batchVariableList) = PrepareSampleProcessing(inFile,idTableList,batchNamesList,batchCovarList,batchVcfFileList,snpinfoFileList,outDir,tempDir)

s = drmaa.Session()
s.initialize()

# Make VCF files for designated samples
thisJobIdList = ExtractVcfDataset(s, thisSampleToBatchList, thisNSampCovarDict, thisIdDictList, outDir, tempDir, batchNamesList, batchVariableList, extractVcfScript)
# Need to merge and move final VCF files
thisMergeJobIdList = MergeVcfFiles(s, thisJobIdList, thisNSampCovarDict, batchVariableList, tempDir, mergeVcfScript)

# Make one merged readme file
MakeMergedReadmeFile(readmeFile, outDir, thisNSampCovarDict)

# Bundle to shared location
BundleToDropbox(s, outDir, tempDir, I_readonly, I_compress, I_md5sum, I_debug, thisMergeJobIdList, bundleDataScript)


