#!/usr/bin/env python

# Purpose:       Get individual imputed genotypes 
# 
# Usage:         python GetIndividualImputedGenotypesMultiArray.py ...
# 
# Author:        et85, etsai@bwh.harvard.edu


from argparse import ArgumentParser
import drmaa, os, re, shutil, subprocess, sys, time, yaml

####################################################################################
#############################   Function Definitions   #############################
####################################################################################
def AssignVariableFromConfig(thisYaml, var1, var2, desc):
    try:
        thisVariable = thisYaml[var1][var2]
        
        # if this is a data file, please make sure it exists
        if(var1 == 'DataFiles'):
            if(os.path.isfile(thisVariable)):
                return thisVariable
            else:
                print "Invalid " + desc + " path: " + thisVariable
                sys.exit()
#         elif(var1 == 'Directories'):        # don't assign this here,
#             if(os.path.exists(thisVariable)):
#                 return thisVariable
        elif(desc.startswith("I_")):
            return thisVariable
        else:
            return thisVariable
    except:
        sys.exit("Please assign " + desc + " by adding this parameter or editing the config file")

def FormatDate(thisBaseName):
    thisUnformDate = thisBaseName.split("_")[2]
    thisFormDate = thisUnformDate[0:4] + '-' + thisUnformDate[4:6] + '-' + thisUnformDate[6:]
    return thisFormDate

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

## start editing here
def MakePedFile(thisIdList, thisIdDict, thisPedPath, thisOutPath, thisBaseName, thisFullPed, I_allsamp):
    if not os.path.exists(thisOutPath):
        os.makedirs(thisOutPath)
    
    if(I_allsamp):    #  just copy the full PED file over
        shutil.copy(thisFullPed,thisOutPath + '/' + thisBaseName + '.ped')
    
    else:   # go through the CSV file one by one
        thisOutFH = open(thisOutPath + '/' + thisBaseName + '.ped','w')
        for thisBiobankId in thisIdList:
            if(thisIdDict.has_key(thisBiobankId)):
                thisSentrixId = thisIdDict[thisBiobankId]
                print thisBiobankId, thisSentrixId
                
                with open(thisPedPath + '/' + thisSentrixId + '.ped', 'r') as thisPedFH:
                    thisOutFH.write(thisPedFH.read())
                
            else:
                print "Cannot find key."
        thisOutFH.close()
    
def CopySnpFiles(thisMapFile, thisSnpInfoFile, thisOutPath, thisBaseName):    
    # check to make sure MAP file exists in the correct directory
    if not os.path.isfile(thisMapFile):
        sys.exit("MAP file does not exist. Please check file path: " + thisMapFile)
    if not os.path.isfile(thisSnpInfoFile):
        sys.exit("SNP Info file does not exist. Please check file path: " + thisSnpInfoFile)
    
    shutil.copy(thisMapFile,thisOutPath + '/' + thisBaseName + '.map')
    shutil.copy(thisSnpInfoFile,thisOutPath + '/' + thisBaseName + '.snpinfo.txt')
    
def MakeMergedReadmeFile(thisReadmeFile, thisOutPath, thisBaseName1, thisBaseName2):
    if not os.path.isfile(thisReadmeFile):
        sys.exit("README file does not exist. Please check file path: " + thisReadmeFile)
    
    thisFormDate = time.strftime('%Y-%m-%d')
    thisOutFH = open(thisOutPath + '/README.txt',"w")
    
    with open(thisReadmeFile, "r") as thisReadmeFH:
        for thisLine in thisReadmeFH.readlines():
            thisLine = thisLine.replace("INSERTFILENAME_MEGAEX", thisBaseName2)
            thisLine = thisLine.replace("INSERTFILENAME_MEGA", thisBaseName1)
            thisLine = thisLine.replace("INSERTDATE", thisFormDate)
            thisOutFH.write(thisLine)
    
    thisOutFH.close()

## stop editing here

def ExtractDataset(thisIdList, thisIdDict, thisImputeDir, thisOutPath, thisTempDir, thisBaseName, I_allsamp, thisCopyChrScript, thisExtractScript):
    
#     thisFH = open(thisInfile, "r")
        
    # Make output dir if it doesn't exist.
    if os.path.exists(thisOutPath):
        # temporarily disable check for empty outdir
        pass
#         if(os.listdir(thisOutPath) != []):
#             sys.exit("This out dir is not empty. Please delete the contents in the out dir before trying again:\n" + thisOutPath)
    else:
        os.makedirs(thisOutPath)
    
    
    # Initialize array
    thisSampleList = []
    
    for thisBiobankID in thisIdList:
        if(thisIdDict.has_key(thisBiobankID)):
            thisSentrixID = thisIdDict[thisBiobankID]
            thisSampleList.append(thisSentrixID + '-' + thisBiobankID)
            print thisSentrixID + '-' + thisBiobankID
        else:
            print "Cannot find key."
    
    # Make temporary outfile for sample names in ascii
#     thisSampListFile = 'biobank_x'+str(len(thisSampleList))
#     thisSampListFH = open(thisTempDir+'/'+thisSampListFile + '.txt','w')
    thisSampListFH = open(thisTempDir+'/'+thisBaseName + '.txt','w')

    for thisSamp in thisSampleList:
        thisSampListFH.write(thisSamp + '\n')
    thisSampListFH.close()
    
    s = drmaa.Session()
    s.initialize()
    print('creating job template for extract_by_chr')
    jt = s.createJobTemplate()
    if(I_allsamp):
        # TODO: Copy all the files over if there is just one data set
        thisBsubOut = thisTempDir + '/copy_by_chr_%J_%I.out'
        thisBsubErr = thisTempDir + '/extract_by_chr_%J_%I.err'
        jt.nativeSpecification = '-R \'rusage[mem=8192]\' ' + \
                                 '-M 12288 ' + \
                                 '-q pcpgmwgs ' + \
                                 '-J test_array ' + \
                                 '-o ' + thisBsubOut + ' ' \
                                 '-e ' + thisBsubErr
        jt.remoteCommand = thisCopyChrScript
        jt.args = [thisBaseName,thisImputeDir,thisOutPath]

        jobid = s.runBulkJobs(jt, 1, 22, 1)
        
    else:
        thisBsubOut = thisTempDir + '/extract_by_chr_%J_%I.out'
        thisBsubErr = thisTempDir + '/extract_by_chr_%J_%I.err'
        jt.nativeSpecification = '-R \'rusage[mem=8192]\' ' + \
                                 '-M 12288 ' + \
                                 '-q pcpgmwgs ' + \
                                 '-J test_array ' + \
                                 '-o ' + thisBsubOut + ' ' \
                                 '-e ' + thisBsubErr
        
        jt.remoteCommand = thisExtractScript
#         jt.args = [thisSampListFile,thisImputeDir,thisOutPath,thisTempDir]
        jt.args = [thisBaseName,thisImputeDir,thisOutPath,thisTempDir]
        jobid = s.runBulkJobs(jt, 1, 22, 1)
        
    print('extract_by_chr has been submitted with ID %s' % jobid)
    print('cleaning up')
    s.deleteJobTemplate(jt)
    
    thisMainJobID = int(str(jobid[0]).split('[')[0])
    
    s.exit()
    
#     return thisMainJobID, thisSampListFile
    return thisMainJobID, thisBaseName

    
def CopyInfoFiles(thisOutPath, thisImputeDir, thisTempDir, thisSampleListName, thisCopyInfoScript):
    s = drmaa.Session()
    s.initialize()
    print('creating job template for copy_info_files')
    jt = s.createJobTemplate()
    thisBsubOut = thisTempDir + '/copy_info_files_%J.out'
    thisBsubErr = thisTempDir + '/copy_info_files_%J.err'
    jt.nativeSpecification = '-R \'rusage[mem=4096]\' ' + \
                             '-M 8192 ' + \
                             '-q pcpgmwgs ' + \
                             '-J copy_info_files ' + \
                             '-o ' + thisBsubOut + ' ' \
                             '-e ' + thisBsubErr
    
    jt.remoteCommand = thisCopyInfoScript
    
    jt.args = [thisOutPath,thisImputeDir,thisSampleListName]
    jobid = s.runJob(jt)
    
    print('copy_info_files has been submitted with id %s' % jobid)
    print('cleaning up')
    s.deleteJobTemplate(jt)
    
    s.exit()
    
    return jobid
    
def BundleToDropbox(thisOutPath, thisBaseName, thisSharePath, thisTempDir, I_compress, I_readonly, I_md5sum, thisUsernameList,thisExtractJobID1,thisCopyJobID1,thisExtractJobID2,thisCopyJobID2):
    
    thisOutDirPrefix = "/".join(thisOutPath.split("/")[:-1])  # keep directory name up to the outer most directory
    thisOutPathBase = thisOutPath.split("/")[-1]

    
    # Make a bsub script that does the file transfers after the previous commands are successfully executed
    thisDropboxScriptFH = open(thisTempDir+'/transfer_files.sh','w')
    thisDropboxScriptFH.write('#!/bin/bash\n\n')
    
    # Set files to read-only (optional)
    if(I_readonly):
        thisDropboxScriptFH.write('chmod 440 '+thisOutPath+'/*\n')    
        thisDropboxScriptFH.write('chmod 440 '+thisTempDir+'/README.txt\n')
    
    # Compress files (optional)
    if(I_compress):
        thisDropboxScriptFH.write('chmod 550 '+thisOutPath+'\n')
        thisDropboxScriptFH.write('tar czvf '+thisOutPath+'.tgz '+'-C '+thisOutDirPrefix+' '+thisOutPathBase+'/\n')
        
        if(I_md5sum):
            thisDropboxScriptFH.write('md5sum '+thisOutPath+' > '+thisOutPath + '.tgz.md5\n')
            
    else: # put the readme in the directory
        thisDropboxScriptFH.write('mv '+thisTempDir+'/README.txt '+thisOutPath+'/README.txt\n')
        thisDropboxScriptFH.write('chmod 550 '+thisOutPath+'\n')
    
    if not os.path.exists(thisSharePath):
        thisDropboxScriptFH.write('mkdir '+thisSharePath+'\n')
    else:
        print 'This dir already exists: ' + thisSharePath
    
    # Copy the files over to the shared directory
    #@ put README into the directory if it's not compressed
    if(I_compress):
        # Check to make sure that these files do not already exist
        if(os.path.exists(thisSharePath + '/' + thisBaseName + '.tgz')):
            print thisSharePath + '/' + thisBaseName + '.tgz' + " already exists. No replacement attempted"
        else:
            thisDropboxScriptFH.write('cp -pr '+thisOutPath+'.tgz '+thisSharePath+'\n')
            thisDropboxScriptFH.write('cp -p '+thisTempDir+'/README.txt '+thisSharePath+'/README.txt\n')
            thisDropboxScriptFH.write('chmod 660 '+thisTempDir+'/README.txt\n')
            
            if(I_md5sum):
                thisDropboxScriptFH.write('cp -pr '+thisOutPath+'.tgz.md5 '+thisSharePath+'\n')
                if(I_readonly):
                    thisDropboxScriptFH.write('chmod 440 '+thisSharePath+'/'+thisBaseName+'.tgz.md5\n')
         
        if(I_readonly):
            # remove the directory
            thisDropboxScriptFH.write('chmod 440 '+thisSharePath+'/'+thisBaseName+'.tgz\n')
            thisDropboxScriptFH.write('chmod 700 '+thisOutPath+'\n')
         
        # remove the directory file
        thisDropboxScriptFH.write('rm -rf '+thisOutPath+'\n')
    else:
        if(os.path.exists(thisSharePath + '/' + thisBaseName)):
            print thisSharePath + '/' + thisBaseName + " already exists. No replacement attempted"
        else:
            thisDropboxScriptFH.write('cp -pr '+thisOutPath+' '+thisSharePath+'\n')
     
    if(thisUsernameList != None):
        for thisName in thisUsernameList:
            thisDropboxScriptFH.write('setfacl -R -m u:'+thisName+':rwX '+thisSharePath+'\n')
                
    # remove the temp directory
    thisDropboxScriptFH.write('rm -rf '+thisTempDir+'\n')
    
    thisDropboxScriptFH.close()
    os.chmod(thisTempDir+'/transfer_files.sh', 0770)
    
    s = drmaa.Session()
    s.initialize()
    print('creating job template for copy_info_files')
    jt = s.createJobTemplate()
    thisBsubOut = thisTempDir + '/transfer_files_%J.out'
    thisBsubErr = thisTempDir + '/transfer_files_%J.err'
    jt.nativeSpecification = '-R \'rusage[mem=4096]\' ' + \
                             '-M 8192 ' + \
                             '-q pcpgmwgs ' + \
                             '-J transfer_files ' + \
                             '-w \'done('+str(thisExtractJobID1)+') && done('+str(thisCopyJobID1) + \
                             ') && done('+str(thisExtractJobID2)+') && done('+str(thisCopyJobID2)+')\' ' + \
                             '-o ' + thisBsubOut + ' ' \
                             '-e ' + thisBsubErr
    
    jt.remoteCommand = thisTempDir+'/transfer_files.sh'
    
    jobid = s.runJob(jt)
    
    print('transfer_files has been submitted with id %s' % jobid)
    print('cleaning up')
    s.deleteJobTemplate(jt)
    s.exit()
 
# def BundleToDropbox(thisOutPath, thisBaseName, thisBaseName1, thisBaseName2, thisSharePath, I_readonly, I_compress, I_md5sum, thisUsernameList):
#     thisRead = thisOutPath.split("/")
#     thisOutDirName = thisRead[-1]
#     
#     if(I_readonly):
#         subprocess.call(['chmod','440',thisOutPath + '/' + thisBaseName1 + '.ped'])
#         subprocess.call(['chmod','440',thisOutPath + '/' + thisBaseName1 + '.map'])
#         subprocess.call(['chmod','440',thisOutPath + '/' + thisBaseName1 + '.snpinfo.txt'])
#         subprocess.call(['chmod','440',thisOutPath + '/' + thisBaseName2 + '.ped'])
#         subprocess.call(['chmod','440',thisOutPath + '/' + thisBaseName2 + '.map'])
#         subprocess.call(['chmod','440',thisOutPath + '/' + thisBaseName2 + '.snpinfo.txt'])
#         subprocess.call(['chmod','440',thisOutPath + '/README.txt'])
#         subprocess.call(['chmod', '550', thisOutPath])
#     
#     if(I_compress):
#         subprocess.call(['tar', 'czvf', thisOutPath + '.tgz', '-C', thisOutPath + '/../', thisOutDirName])
#         if(I_md5sum):
#             thisMd5FH = open(thisOutPath + '.tgz.md5', 'w')
#             subprocess.call(['md5sum', thisOutPath + '.tgz'], stdout=thisMd5FH)
#             thisMd5FH.close()
#     
#     if not os.path.exists(thisSharePath):
#         subprocess.call(['mkdir', thisSharePath])
#     else:
#         print 'This dir already exists: ' + thisSharePath
#     
#     if(I_compress):
#         # Check to make sure that these files do not already exist
#         if(os.path.exists(thisSharePath + '/' + thisBaseName + '.tgz')):
#             print thisSharePath + '/' + thisBaseName + '.tgz' + " already exists. No replacement attempted"
#         else:
#             subprocess.call(['cp', '-pr', thisOutPath + '.tgz', thisSharePath])
#             if(I_md5sum):
#                 subprocess.call(['cp', '-pr', thisOutPath + '.tgz.md5', thisSharePath])
#         
#         if(I_readonly):
#             # remove the directory
#             subprocess.call(['chmod', '700', thisOutPath])
#         
#         # remove the directory file
#         subprocess.call(['rm', '-rf', thisOutPath])
#     else:
#         if(os.path.exists(thisSharePath + '/' + thisBaseName)):
#             print thisSharePath + '/' + thisBaseName + " already exists. No replacement attempted"
#         else:
#             subprocess.call(['cp', '-pr', thisOutPath, thisSharePath])
#     
#     if(thisUsernameList != None):
#         for thisName in thisUsernameList:
#             subprocess.call(['setfacl', '-R', '-m', 'u:' + thisName + ':rwX', thisSharePath])

def setDefaultOutDir(thisBaseDir,thisBaseName):
    if (thisBaseDir == ''):
        thisBaseDir = '/data/ppm/biobank/release'
    
    thisUsername = FormatUsername(thisBaseName)

    #Check to make sure you're on panasas
    if(os.path.exists(thisBaseDir)):
        subprocess.call(['mkdir', '-p', thisBaseDir + '/' + thisUsername])
    else:
        sys.exit('Cannot find path: '+thisBaseDir)
    return thisBaseDir + '/' + thisUsername + '/' + thisBaseName
    
def setDefaultTransferDir(thisBaseName):
    thisSharePath = '/pub/dropbox/' + time.strftime('%Y-%m-%d') + '_for_' + FormatUsername(thisBaseName)
        
    #Check to make sure you're on panasas
    if(os.path.exists('/pub/dropbox/')):
        return thisSharePath
    else:
        sys.exit('Cannot find path: /pub/dropbox/')

def checkSampleArray(thisFile,thisIdDict1,thisIdDict2):
    thisV1List = []
    thisV2List = []
    
    thisFH = open(thisFile,'r')
    # Decode Byte Order Mark
    s = thisFH.read()
    u = s.decode("utf-16")
    s = u.encode("utf-8")
    s = s.rstrip()
    thisIdList = s.split('\r\n')[2:]     # skip the first two header rows
    for thisBiobankID in thisIdList:
        if(thisIdDict1.has_key(thisBiobankID)):
            thisV1List.append(thisBiobankID)
        elif(thisIdDict2.has_key(thisBiobankID)):
            thisV2List.append(thisBiobankID)
        else:
            print thisBiobankID+': Cannot find sentrix id for this subject id'
    return (thisV1List, thisV2List)

####################################################################################
################################   Argument Parser   ###############################
####################################################################################
parser = ArgumentParser()

parser.add_argument("-i","--in",
                    type=str,
                    dest="infile",
                    help="csv-version of requested samples")

parser.add_argument("--impute1",
                    type=str,
                    dest="imputedir1",
                    help="directory containing imputation data for the entire dataset")

parser.add_argument("--impute2",
                    type=str,
                    dest="imputedir2",
                    help="directory containing imputation data for the entire dataset")

parser.add_argument("-o","--out",
                    type=str,
                    dest="outdir",
                    help="output directory for data")

parser.add_argument("-t","--transfer",
                    type=str,
                    dest="transferdir",
                    help="transfer directory for data")

parser.add_argument("--idtable1",
                    type=str,
                    dest="idtable1",
                    help="Concordance table between biobank ID and Illumina MEGA array sentrix ID")

parser.add_argument("--idtable2",
                    type=str,
                    dest="idtable2",
                    help="Concordance table between biobank ID and Illumina MEGAEX array sentrix ID")

parser.add_argument("--readme1",
                    type=str,
                    dest="readmefile1",
                    help="corresponding README file for the final vcf files of MEGA array")

parser.add_argument("--readme2",
                    type=str,
                    dest="readmefile2",
                    help="corresponding README file for the final vcf files of MEGAEX array")

parser.add_argument("--readme",
                    type=str,
                    dest="readmefile",
                    help="corresponding README file for the final vcf files (MEGA and MEGAEX arrays)")

parser.add_argument("-y","--yaml",
                    type=str,
                    dest="yamlfile",
                    required=True,
                    help="corresponding YAML file for the final ped files")

parser.add_argument("-c","--compress",
                    action="store_true",
                    dest="I_compress",
                    help="flag for compression and deposit to dropbox")
 
parser.add_argument("-5","--md5sum",
                    action="store_true",
                    dest="I_md5sum",
                    help="flag for md5sum creation")

parser.add_argument("--yaml1",
                    type=str,
                    dest="yaml1",
                    help="config file for MEGA array")

parser.add_argument("--yaml2",
                    type=str,
                    dest="yaml2",
                    help="config file for MEGAEX array")

parser.add_argument("--allsamp1",
                    action="store_true",
                    dest="I_allsamp1",
                    help="set string for delivering the designated fullped file for MEGA array")

parser.add_argument("--allsamp2",
                    action="store_true",
                    dest="I_allsamp2",
                    help="set string for delivering the designated fullped file for MEGAEX array")

parser.add_argument("-u","--usernamelist",
                    type=str,
                    dest="usernamelist",
                    help="designate a list of partners userids to allow access, comma-separated, no spaces!")

parser.add_argument("--bin",
                    type=str,
                    dest="bindir",
                    help="designate the path to the biobank bin directory")
 
args = parser.parse_args()


####################################################################################
################################    Main Function    ###############################
####################################################################################

# Check between input from config file and commandline parameter. Use the parameter from command line to overwrite 
# any variables assigned through the config file

yamlFile = args.yamlfile

with open(yamlFile, 'r') as stream:
    yml = yaml.load(stream)

# Assign idTable
if(args.idtable1 != None):
    idTable1 = args.idtable1
else:
    idTable1 = AssignVariableFromConfig(yml, 'DataFiles', 'ConcordanceTable1', "concordance table1")
if(args.idtable2 != None):
    idTable2 = args.idtable2
else:
    idTable2 = AssignVariableFromConfig(yml, 'DataFiles', 'ConcordanceTable2', "concordance table2")
 
# Assign Readme
if(args.readmefile != None):
    readmeFile = args.readmefile
else:
    readmeFile = AssignVariableFromConfig(yml, 'DataFiles', 'MergedReadme', "readme")

if(args.readmefile1 != None):
    readmeFile1 = args.readmefile1
else:
    readmeFile1 = AssignVariableFromConfig(yml, 'DataFiles', 'Readme1', "readme1")
if(args.readmefile2 != None):
    readmeFile2 = args.readmefile2
else:
    readmeFile2 = AssignVariableFromConfig(yml, 'DataFiles', 'Readme2', "readme2")

# Assign Individual Yaml Config file
if(args.yaml1 != None):
    yamlFile1 = args.yaml1
else:
    yamlFile1 = AssignVariableFromConfig(yml, 'DataFiles', 'ConfigFile1', "yaml1")
if(args.yaml2 != None):
    yamlFile2 = args.yaml2
else:
    yamlFile2 = AssignVariableFromConfig(yml, 'DataFiles', 'ConfigFile2', "yaml2")

# Deliver FullPed file
if(args.I_allsamp1 != None):
    I_allsamp1 = args.I_allsamp1
else:
    I_allsamp1 = AssignVariableFromConfig(yml, 'Settings', 'FullDataset1', "I_allsamp")
if(args.I_allsamp2 != None):
    I_allsamp2 = args.I_allsamp2
else:
    I_allsamp2 = AssignVariableFromConfig(yml, 'Settings', 'FullDataset2', "I_allsamp")

################
 
 
# Directories

# Assigning impute directory
if(args.imputedir1 != None):
    imputeDir1 = args.imputedir1
else:
    imputeDir1 = AssignVariableFromConfig(yml, 'Directories', 'ImputeDir1', "impute dir1")

if(args.imputedir2 != None):
    imputeDir2 = args.imputedir2
else:
    imputeDir2 = AssignVariableFromConfig(yml, 'Directories', 'ImputeDir2', "impute dir2")

# infile
if(args.infile != None):
    inFile = args.infile
else:
    inFile = AssignVariableFromConfig(yml, 'Directories', 'InFile', "infile")
    if(inFile == None):
        sys.exit("Need to specify an inFile")
 
## define baseName
baseName = os.path.basename(inFile).replace('.csv','')


# @TODO: basedir is currently unused. Not sure if we need this?
baseDir = AssignVariableFromConfig(yml, 'Directories', 'BaseDir', "basedir")

# outfile
if(args.outdir != None):
    outDir = args.outdir
else:
    outDir = AssignVariableFromConfig(yml, 'Directories', 'OutDir', "outdir")
    if(outDir == None):
        if(baseName == ''):
            print "something wrong, no basename"
        elif(baseDir != None):
            outDir = setDefaultOutDir(baseDir,baseName)
        else:
            outDir = setDefaultOutDir('',baseName)


# transferdir
if(args.transferdir != None):
    transferDir = args.transferdir
else:
    try:
        transferDir = yml['Directories']['TransferDir']
        if(transferDir == None):
            transferDir = setDefaultTransferDir(baseName)
    except:
        transferDir = setDefaultTransferDir(baseName)

# temppath/tempdir
try:
    tempPath = yml['Directories']['TempPath']
    if(tempPath == None):
        tempPath = outDir
except:
    tempPath = outDir
# create tempdir using basename
tempDir = tempPath+'/'+baseName
if not os.path.isdir(tempDir):
    os.makedirs(tempDir)
else:   # currently we are just going to continue saving to this directory. We can decide to throw an error if required alter
#     sys.exit('tempDir may not be empty. Please delete this directory and try again: '+tempDir)
    pass


### Other Settings
if(args.I_compress != None):
    I_compress = args.I_compress
else:
    I_compress = AssignVariableFromConfig(yml, 'Settings', 'CompressDir', "I_compress")
 
if(args.I_md5sum != None):
    I_md5sum = args.I_md5sum
else:
    I_md5sum = AssignVariableFromConfig(yml, 'Settings', 'Md5sum', "I_md5sum")

if(args.bindir != None):
    binDir = args.bindir
else:
    binDir = AssignVariableFromConfig(yml, 'Settings', 'BiobankBinDir', "biobank bin dir")



# @TODO: currently, the stuff on the bottom don't have commandline parameters
# ReadOnlyPermissions
I_readonly = AssignVariableFromConfig(yml, 'Settings', 'ReadOnlyPermissions', "I_readonly")
 
# Move to dropbox
I_dropbox = AssignVariableFromConfig(yml, 'Settings', 'MoveToDropbox', "I_dropbox")
 
# username list
if(args.usernamelist != None):
    username = args.usernamelist
    usernameList = args.usernamelist.split(',')
else:
    username = AssignVariableFromConfig(yml, 'Settings', 'AllowUsernames', "usernameList")
    if(username != None):
        usernameList = username.split(',')
    else:
        usernameList = None

###  Templates for bsub
extractchrScript = AssignVariableFromConfig(yml, 'Scripts', 'ExtractByChr', "extractchrScript")
copychrScript = binDir + '/copy_by_chr.sh'
copyinfoScript = AssignVariableFromConfig(yml, 'Scripts', 'CopyInfoFiles', "copyinfoScript")

# check to make sure these scripts are in place
if(extractchrScript == None):
    sys.exit('Please assign location of ExtractByChr script in the config file\n')
elif(copyinfoScript == None):
    sys.exit('Please assign location of CopyInfoFiles script in the config file\n')
else:   # check to make sure this script actually exists    
    if (not os.path.isfile(extractchrScript)):
        sys.exit('Please make sure location of ExtractByChr script is correct: ' + extractchrScript + '\n')
    elif (not os.path.isfile(copyinfoScript)):
        sys.exit('Please make sure location of CopyInfoFiles script is correct: ' + copyinfoScript + '\n')



#### Make sure the directories exist

# make sure outDir exists
if not os.path.isdir(outDir):
    # check to see if the directory right above it exists
    if not os.path.isdir('/'.join(outDir.split("/")[:-1])):
        sys.exit("Cannot find path to out dir: " + outDir)
    else:
        # make directory if only that first 
        os.makedirs(outDir)
 
# set defaults for I_readonly, I_dropbox
if(I_readonly == None):
    I_readonly = False
if(I_dropbox == None):
    I_dropbox = False
 

######################### End of setting config variables, start of official script
IdDict1 = MakeConcordanceTable(idTable1)
IdDict2 = MakeConcordanceTable(idTable2)

(thisV1List, thisV2List) = checkSampleArray(inFile,IdDict1,IdDict2)

if(thisV1List == []):   # only v2 samples requested, call script to generate this
    thisCmd = ['python',binDir+'/GetImputedGenotypes.py','-i',inFile,'-o',outDir,'-y',yamlFile2,
                    '--idtable',idTable2,'--impute',imputeDir2,
                    '--readme',readmeFile2,'--transfer',transferDir]
    if(I_compress):
        thisCmd.append('--compress')
    if(I_md5sum):
        thisCmd.append('--md5sum')
    if(I_allsamp2):
        thisCmd.append('--allsamp')
    if(username != None):
        thisCmd.append('--usernamelist')
        thisCmd.append(username)
    
    print "MEGAEX Array samples only. Calling the following command:"
    print " ".join(thisCmd)
    subprocess.call(thisCmd)
elif(thisV2List == []): # only v1 samples requested, call script to generate this
    
    thisCmd = ['python',binDir+'/GetImputedGenotypes.py','-i',inFile,'-o',outDir,'-y',yamlFile1,
                    '--idtable',idTable1,'--impute',imputeDir1,
                    '--readme',readmeFile1,'--transfer',transferDir]
    
    if(I_compress):
        thisCmd.append('--compress')
    if(I_md5sum):
        thisCmd.append('--md5sum')
    if(I_allsamp1):
        thisCmd.append('--allsamp')
    if(username != None):
        thisCmd.append('--usernamelist')
        thisCmd.append(username)
    
    print "MEGA Array samples only. Calling the following command:"
    print " ".join(thisCmd)
    subprocess.call(thisCmd)

### These values are not assigned
#     baseDir
#     I_dropbox
#     I_readonly
#     binDir
    
else:   # mix of v1 and v2 samples
    # Create new BaseNames
    baseName1 = '_'.join(baseName.split('_')[:-1]) + '_MEGA_x' + str(len(thisV1List))
    baseName2 = '_'.join(baseName.split('_')[:-1]) + '_MEGAEX_x' + str(len(thisV2List))
    
    # Make impute directories for MEGA and MEGAEX
    (thisExtractJobID1,thisSampleListName1) = ExtractDataset(thisV1List, IdDict1, imputeDir1, outDir, tempDir, baseName1, I_allsamp1, copychrScript, extractchrScript)
    thisCopyJobID1 = CopyInfoFiles(outDir, imputeDir1, tempDir, thisSampleListName1, copyinfoScript)
    (thisExtractJobID2,thisSampleListName2) = ExtractDataset(thisV2List, IdDict2, imputeDir2, outDir, tempDir, baseName2, I_allsamp2, copychrScript, extractchrScript)
    thisCopyJobID2 = CopyInfoFiles(outDir, imputeDir2, tempDir, thisSampleListName2, copyinfoScript)
    
    # Make one merged readme file
    MakeMergedReadmeFile(readmeFile, outDir, baseName1, baseName2)
    
    if(I_dropbox):
#         BundleToDropbox(thisOutPath, thisBaseName, thisSharePath, thisTempDir, I_compress, I_readonly, I_md5sum, thisUsernameList,thisExtractJobID,thisCopyJobID):
        BundleToDropbox(outDir, baseName, transferDir, tempDir, I_compress, I_readonly, I_md5sum, usernameList,thisExtractJobID1,thisCopyJobID1,thisExtractJobID2,thisCopyJobID2)

#         BundleToDropbox(outDir, baseName, baseName1, baseName2, transferDir, I_readonly, I_compress, I_md5sum, usernameList)
