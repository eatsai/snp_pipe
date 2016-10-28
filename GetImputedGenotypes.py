#!/usr/bin/env python

# Purpose:       Get the imputed genotypes for a selected number of individuals
# 
# Usage:         python GetImputedGenotypes.py ...
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

def MakeConcordanceTable(idtablepath, thisDict):
    thisFH = open(idtablepath, 'r')
    
    for thisLine in thisFH.readlines()[1:]:      # skip first line (header)
        thisLine = thisLine.rstrip()
        thisRead = thisLine.split('\t')
        thisDict[thisRead[0]] = thisRead[1]
    
def ExtractDataset(thisInfile, thisIdDict, thisImputeDir, thisOutPath, thisTempDir, thisBaseName, thisExtractScript):
    thisFH = open(thisInfile, "r")
        
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
    
    # Decode Byte Order Mark
    s = thisFH.read()
    u = s.decode("utf-16")
    s = u.encode("utf-8")
    s = s.rstrip()
    thisIdList = s.split('\r\n')[2:]     # skip the first two header rows
    for thisBiobankID in thisIdList:
        if(thisIdDict.has_key(thisBiobankID)):
            thisSentrixID = thisIdDict[thisBiobankID]
            thisSampleList.append(thisSentrixID + '-' + thisBiobankID)
            print thisSentrixID + '-' + thisBiobankID
        else:
            print "Cannot find key."
    
    # Make temporary outfile for sample names in ascii
    thisSampListFile = 'biobank_x'+str(len(thisSampleList))
    thisSampListFH = open(thisTempDir+'/'+thisSampListFile + '.txt','w')
    for thisSamp in thisSampleList:
        thisSampListFH.write(thisSamp + '\n')
    thisSampListFH.close()
    
    s = drmaa.Session()
    s.initialize()
    print('creating job template for extract_by_chr')
    jt = s.createJobTemplate()
    thisBsubOut = thisTempDir + '/extract_by_chr_%J_%I.out'
    thisBsubErr = thisTempDir + '/extract_by_chr_%J_%I.err'
    jt.nativeSpecification = '-R \'rusage[mem=8192]\' ' + \
                             '-M 12288 ' + \
                             '-q pcpgmwgs ' + \
                             '-J test_array ' + \
                             '-o ' + thisBsubOut + ' ' \
                             '-e ' + thisBsubErr
    
    jt.remoteCommand = thisExtractScript
    jt.args = [thisSampListFile,thisImputeDir,thisOutPath,thisTempDir]
    jobid = s.runBulkJobs(jt, 1, 22, 1)
    
    print('extract_by_chr has been submitted with ID %s' % jobid)
    print('cleaning up')
    s.deleteJobTemplate(jt)
    
    thisMainJobID = int(str(jobid[0]).split('[')[0])
    
    s.exit()
    
    return thisMainJobID, thisSampListFile
    
    
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
        
def MakeReadmeFile(thisReadmeFile, thisOutPath, thisTempDir, thisSampleListName):
    if not os.path.isfile(thisReadmeFile):
        sys.exit("README file does not exist. Please check file path: " + thisReadmeFile)
     
    thisFormDate = time.strftime('%Y-%m-%d')
    thisOutFH = open(thisTempDir + '/README.txt',"w")
     
    with open(thisReadmeFile, "r") as thisReadmeFH:
        for thisLine in thisReadmeFH.readlines():
            thisLine = thisLine.replace("INSERTFILENAME", thisSampleListName)
            thisLine = thisLine.replace("INSERTDATE", thisFormDate)
            thisOutFH.write(thisLine)
     
    thisOutFH.close()
    
    
def BundleToDropbox(thisOutPath, thisBaseName, thisSharePath, thisTempDir, I_compress, I_readonly, I_md5sum, thisUsernameList,thisExtractJobID,thisCopyJobID):
    
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
                             '-w \'done('+str(thisExtractJobID)+') && done('+str(thisCopyJobID)+')\' ' + \
                             '-o ' + thisBsubOut + ' ' \
                             '-e ' + thisBsubErr
    
    jt.remoteCommand = thisTempDir+'/transfer_files.sh'
    
    jobid = s.runJob(jt)
    
    print('transfer_files has been submitted with id %s' % jobid)
    print('cleaning up')
    s.deleteJobTemplate(jt)
    s.exit()
    
    
def setDefaultOutDir(thisBaseName):
    thisUsername = FormatUsername(thisBaseName)

    #Check to make sure you're on panasas
    if(os.path.exists('/data/ppm/biobank/release/')):
        subprocess.call(['mkdir', '-p', '/data/ppm/biobank/release/' + thisUsername])
    else:
        sys.exit('Cannot find path: /data/ppm/biobank/release/')
    return '/data/ppm/biobank/release/' + thisUsername + '/' + thisBaseName + '_impute'
    
def setDefaultTransferDir(thisBaseName):
    thisSharePath = '/pub/dropbox/' + time.strftime('%Y-%m-%d') + '_for_' + FormatUsername(thisBaseName)
        
    #Check to make sure you're on panasas
    if(os.path.exists('/pub/dropbox/')):
        return thisSharePath
    else:
        sys.exit('Cannot find path: /pub/dropbox/')


####################################################################################
################################   Argument Parser   ###############################
####################################################################################
parser = ArgumentParser()

parser.add_argument("-i","--in",
                    type=str,
                    dest="infile",
                    help="csv-version of requested samples")

parser.add_argument("-m","--impute",
                    type=str,
                    dest="imputedir",
                    help="directory containing imputation data for the entire dataset")

parser.add_argument("-o","--out",
                    type=str,
                    dest="outdir",
                    help="output directory for data")

parser.add_argument("-t","--transfer",
                    type=str,
                    dest="transferdir",
                    help="transfer directory for data")

parser.add_argument("-a","--idtable",
                    type=str,
                    dest="idtable",
                    help="Concordance table between biobank ID and Illumina Sentrix ID")

parser.add_argument("-r","--readme",
                    type=str,
                    dest="readmefile",
                    help="corresponding README file for the final ped files")

parser.add_argument("-y","--yaml",
                    type=str,
                    dest="yamlfile",
                    required=True,
                    help="corresponding YAML file for the final ped files")

parser.add_argument("-5","--md5sum",
                    action="store_true",
                    dest="I_md5sum",
                    help="flag for md5sum creation")

parser.add_argument("-u","--usernamelist",
                    type=str,
                    dest="usernamelist",
                    help="designate a list of partners userids to allow access, comma-separated, no spaces!")

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
if(args.idtable != None):
    idTable = args.idtable
else:
    idTable = AssignVariableFromConfig(yml, 'DataFiles', 'ConcordanceTable', "concordance table")

# Assign Readme
if(args.readmefile != None):
    readmeFile = args.readmefile
else:
    readmeFile = AssignVariableFromConfig(yml, 'DataFiles', 'Readme', "readme")

################


# Directories
# infile
if(args.infile != None):
    inFile = args.infile
else:
    inFile = AssignVariableFromConfig(yml, 'Directories', 'InFile', "infile")
    if(inFile == None):
        sys.exit("Need to specify an inFile")

## define baseName
baseName = os.path.basename(inFile).replace('.csv','')

# imputedir
if(args.imputedir != None):
    imputeDir = args.imputedir
else:
    imputeDir = AssignVariableFromConfig(yml, 'Directories', 'ImputeDir', "imputedir")
    if(imputeDir == None):
        sys.exit("Need to specify an imputeDir")

# outfile
if(args.outdir != None):
    outDir = args.outdir
else:
    outDir = AssignVariableFromConfig(yml, 'Directories', 'OutDir', "outdir")
    if(outDir == None):
        if(baseName == ''):
            print "something wrong, no basename"
        else:
            outDir = setDefaultOutDir(baseName)

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
if(args.I_md5sum != None):
    I_md5sum = args.I_md5sum
else:
    I_md5sum = AssignVariableFromConfig(yml, 'Settings', 'Md5sum', "I_md5sum")


# @TODO: currently, the stuff on the bottom don't have commandline parameters
# CompressDirectory
I_compress = AssignVariableFromConfig(yml, 'Settings', 'CompressDir', "I_compress")

# ReadOnlyPermissions
I_readonly = AssignVariableFromConfig(yml, 'Settings', 'ReadOnlyPermissions', "I_readonly")

# Move to dropbox
I_dropbox = AssignVariableFromConfig(yml, 'Settings', 'MoveToDropbox', "I_dropbox")

# username list
if(args.usernamelist != None):
    templist = args.usernamelist.replace(' ', '')
    usernameList = args.usernamelist.split(',')
else:
    username = AssignVariableFromConfig(yml, 'Settings', 'AllowUsernames', "usernameList")
    if(username != None):
        usernameList = username.split(',')
    else:
        usernameList = None


###  Templates for bsub
extractchrScript = AssignVariableFromConfig(yml, 'Scripts', 'ExtractByChr', "extractchrScript")
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


# make sure outDir exists
if not os.path.isdir(outDir):
    # check to see if the directory right above it exists
    if not os.path.isdir('/'.join(outDir.split("/")[:-1])):
        sys.exit("Cannot find path to out dir: " + outDir)
    else:
        # make directory if only that first 
        os.makedirs(outDir)


# set defaults for I_readonly, I_dropbox
if(I_compress == None):
    I_compress = False
if(I_readonly == None):
    I_readonly = False
if(I_dropbox == None):
    I_dropbox = False


# put this in concordance table
IdDict = {}

MakeConcordanceTable(idTable, IdDict)
(thisExtractJobID,thisSampleListName) = ExtractDataset(inFile, IdDict, imputeDir, outDir, tempDir, baseName, extractchrScript)
thisCopyJobID = CopyInfoFiles(outDir, imputeDir, tempDir, thisSampleListName, copyinfoScript)
MakeReadmeFile(readmeFile, outDir, tempDir, thisSampleListName)

if(I_dropbox):
    BundleToDropbox(outDir, baseName, transferDir, tempDir, I_compress, I_readonly, I_md5sum, usernameList,thisExtractJobID,thisCopyJobID)

