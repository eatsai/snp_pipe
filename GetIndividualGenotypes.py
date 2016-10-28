#!/usr/bin/env python

# Purpose:       Get the individual genotypes 
# 
# Usage:         python GetIndividualGenotypes.py ...
# 
# Author:        et85, etsai@bwh.harvard.edu


from argparse import ArgumentParser
import os, re, shutil, subprocess, sys, time, yaml

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
    
def MakePedFile(thisInfile, thisIdDict, thisPedPath, thisOutPath, thisBaseName, thisFullPed, I_allsamp):
    thisFH = open(thisInfile, "r")
        
    # Make output dir if it doesn't exist.
    if os.path.exists(thisOutPath):
        if(os.listdir(thisOutPath) != []):
            sys.exit("This out dir is not empty. Please delete the contents in the out dir before trying again:\n" + thisOutPath)
    else:
        os.makedirs(thisOutPath)
    
    
    if(I_allsamp):    #  just copy the full PED file over
        shutil.copy(thisFullPed,thisOutPath + '/' + thisBaseName + '.ped')
    
    else:   # go through the CSV file one by one
        thisOutFH = open(thisOutPath + '/' + thisBaseName + '.ped','w')
        # Decode Byte Order Mark
        s = thisFH.read()
        u = s.decode("utf-16")
        s = u.encode("utf-8")
        s = s.rstrip()
        thisIdList = s.split('\r\n')[2:]     # skip the first two header rows
        for thisBiobankID in thisIdList:
            if(thisIdDict.has_key(thisBiobankID)):
                thisSentrixID = thisIdDict[thisBiobankID]
                print thisBiobankID, thisSentrixID
                
                with open(thisPedPath + '/' + thisSentrixID + '.ped', 'r') as thisPedFH:
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
    
def MakeReadmeFile(thisReadmeFile, thisOutPath, thisBaseName):
    if not os.path.isfile(thisReadmeFile):
        sys.exit("README file does not exist. Please check file path: " + thisReadmeFile)
    
    thisFormDate = time.strftime('%Y-%m-%d')
    thisOutFH = open(thisOutPath + '/README.txt',"w")
    
    with open(thisReadmeFile, "r") as thisReadmeFH:
        for thisLine in thisReadmeFH.readlines():
            thisLine = thisLine.replace("INSERTFILENAME", thisBaseName)
            thisLine = thisLine.replace("INSERTDATE", thisFormDate)
            thisOutFH.write(thisLine)
    
    thisOutFH.close()
    
def BundleToDropbox(thisOutPath, thisBaseName, thisSharePath, I_readonly, I_compress, I_md5sum, thisUsernameList):
    thisRead = thisOutPath.split("/")
    thisOutDirName = thisRead[-1]
    
    if(I_readonly):
        os.chmod(thisOutPath + '/' + thisBaseName + '.ped', 0440)
        os.chmod(thisOutPath + '/' + thisBaseName + '.map', 0440)
        os.chmod(thisOutPath + '/' + thisBaseName + '.snpinfo.txt', 0440)
        os.chmod(thisOutPath + '/README.txt', 0400)
        subprocess.call(['chmod', '500', thisOutPath])
    
    if(I_compress):
        subprocess.call(['tar', 'czvf', thisOutPath + '.tgz', '-C', thisOutPath + '/../', thisOutDirName])
        if(I_md5sum):
            thisMd5FH = open(thisOutPath + '.tgz.md5', 'w')
            subprocess.call(['md5sum', thisOutPath + '.tgz'], stdout=thisMd5FH)
            thisMd5FH.close()
    
    if not os.path.exists(thisSharePath):
        subprocess.call(['mkdir', thisSharePath])
    else:
        print 'This dir already exists: ' + thisSharePath
    
    if(I_compress):
        # Check to make sure that these files do not already exist
        if(os.path.exists(thisSharePath + '/' + thisBaseName + '.tgz')):
            print thisSharePath + '/' + thisBaseName + '.tgz' + " already exists. No replacement attempted"
        else:
            subprocess.call(['cp', '-pr', thisOutPath + '.tgz', thisSharePath])
            if(I_md5sum):
                subprocess.call(['cp', '-pr', thisOutPath + '.tgz.md5', thisSharePath])
        
        if(I_readonly):
            # remove the directory
            subprocess.call(['chmod', '700', thisOutPath])
        
        # remove the directory file
        subprocess.call(['rm', '-rf', thisOutPath])
    else:
        if(os.path.exists(thisSharePath + '/' + thisBaseName)):
            print thisSharePath + '/' + thisBaseName + " already exists. No replacement attempted"
        else:
            subprocess.call(['cp', '-pr', thisOutPath, thisSharePath])
    
    if(thisUsernameList != None):
        for thisName in thisUsernameList:
            subprocess.call(['setfacl', '-R', '-m', 'u:' + thisName + ':rwX', thisSharePath])

def setDefaultOutDir(thisBaseName):
    thisUsername = FormatUsername(thisBaseName)

    #Check to make sure you're on panasas
    if(os.path.exists('/data/ppm/biobank/release/')):
        subprocess.call(['mkdir', '-p', '/data/ppm/biobank/release/' + thisUsername])
    else:
        sys.exit('Cannot find path: /data/ppm/biobank/release/')
    return '/data/ppm/biobank/release/' + thisUsername + '/' + thisBaseName
    
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

parser.add_argument("-p","--ped",
                    type=str,
                    dest="peddir",
                    help="ped directory for data location")

parser.add_argument("-m","--map",
                    type=str,
                    dest="mapfile",
                    help="corresponding map file for the final ped files")

parser.add_argument("-s","--snpinfo",
                    type=str,
                    dest="snpinfofile",
                    help="corresponding snpinfo file for the final ped files")

parser.add_argument("-r","--readme",
                    type=str,
                    dest="readmefile",
                    help="corresponding README file for the final ped files")

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

parser.add_argument("-f","--fullped",
                    type=str,
                    dest="fullped",
                    help="corresponding full PED file for data release")

parser.add_argument("-z","--allsamp",
                    action="store_true",
                    dest="I_allsamp",
                    help="set string for delivering the designated fullped file")

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

# Assign mapFile
if(args.mapfile != None):
    mapFile = args.mapfile
else:
    mapFile = AssignVariableFromConfig(yml, 'DataFiles', 'Map', "map")

# Assign snpInfoFile
if(args.snpinfofile != None):
    snpInfoFile = args.snpinfofile
else:
    snpInfoFile = AssignVariableFromConfig(yml, 'DataFiles', 'SnpInfo', "snp info")

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

# Assign FullPed file
if(args.fullped != None):
    fullPedFile = args.fullped
else:
    fullPedFile = AssignVariableFromConfig(yml, 'DataFiles', 'FullPed', "fullped")

# Deliver FullPed file
if(args.I_allsamp != None):
    I_allsamp = args.I_allsamp
else:
    I_allsamp = AssignVariableFromConfig(yml, 'Settings', 'FullDataset', "I_allsamp")

if(I_allsamp):
    if(fullPedFile == None):
        sys.exit('Cannot find path to full ped file')
else:
    fullPedFile = ''
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

# peddir
if(args.peddir != None):
    pedDir = args.peddir
else:
    pedDir = AssignVariableFromConfig(yml, 'Directories', 'PedDir', "pedir")

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

# @TODO: basedir is currently unused. Not sure if we need this?


### Other Settings

# check inFile exists
if(not I_allsamp):
    # make sure infile exists
    if not os.path.isfile(inFile):      # don't need to rely on the inFile for full dataset delivery
        sys.exit("Cannot find input file: " + inFile)


# there is still a dummy input file designated, so ignore this for now
# if(I_allsamp and not os.path.isfile(fullPedFile)):
#     sys.exit('Cannot find fullped file to copy to output directory')

if(args.I_compress != None):
    I_compress = args.I_compress
else:
    I_compress = AssignVariableFromConfig(yml, 'Settings', 'CompressDir', "I_compress")

if(args.I_md5sum != None):
    I_md5sum = args.I_md5sum
else:
    I_md5sum = AssignVariableFromConfig(yml, 'Settings', 'Md5sum', "I_md5sum")


# @TODO: currently, the stuff on the bottom don't have commandline parameters
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


#### Make sure the directories exist    


# make sure outDir exists
if not os.path.isdir(outDir):
    # check to see if the directory right above it exists
    if not os.path.isdir('/'.join(outDir.split("/")[:-1])):
        sys.exit("Cannot find path to out dir: " + outDir)
    else:
        # make directory if only that first 
        os.makedirs(outDir)

# make sure ped path exists
if not os.path.isdir(pedDir):
    sys.exit("Cannot find path to ped dir: " + pedDir)

# set defaults for I_readonly, I_dropbox
if(I_readonly == None):
    I_readonly = False
if(I_dropbox == None):
    I_dropbox = False


# put this in concordance table
IdDict = {}

MakeConcordanceTable(idTable, IdDict)
MakePedFile(inFile, IdDict, pedDir, outDir, baseName, fullPedFile, I_allsamp)
CopySnpFiles(mapFile, snpInfoFile, outDir, baseName)
MakeReadmeFile(readmeFile, outDir, baseName)
  
 
if(I_dropbox):
    BundleToDropbox(outDir, baseName, transferDir, I_readonly, I_compress, I_md5sum, usernameList)

    




