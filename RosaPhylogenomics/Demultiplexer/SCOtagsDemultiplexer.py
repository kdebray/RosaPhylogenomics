#! /bin/python



# SCOtagsDemultiplexer.py
# Written by K Debray
# Python v3.6.3

""" Aims to demultiplex genes from cleaned reads provided by the EPGV sequencing platform """

# Packages import
#--------------------------------------------------------------------------------------------------
import os, sys, glob, argparse
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Levenshtein import distance

parser = argparse.ArgumentParser()
# Required arguments
parser.add_argument("-d", "--maindir", dest="d", required=True, type=str, help="Path to maindir")
parser.add_argument("-p", "--primers", dest="p", required=True, type=str, help="Path to 3-column primer file")
# Optional arguments
parser.add_argument("-e", "--error", dest="e", required=False, type=int, default=0, \
                    help="Levenshtein threshold for primer match")
parser.add_argument("-f", "--firmend", dest="f", required=False, type=int, default=0, help="Length threshold for firm end")

args = parser.parse_args()

MainDir = args.d
Primers = args.p
SpecPrimErrThr = args.e
EndFirmThr = args.f

# Function Parts
#--------------------------------------------------------------------------------------------------

def RemoveUnpaired(file1, file2, outdir):
    # Retrieve reads id in file1 and 2
    FwdReadListSorted = [f.id.split('/')[0] for f in sorted(SeqIO.parse(file1, "fastq"), \
        key=lambda x : x.id)]
    RevReadListSorted = [f.id.split('/')[0] for f in sorted(SeqIO.parse(file2, "fastq"), \
        key=lambda x : x.id)]

    # Compare
    Union = set(FwdReadListSorted) & set(RevReadListSorted)

    # Save F file
    OutFname = outdir + file1.split('/')[-1].rsplit('.')[0] + '.paired.fastq'
    with open(OutFname, 'w') as handleF:
        for title, seq, qual in FastqGeneralIterator(open(file1)):
            if title.split('/')[0] in Union:
                handleF.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

    # Save R file
    OutRname = outdir + file2.split('/')[-1].rsplit('.')[0] + '.paired.fastq'
    with open(OutRname, 'w') as handleR:
        for title, seq, qual in FastqGeneralIterator(open(file2)):
            if title.split('/')[0] in Union:
                handleR.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))


# Main
#--------------------------------------------------------------------------------------------------

# Check that requirements are fulfilled
#--------------------------------------------------------------------------------------------------

print("\nChecking requirements...")

trouble = 0
# Does primer file have 3 tab-separated columns?
with open(Primers, "r") as PrimerFile:
    # Check that primer file contain 3 tab-separated columns
    for line in PrimerFile:
        if len(line.split("\t")) != 3:
            trouble += 1
            print("The primer file has not 3 columns. Please check and restart.")
            sys.exit()

# Does the main directory provided by user ends with a "/"
if MainDir[-1] != "/":
    trouble += 1
    print("Error: the main directory you provided does not end with a '/'. Please restart.")
    sys.exit()

if trouble == 0:
    print("... OK: Requirements are fulfilled\n")


# Demultiplexing
#-------------------------------------------------------------------------------------------------

print("\nDemultiplexing...")

# Create individual directories
l = []
for file in glob.glob(MainDir+"*.fastq"):
    l.append(file.split("/")[-1].split("_")[0])
IndList = sorted(set(l))

# Demultiplex
for ind in IndList:
    print("...Demultiplexing for individual %s" % ind)
    if not os.path.exists(ind):
        os.makedirs(ind)

    for PrimOrientation in ["F", "R"]:
        # Get the reads file
        ReadsFile = MainDir + ind + "_" + PrimOrientation + ".fastq"
    
        # Loop through the primers files
        with open(Primers, 'r') as PrimersFile:
            for line in PrimersFile:
                PairName = line.strip().split()[0]
                if PrimOrientation == "F":
                    PrimerSeq = line.strip().split()[1]
                elif PrimOrientation == "R":
                    PrimerSeq = line.strip().split()[2]
                PrimerLen = len(PrimerSeq)
            
                # Loop through fastq file
                with open(ReadsFile) as reads:
                    for title, seq, qual in FastqGeneralIterator(reads):
                        # extract the corresponding primer length
                        ReadExtract = seq[:PrimerLen]
                        # calculate the Levenshtein distance
                        Ldist = distance(ReadExtract, PrimerSeq)
                        ## Assign read to a gene according to the 2 thresholds
                        #-------------------------------------------------------------
                        # Perform the first test (Levenshtein distance)
                        if Ldist <= SpecPrimErrThr:
                            # Perform the second test (3' firm end)
                            ReadEnd = ReadExtract[-EndFirmThr:]
                            PrimEnd = PrimerSeq[-EndFirmThr:]
                            if ReadEnd == PrimEnd:
                                # Create an outfile if doesnt exist
                                if not os.path.exists('./'+ind+'/'+ ind + '_' + PairName + '_' + PrimOrientation + \
                                    '.fastq'):
                                    fpi = open('./'+ind+'/'+ ind + '_' + PairName + '_' + PrimOrientation + \
                                        '.fastq', 'w')
                                    fpi.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                                    fpi.close()
                                # Otherwise append to existing file
                                else:
                                    with open('./'+ind+'/'+ ind + '_' + PairName + '_' + PrimOrientation + '.fastq', 'a') as out:
                                        out.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

# Cleaning
#----------------------------------------------------------------------------------------

print("\nCleaning...")

# Create a clean directory
if not os.path.exists('Cleaning'):
    os.makedirs('Cleaning')

for ind in IndList:
    print("...Cleaning individual %s" % ind)
    # Create an individual directory in Cleaning/
    if not os.path.exists('./Cleaning/'+ind):
        os.makedirs('./Cleaning/'+ind)

    # Get the genes files
    genes = sorted(list(set([f[6:13] for f in os.listdir('./'+ind+'/')])))
    #print(genes)
    for gene in genes:
        # Check if 2 files are present
        Pairs = sorted(glob.glob('./'+ind+'/*'+gene+'*'))
        # If OK, remove unpaired reads
        if len(Pairs) == 2 and 'F.fastq' in Pairs[0] and 'R.fastq' in Pairs[1]:
            RemoveUnpaired(Pairs[0], Pairs[1], './Cleaning/'+ind+'/')

