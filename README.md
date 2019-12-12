# RosaPhylogenomics
Some custom scripts related to the article Debray et al. 2020 in prep

# Demultiplexer
This folder contains a script (SCOtagsDemultiplexer.py) that aims to demultiplex SCOtag loci from paired-end reads provided by the EPGV sequencing platform

Options:<br/>
- Required arguments:<br/>
  * `-d STR, --maindir STR`<br/>
                        Path to the main directory<br/>
  * `-p STR, --primers STR`<br/>
                        Path to a 3-column primer file<br/>
  * `-e INT, --error STR`<br/>
                        Levenshtein threshold for primer match<br/>
  * `-f INT, --firmend STR`<br/>
                        Length threshold for firm ends<br/>

- Optional arguments:<br/>
  * `-h, --help`
                        Show the help message and exit<br/>


Usage example:<br/>
`python SCOtagsDemultiplexer.py -d "./" -p SCOtagsPrimerPairs.txt -e 2 -f 4`

Notes:<br/>
- The main directory contains multiplexed, cleaned paired-end fastq primer files. Each accession is coded with three letters and two digits. <br/>
- The primer pairs file is a txt file with 3 columns: first column = SCOtag identifier; second column = forward primer sequence; third column = reverse primer sequence
- The script output a folder for each accession that contain demultiplexed fastq files. In addition, a "Cleaning" folder is created and contains a folder for each accessions that contains only paired-end reads (i.e. reads that whose fwd and rev mates were correctly recovered accordinf to the e and f parameters).



# HybridizationNetworks
This folder contains a script (HybridMapper.py) that aims to create a nexus file that can serve as input for PhyloNet

Options:<br/>
- Required arguments:<br/>
  * `-l STR, --list STR`<br/>
                        Path to a txt file listing the putative parental lineages<br/>
  * `-t STR, --genetrees STR`<br/>
                        Path to a newick style file of gene trees<br/>
  * `-H STR, --hybrid STR`<br/>
                        Name of the putative hybrid<br/>

- Optional arguments:<br/>
  * `-h, --help`
                        Show the help message and exit<br/>

Usage example:<br/>
`python HybridMapper.py -l PutativeParentalLineages.txt -t GeneTrees.txt -H PIM03`

Notes:<br/>
- The last line of HybridMapper.py contain the string that will be used as a command for PhyloNet and can be easily modified to fit one's parameters and options <br/>
- For computational tractability, the list of putative parental lineages contains a few ten accesssion names
- The resulting nexus file can be used in PhyloNet (3.7.1) using the folowing command `java -jar /PATH/TO/PHYLONET/PhyloNet_3.7.1.jar PIM03.surroundings.nexus`
