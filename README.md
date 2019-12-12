# SCOtags
Some custom scripts related to the article Debray et al. 2020 in prep

# SCOtagsDemultiplexer.py
This script demultiplex SCOtag loci from paired-end reads provided by the EPGV sequencing platform

Options:<br/>
- Required arguments:<br/>
  * `-d STR, --maindir STR`<br/>
                        Path to the main directory<br/>
  * `-p STR, --primers STR`<br/>
                        Path to a 3-column primer file<br/>
  * `-e STR, --error STR`<br/>
                        Levenshtein threshold for primer match<br/>
  * `-f STR, --firmend STR`<br/>
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
