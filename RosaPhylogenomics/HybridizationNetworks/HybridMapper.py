#! /bin/python

# HybridMapper.py
# Written by K Debray
# Python v3.6.3

""" Aims to create nexus file input for PhyloNet 3.7.1"""

# Packages import
#--------------------------------------------------------------------------------------------------
import re
import argparse
from ete3 import Tree
from collections import defaultdict


# Arguments section
#--------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
# List of BB 2x NH species
parser.add_argument("-l", "--taxalist", dest="list", required=True, type=str, \
    help="Path to a txt file listing the putative parental lineages")

# Set of genetrees in a file
parser.add_argument("-t", "--genetrees", dest="trees", required=True, type=str, \
    help="Path to a newick file of gene trees")

# Name of the hybrid accession
parser.add_argument("-H", "--hybrid", dest="hybrid", required=True, type=str, \
    help="Name of the putative hybrid")

args = parser.parse_args()

BB = args.list #BackBone
GTfile = args.trees #Gene Trees file
Hyb = args.hybrid #Hybrid


# Main
#--------------------------------------------------------------------------------------------------

# Step 1: Pruning
#--------------------------------------------------------------------------------------------------
Leaves2Keep = []
with open(BB, 'r') as bb:
    for l in bb:
        Leaves2Keep.append(l.strip())

i = 0
with open(GTfile, 'r') as gtfile:
    for l in gtfile:
        t = Tree(l.strip(), format=0)

        MyLeaves = []
        MyTest = False

        for leaf in t:
            if leaf.name.startswith(Hyb) and leaf.name not in MyLeaves:
                MyTest = True
                MyLeaves.append(leaf.name)
            elif leaf.name[:5] in Leaves2Keep and leaf.name not in MyLeaves:
                MyLeaves.append(leaf.name)
        if MyTest == True:
            t.prune(MyLeaves, preserve_branch_length=True)

            with open(GTfile.rsplit('.',1)[0]+'.pruned.txt', 'a+') as out:
                out.write(t.write(format=0))
                out.write('\n')


# Step 2: Preparing PhyloNet nexus file
#--------------------------------------------------------------------------------------------------
TreeDico = {}
AlleleDico = defaultdict(list)

with open(GTfile.rsplit('.',1)[0]+'.pruned.txt', 'r') as inprune:
    for i, l in enumerate(inprune):
        t = Tree(l.strip(), format=0)
        for leaf in t:
            if not re.match(r"[A-Z]{3}[0-9]{2}_[0-9]{1}", leaf.name):
                leaf.name = leaf.name + '_1'

            if leaf.name not in AlleleDico[leaf.name[:5]]:
                 AlleleDico[leaf.name[:5]].append(leaf.name)

        TreeDico['Tree g'+str(i)] = t.write(format=0)

    
TaxaMap = '-a <'
for key in AlleleDico:
    s = ''
    s += key+':'
    for i in AlleleDico[key]:
        s += i+','
    s = s[:-1] + '; '
    
    TaxaMap += s
TaxaMap = TaxaMap[:-2] + '>'

# Build nexus file
with open(Hyb + '.surroundings.nexus', 'w') as outfile:
    outfile.write("#NEXUS\n\nBEGIN TREES;\n\n")
    for t in TreeDico:
        outfile.write("%s = %s\n" % (t, TreeDico[t]))
    outfile.write("\n\nEND;\n\nBEGIN PHYLONET;\n\nInferNetwork_MPL (all) 1 -h {%s} -x 5 -o -pl 4 -di %s;\n\nEND;" % (Hyb, TaxaMap))


