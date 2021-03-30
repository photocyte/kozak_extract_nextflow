#! /usr/bin/env python
import argparse
import sys
import pkg_resources
import Bio
import Bio.SeqIO
import re
import gzip
import os

try:
    gffutils_version = pkg_resources.get_distribution("gffutils").version
    import gffutils
except:
    print("This script requires gffutils version >= 0.9. See http://daler.github.io/gffutils")
    exit(1)



#if float(gffutils_version) >= 0.899:
#        import gffutils
#else:
#        print("This script requires gffutils version >= 0.9. See http://daler.github.io/gffutils")
#        exit(1)


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-g")
parser.add_argument("--i",default=20)
parser.add_argument("--in_memory",default=False,action='store_true')
parser.add_argument("--transdecoder_mode",default=False,action='store_true')
parser.add_argument("--sqlite_name",default="gffutils.sqlite")
parser.add_argument("--gff3_compliant",default=False,action='store_true')
args = parser.parse_args()

db_path=args.g+".gffutils.db"
sys.stderr.write("Reading GFF3 file: "+args.g+"\n")
sys.stderr.flush()

from_string_switch = False
if args.g[-3:].lower() == ".gz":
    sys.stderr.write("Guessing it is a GZIP file based on the suffix. Decompressing into memory...\n")
    from_string_switch = True
    handle = gzip.open(args.g,"r")
    args.g = handle.read()
    sys.stderr.write("**NOT IMPLEMENTED** exiting...")
    exit(1)

if args.in_memory == True:
    sys.stderr.write("Coverting to in memory gffutils sqlite database\n")
    db = gffutils.create_db(args.g,":memory:", force=True,merge_strategy="create_unique",from_string=from_string_switch)
else:
    if not os.path.isfile(args.sqlite_name):
        sys.stderr.write("Coverting to "+args.sqlite_name+"gffutils sqlite database\n")
        db = gffutils.create_db(args.g,args.sqlite_name, merge_strategy="create_unique",from_string=from_string_switch)
    else:
        sys.stderr.write("sqlite database seems to be already written at "+args.sqlite_name+" \n")
        db = gffutils.FeatureDB(args.sqlite_name)
sys.stderr.write("Done converting. Now printing modified GFF3 to stdout...\n")
sys.stderr.flush()
i=0
if args.gff3_compliant:
   sys.stdout.write("##gff-version 3\n")
for g in db.features_of_type("gene"):
    if args.transdecoder_mode and "complete" not in g["Name"][0]:
        continue
    if args.gff3_compliant:
        print(g)
    for m in db.children(g.id,featuretype='mRNA'):
        if args.gff3_compliant:
            print(m)
    for c in db.children(g.id,featuretype='CDS'):
        if c.strand == "+":
            new_start = c.start - args.i
            new_stop = c.start + 14
        if c.strand == "-":
            new_start = c.stop - 14
            new_stop = c.stop + args.i
        if new_start < 1:
            continue
        if new_stop > g.stop:
            continue
        
        new_attrs = []
        for a in c.attributes:
            ##sys.stderr.write(a+";;"+str(c.attributes[a])+'\n')
            if len(c.attributes[a]) == 0:
                continue
            new_attrs.append(a+"="+c.attributes[a][0])
        new_attr_string = ";".join(new_attrs)
        gene_string = '\t'.join([c.chrom,c.source,"kozakseq",str(new_start),str(new_stop),c.score,c.strand,c.frame,new_attr_string])
        sys.stdout.write(gene_string+"\n")

