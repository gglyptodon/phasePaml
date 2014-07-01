import getopt
import re
import sys
import sqlite3
import os
import csv
HEADERSTART = "d.id"

class FastaParser(object):
    def read_fasta(self, fasta, delim = None, asID = 0):
        """read from fasta fasta file 'fasta'
        and split sequence id at 'delim' (if set)\n
        example:\n
        >idpart1|idpart2\n
        ATGTGA\n
        and 'delim="|"' returns ("idpart1", "ATGTGA")
        """
        name = ""
        print(delim,fasta)
        fasta = open(fasta, "r")
        while True:
            line = name or fasta.readline()
            if not line:
                break
            seq = []
            while True:
                name = fasta.readline()
                name = name.rstrip()
                if not name or name.startswith(">"):
                    break
                else:
                    seq.append(name)
            joinedSeq = "".join(seq)
            line = line[1:]
            if delim:
                line = line.split(delim)[asID]
            yield (line.rstrip(), joinedSeq.rstrip())
        fasta.close()


def readCsvExportFasta(f,outdir,dct, prefix = "grp"):
    #todo add id
    if not (os.path.isdir(outdir)):
        os.makedirs(outdir)
    try:
        infile = open(f,"r")
        reader = csv.reader(infile, delimiter="|")
        for r in reader:
            if r[0].startswith(HEADERSTART):
                continue
            outfile = open(outdir+os.sep+prefix+str(r[0])+".fa",'w+')
            for i in r:
                try:
                    outfile.write(">{}\n{}\n".format(i,dct[i]))
                except KeyError as e:
                    if i != r[0]:
                        sys.stderr.write("{} not in Fasta dictionary \n".format(str(i)))
    finally:
        infile.close()
    return(None)
def checkCsv(f):
    fasta2save = set()
    try:
        infile = open(f,"r")
        reader = csv.reader(infile, delimiter="|")
        for r in reader:
            #print(r)
            if r[0].startswith(HEADERSTART):
                continue
            res = [x for x in r if x.startswith("GR") or x.startswith("A") or x.startswith("E")]
            for x in res:
                fasta2save.add(x)

    finally:
        infile.close()
    return(fasta2save.copy())
#read csv, note relevant headers
#read fasta, put relevant seq in dct
#read csv export files to outdir

def usage():
    print ("""
    ##############################
    #
    ##############################
    -i, --infile=XXX.csv
    -h, --help
    -f, --fasta=FASTA
    -o, --out=DIR
    [-p, --prefix=PRE] default: grp
    """)
    sys.exit(2)

def main():
    outdir = "."
    fasta = None
    infile = None
    prefix = "grp"
    global HEADERSTART

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "i:f:o:p:H:h", ["infile=","fasta=", "out=","prefix=","headerstart=","help"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-i", "--infile"):
            infile = a
        elif o in ("-f", "--fasta"):
            fasta = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-o", "--out"):
            outdir = a
        elif o in ("-p", "--prefix"):
            prefix = a
        elif o in ("-H", "--headerstart"):
            HEADERSTART = a
        else:
            assert False, "unhandled option"
    if not fasta:
        usage()
    if not infile:
        usage()

    fasta2save = checkCsv(infile)
    fastaseq = {}
    for f in FastaParser().read_fasta(fasta, None, None):
        if f[0] in fasta2save:
            fastaseq[f[0]] = f[1]
            print(f[0])
    readCsvExportFasta(f=infile, outdir=outdir, dct=fastaseq, prefix=prefix)

####################################
if __name__ == "__main__":
    main()
####################################
