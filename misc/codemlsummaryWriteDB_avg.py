import getopt
import re
import sys
import sqlite3
import os
def checkDB(dbname):
    pass
    db = sqlite3.connect("scythe_codeml_summary.db")
    c = db.cursor()
    #Key:GrpDatasetGeneAGeneB
    cmd = 'CREATE TABLE IF NOT EXISTS codeml_w_summary_prank (key TEXT PRIMARY KEY NOT NULL, dataset TEXT, grp TEXT, w_longest_min FLOAT, w_longest_max FLOAT, w_longest_avg FLOAT, w_def_min FLOAT, w_def_max FLOAT, w_def_avg FLOAT, w_max_min FLOATw_max_max FLOAT, w_max_avg FLOAT, w_sum_min FLOAT, w_sum_max FLOAT, w_sum_avg FLOAT);'
    #print(cmd)
    c.execute(cmd)
    db.commit()
    #cmd = 'CREATE TABLE IF NOT EXISTS codeml_w_clustal (key TEXT PRIMARY KEY NOT NULL, dataset TEXT, grp TEXT, GeneA TEXT, GeneB TEXT, w_longest FLOAT, w_def FLOAT, w_max FLOAT, w_sum FLOAT);'

    #c.execute(cmd)
    #print(cmd)
    #db.commit()
    
    db.close()
def readSummaryTable_codeml_w(infile):
    codeml_w = None
    grp =None
    geneA = None
    geneB = None
    infileh = open(infile, 'r')
    for ln in infileh:
        ln.rstrip()
        if ln.startswith("file"):
            pass
        else:
            l = ln.split("\t")
            #print(l)
            grp = l[0].split(".")[0]
            geneA = l[1]
            geneB = l[2]
            codeml_w = l[7]
            yield(grp, geneA, geneB, codeml_w)

def insertDBCodemlW(dataset, grp, geneA, geneB, aligner, method, value):
    #method eg w_longest, w_def
    db = sqlite3.connect("scythe_codeml.db")
    c = db.cursor()
    geneAt = geneA.replace("|","_")
    geneBt= geneB.replace("|","_")
    #print(grp)
    k = "".join([dataset,grp,geneAt,geneBt])
    dataset = q(dataset)
    grp = q(grp)
    k.replace("|","_")
    k = q(k)
    method =q(method)
    group = grp 
    #print(k)
    if aligner == "prank":
        #try:
        cmd = 'INSERT OR IGNORE INTO codeml_w_prank (key, dataset, grp, {}) VALUES ({}, {}, {}, {});'.format(method, k,dataset,group,value)
        #print(cmd)
        c.execute(cmd)
        db.commit()
        #except sqlite3.IntegrityError:
        cmd = 'UPDATE codeml_w_prank SET dataset={}, grp={}, {}={} WHERE key={};'.format(dataset,group,method,value,k)
        #print(cmd)
        #db.commit()
        c.execute(cmd)
        db.commit()
    elif (aligner=="clustal"):
        cmd = 'INSERT OR IGNORE INTO codeml_w_clustal  (key, dataset, grp, {}) VALUES ({},{}, {}, {},{});'.format(method,k, dataset,group,value)
        #print(cmd)
        c.execute(cmd)
        db.commit()
        cmd = 'UPDATE codeml_w_clustal SET dataset={}, grp={}, {}={} WHERE key={};'.format(dataset,group,method,value,k)
        #print(cmd)
        #db.commit()
        c.execute(cmd)
        db.commit()
    else:
        print("No db for ", aligner)
    db.close()
def q(s):
    return('"'+s+'"')
def usage():
    print ("""
    ##############################
    #
    ##############################
    -i, --infile=ALIGNMENT
    -h, --help
    -D, --dataset=NAME
    -M, --method=[w_longest, w_sum, w_max, w_def]
    """)
    sys.exit(2)

#def readFasta(fasta,stfu=True):
#    alignment = ""
#    alignment_length = 0
#    number_species = 0
#    number_gaps = None
#    infile = open(fasta,'r')
#    tmp = infile.read()
#    for ln in tmp:
#        if ln.startswith(">"):
#            number_species+=1
#        else:
#            ln = ln.strip()
#            #print(len(ln))
#            alignment += ln
#    alignment_length= len(alignment)
#    number_gaps = alignment.count("-")
#    perc_gaps = number_gaps/alignment_length*100
#    if not stfu:
#        print(number_gaps)
#    return(perc_gaps)

#def guessFormat(infile, stfu=True):
#    if not stfu:
#        print("guessing format: believe it's... ", end="")
#    infile = open(infile,'r')
#    ln = infile.readline()
#    infile.close()
#    if ln.startswith("#NEXUS"):
#        if not stfu:
#            print("nexus.\n")
#        return ("nexus")
#    elif ln.startswith(">"):
#        if not stfu:
#            print("fasta.\n")
#        return("fasta")
#    else:
#        sys.stderr.write("unknown format\n")
#        return(None)


def main():

    ###################################
    infile = None
    stfu = True
    verbose = False
    dataset = None
    dbname = "scythe_codeml.db"
    grp = None
    info = None
    aligner = None
    method = None
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "i:D:M:vVh", ["infile=","dataset=","method=","verbose","Verbose","help"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-i", "--infile"):
            infile = a
        elif o in ("-M", "--Method"):
            method = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-v", "--verbose"):
            stfu=False
        elif o in ("-V", "--Verbose"):
            verbose=True
        elif o in ("-D","--dataset"):
            dataset = a
        
        elif o in ("-I","--info"):
            info = a
        else:
            assert False, "unhandled option"

    if not infile or not dataset or not method:
        usage()
    if method not in ["w_longest", "w_def", "w_max", "w_sum"]:
        usage()
    else:
        aligner="prank"
        if dbname:
            grp = infile.split(os.sep)[-1].split(".")[0]
            checkDB(dbname)
            for i in readSummaryTable_codeml_w(infile):
                grp, geneA,geneB, codeml_w = i[0],i[1],i[2],i[3]
                print(aligner, method, dataset, grp, codeml_w, geneA, geneB)

                insertDBCodemlW(aligner=aligner, method=method,dataset = dataset, grp = grp,value=codeml_w, geneA=geneA, geneB=geneB)
####################################
if __name__ == "__main__":
    main()
####################################
