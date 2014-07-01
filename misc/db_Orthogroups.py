import getopt
import re
import sys
import sqlite3
import os
import csv
DBNAME = None
ORTHO = "Ortho5G"


HEADERSTART = "###"

def checkExists():
    global DBNAME
    if os.path.isfile(DBNAME):
        return(True)
    else:
        return(False)


# throw out all tandem duplicates for this table
# try to find GRMZ id as representative
def createOrtho(table=ORTHO):
    global DBNAME
    db = sqlite3.connect(DBNAME)    
    cmd = "\
    CREATE TABLE IF NOT EXISTS {}\
    (\
    id INTEGER PRIMARY KEY,\
    group TEXT ,\
     ,\
    maizeB1 TEXT ,\
    maizeB2 TEXT \
    );\
    ".format(table)
    c = db.cursor()
    c.execute(cmd)
    db.commit()
    db.close()
    
    
def insertDBZmDup(db, cursor,res, table):
    A1 = res[0]
    A2 = res[1]
    B1 = res[2]
    B2 = res [3]
    
    #db = sqlite3.connect(DBNAME)
    #c = db.cursor()
    cmd = "INSERT OR IGNORE INTO {} ({}, {}, {}, {}) VALUES ({}, {},{},{})".format(table, "maizeA1","maizeA2","maizeB1", "maizeB2",q(A1),q(A2),q(B1),q(B2))
    print(cmd)
    cursor.execute(cmd)
    #db.commit()
    #return(None)

def fetchGRM(l):
    for i in l:
        if(i.startswith("GRZM")):
            return(i)
    for i in l:
        if (i.startswith("AC")):
            return(i)
    return(l[0])
def readCsvIntoDB(f,table):
    db = sqlite3.connect(DBNAME)
    c = db.cursor()
    try:
        infile = open(f,"r")
        reader = csv.reader(infile)
        for r in reader:
            if r[0].startswith(HEADERSTART):
                continue
            #print(r)
            #take GRMZM preferably
            res = [fetchGRM(x.split("||")) if len(x.split("||"))>1 else x for x in r ]
            res = ["NULL" if x=="" else x for x in res ]
            #print(res)
            insertDBZmDup(cursor=c, db=db, res=res, table=table)
    finally:
        infile.close()
        db.commit()
        db.close()
    return(None)

def q(s):
    if s=="NULL":
        return(s)
    return('"'+s+'"')

def usage():
    print ("""
    ##############################
    #
    ##############################
    -i, --infile=evs900_zm.csv
    -h, --help
    [-t, --table=TABLENAME] default:"ZmDup"
    [-n, --new=DBNAME or   # creates a new database DBNAME
    -o, --old=DBNAME ]     # inserts into existing database DBNAME
    
    """)
    sys.exit(2)
#select a2.phytozome9_0__5a_59,b2.phytozome9_0__5a_59 from 
#ZmIDMap a2,ZmIDMap b2, ancDupA2_B2 d where a2.maizeseq__4a_53=d.maizeA2 
#AND b2.maizeseq__4a_53=d.maizeB2;
def main():
    global DBNAME
    global DUP
    createNew=False
    fillOld = False
    table=None
    #if fillOld: try CREATE TABLE IF NOT EXISTS 
    # if create same
    
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "i:t:n:o:h", ["infile=","table=", "new=", "old=","help"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-i", "--infile"):
            infile = a
        elif o in ("-t", "--table"):
            DUP = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-n", "--new"):
            DBNAME = a
            createNew=True;
            fillOld = False
        elif o in ("-o", "--old"):
            fillOld=True
            DBNAME=a
            creatNew=False
        else:
            assert False, "unhandled option"
    print(DBNAME, fillOld, createNew)
    if not DBNAME:
        usage()
    if not createNew and not fillOld:
        usage()
    if not infile:
        usage()

        
     ##################################
    if fillOld:
        if not checkExists():
            sys.stderr.write("There is no such file {} \n".format(DBNAME))
        else:
            print("#ok\n")
            createTableZmDup(DUP)
            readCsvIntoDB(infile, DUP)
            createViews()
    if createNew:
        if checkExists():
            sys.stderr.write("{} already exists.\n".format(DBNAME))
            exit(1)
        else:
            print("#ok\n")
            createTableZmDup(DUP)
            readCsvIntoDB(infile, DUP)
            createViews()
    #readCsv(infile)
    pass
####################################
if __name__ == "__main__":
    main()
####################################
