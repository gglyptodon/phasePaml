import getopt
import re
import sys
import sqlite3
import os
import csv
DBNAME = None;
ID0 = "phytozome9_0__5a_59"
ID1 = "maizeseq__4a_53"
HEADERSTART = "5a.59_"

def checkExists():
    global DBNAME
    if os.path.isfile(DBNAME):
        return(True)
    else:
        return(False)

def createTable(table):
    global DBNAME
    db = sqlite3.connect(DBNAME)
    cmd = "\
    CREATE TABLE IF NOT EXISTS {}\
    (\
    id INTEGER PRIMARY KEY,\
    phytozome9_0__5a_59 TEXT UNIQUE,\
    maizeseq__4a_53 TEXT,\
    confidence INTEGER\
    );\
    ".format(table)
    #print(cmd)
    c = db.cursor()
    c.execute(cmd)
    db.commit()
    db.close()
    
def insertDB(db, cursor,res, table):
   
    ID0v = res[0]
    ID1v =res[1]
    confv = res[2]
    print(confv)
    #db = sqlite3.connect(DBNAME)
    #c = db.cursor()
    cmd = "INSERT OR IGNORE INTO {} ({}, {}, {}) VALUES ({}, {},{})".format(table, ID0, ID1,"confidence",q(ID0v),q(ID1v),confv)
    #print(cmd)
    cursor.execute(cmd)
    #db.commit()
    #return(None)

def readCsvIntoDB(f,table):
    db = sqlite3.connect(DBNAME)
    c = db.cursor()
    try:
        infile = open(f,"r")
        reader = csv.reader(infile)
        for r in reader:
            if r[0].startswith(HEADERSTART):
                continue
            res = ["NULL" if x=="-" else x for x in r ]
           
            insertDB(cursor=c, db=db, res=res, table=table)
    finally:
        infile.close()
        db.commit()
        db.close()
        return(None)

def q(s):
    return('"'+s+'"')

def usage():
    print ("""
    ##############################
    #
    ##############################
    -i, --infile=tableofdoom.csv
    -h, --help
    -t, --table=TABLENAME
    [-n, --new=DBNAME or   # creates a new database DBNAME
    -o, --old=DBNAME ]     # inserts into existing database DBNAME
    
    """)
    sys.exit(2)

def main():
    global DBNAME;
    createNew=False;
    fillOld = False;
    table=None;
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
            table = a
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
    if not table:
        table="ZmIDMap"
     ##################################
    if fillOld:
        if not checkExists():
            sys.stderr.write("There is no such file {} \n".format(DBNAME))
        else:
            print("#ok\n")
            createTable(table)
            readCsvIntoDB(infile, table)
    if createNew:
        if checkExists():
            sys.stderr.write("{} already exists.\n".format(DBNAME))
            exit(1)
        else:
            print("#ok\n")
            createTable(table)
    #readCsv(infile)
    pass
####################################
if __name__ == "__main__":
    main()
####################################
