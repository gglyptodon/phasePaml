import getopt
import re
import sys
import sqlite3
import os
import csv
DBNAME = None
TANDEM = "ZmTandemDup"


HEADERSTART = "Maize1"

def checkExists():
    global DBNAME
    if os.path.isfile(DBNAME):
        return(True)
    else:
        return(False)
def createTableZmTandemDup(tablename):
    db = sqlite3.connect(DBNAME)    
    cmd = "\
    CREATE TABLE IF NOT EXISTS {}\
    (\
    id INTEGER PRIMARY KEY,\
    tandemID INTEGER ,\
    tandem1 TEXT ,\
    tandem2 TEXT ,\
    info TEXT ,\
    tandemcount INTEGER \
    );\
    ".format(tablename)
    c = db.cursor()
    c.execute(cmd)
    db.commit()
    db.close()
def insertPairwise(db, cursor, pairslist, info, table,counter):
    for i in range(0, len(pairslist)-1):
        for j in range(i+1, len(pairslist)):
            print(i,j)
            cmd = "INSERT OR IGNORE INTO {} ({}, {}, {},{}) VALUES ({}, {}, {}, {})".format(table, "tandem1", "tandem2","info", "tandemcount", q(pairslist[i]),q(pairslist[j]),q(info), len(pairslist))
            cmd2 = "INSERT OR IGNORE INTO {} ({}, {}, {},{}) VALUES ({}, {}, {}, {})".format(table, "tandem1", "tandem2","info", "tandemcount", q(pairslist[j]),q(pairslist[i]),q(info),len(pairslist))
            print(cmd)
            cursor.execute(cmd)
            print(cmd2)
            cursor.execute(cmd2)
def insertPairwise2(db, cursor, pairslist, info, table, counter):
    for i in range(0, len(pairslist)):
        for j in range(0, len(pairslist)):
            print(i,j)
            cmd = "INSERT OR IGNORE INTO {} ({}, {}, {}, {},{}) VALUES ({},{}, {}, {}, {})".format(table, "tandemID","tandem1", "tandem2","info", "tandemcount", counter, q(pairslist[i]),q(pairslist[j]),q(info), len(pairslist))
            #cmd2 = "INSERT OR IGNORE INTO {} ({}, {}, {},{}) VALUES ({}, {}, {}, {})".format(table, "tandem1", "tandem2","info", "tandemcount", q(pairslist[j]),q(pairslist[i]),q(info),len(pairslist))
            print(cmd)
            cursor.execute(cmd)
            #print(cmd2)
            #cursor.execute(cmd2)
            
def readCsvIntoDB_tandem2(f,table):
    db = sqlite3.connect(DBNAME)
    c = db.cursor()
    counter=1
    try:
        infile = open(f,"r")
        reader = csv.reader(infile)
        for r in reader:
            if r[0].startswith(HEADERSTART):
                continue
            
            #print(r)
            tmpA1, tmpA2, tmpB1, tmpB2 = r[0],r[1],r[2],r[3]
            A1=tmpA1.split("||") #tandem duplicates are in the same field, separated by doube pipes
            A2=tmpA2.split("||")
            B1=tmpB1.split("||")
            B2=tmpB2.split("||")
            # insert pairwise into table, even with redundancy
            insertPairwise2(pairslist = A1, info = "A1", db=db, table=TANDEM, cursor=c, counter=counter) 
            insertPairwise2(pairslist = A2, info = "A2", db=db, table=TANDEM, cursor=c, counter=counter) 
            insertPairwise2(pairslist = B1, info = "B1", db=db, table=TANDEM, cursor=c,counter=counter) 
            insertPairwise2(pairslist = B2, info = "B2", db=db, table=TANDEM, cursor=c,counter=counter) 
            counter+=1
    finally:
        infile.close()
        db.commit()
        db.close()
    return(None)

def q(s):
    if s=="NULL":
        return(s)
    if s=="":
        return(q("NULL"))
    return('"'+s+'"')

def usage():
    print ("""
    ##############################
    #
    ##############################
    -i, --infile=evs900_zm.csv
    -h, --help
    [-t, --table=TABLENAME] default:"ZmTandemDup"
    [-n, --new=DBNAME or   # creates a new database DBNAME
    -o, --old=DBNAME ]     # inserts into existing database DBNAME
    
    """)
    sys.exit(2)

def main():
    global DBNAME
    global TANDEM
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
            TANDEM = a
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
            createTableZmTandemDup(TANDEM)
            readCsvIntoDB_tandem2(infile, TANDEM)

    if createNew:
        if checkExists():
            sys.stderr.write("{} already exists.\n".format(DBNAME))
            exit(1)
        else:
            print("#ok\n")
            createTableZmTandemDup(TANDEM)
            readCsvIntoDB_tandem2(infile, DUP)
            #createViews()
    #readCsv(infile)
    pass
####################################
if __name__ == "__main__":
    main()
####################################
