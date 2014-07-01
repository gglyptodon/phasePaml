import getopt
import imp
import sys
import os.path
import re
import glob


"""
Read yn00 output for pairwise comparisons 
stored in single files in a directory and 
return a table for this directory.
"""

def usage():
    print ("""
    #########################################
    # python parseYN00.py -D dir -o outfile 
    #########################################
    
    general options:
    -D, --dir=DIR    directory with yn00 output files 
    -o, --out=OUTFILE    file result will be written to
    -p, --pattern=PATTERN    eg *.yn00_out
   [ -m, --method=METHOD    one of [all, LWL85, LWL85m, LPB93] ]
    -h, --help  prints this
     
    """)
    sys.exit(2)



class YN00Parser(object):
    import re
    #def __init__(self):
    #    self._patterns = {}
    #    self._res = {}
   
    def retrievePairsNG86(self, text, infile):
        
        res = {}
        res["file"] = infile
        text = text
        m2 = re.findall('ML pairwise comparison.\)\n\n(.*?)\n\np', text, re.DOTALL)
        #print(m2)
        m2string = "".join(m2)
        numGenes = len(m2string.split("\n"))
        #print(numGenes, "numGenes")
        for i in range(0, numGenes-1):
            #print(m2string.split("\n")[i].strip().split(" ")[0])
           
            for j in range(i+1, numGenes):
                res["GeneA"] =m2string.split("\n")[i].strip().strip().split(" ")[0].strip()
                #print(res["GeneA"])
                res["GeneB"] =m2string.split("\n")[j].strip().split(" ")[0].strip()
                #print(res["GeneB"])
                res["NG86_w"] = m2string.split("\n")[j].split(")")[i].split(" ")
                
                t=res["NG86_w"]
                t=[e for e in t if e !=""]
                res["NG86_w"] =t
                #print(res["NG86_w"])
                tmp = res["NG86_w"]
                res["NG86_w"] = tmp[-3:][0]
                res["NG86_dN"] = tmp[-3:][1].strip("(")
                res["NG86_dS"] = tmp[-3:][2]
                #print(res["NG86_w"])
                #print(res["NG86_dN"])
                #print(res["NG86_dS"])    
                #print(res)
                yield res
                
        yield(res)
    def writeNG86(self, text, infilename, outfilename):
        """" text is the whole paml output"""""
        written = set()
        outfile = open(outfilename,'a')
        warnfile =open(outfilename+".error",'a')
        print(infilename,"infile")
        #m = re.findall( '\(A\)(.*?)\(B\)', text, re.DOTALL)
        #text="".join(m)
        #print(text)
        #tableheaders = ["GeneA", "GeneB","NG86_w", "NG86_dN", "NG86_dS"]
        #outfile.write("file"+"\t"+"\t".join(tableheaders)+"\n")
        #print(m)
        for i in self.retrievePairsNG86(text, infilename):
           # print(i)
            try:
                ok = i["GeneA"]
            except KeyError,e:
                err = "WARNING: "+str(e)+" "+infilename+"\n"
                print(err)
                warnfile.write(err)
                continue
            try:    
                if (i["GeneA"], i["GeneB"]) not in written:
                        outfile.write("\t".join([infilename, i["GeneA"], i["GeneB"],i["NG86_w"], i["NG86_dN"], i["NG86_dS"]]))
                        outfile.write("\n")
                written.add((i["GeneA"], i["GeneB"] ))
            except TypeError, e:
                print "Type",e
            
        outfile.close()
    
    
   
    def retrievePairsCODEML(self, text, infilename):
            
            res = {}
            
            tableheaders = ["file","GeneA", "GeneB"]
                            #"LWL85_dS","LWL85_dN","LWL85_w",
                            #"LWL85_S","LWL85_N",
                            #"LWL85m_dS","LWL85m_dN",
                            #"LWL85m_w","LWL85m_rho",
                            #"LPB93_dS","LPB93_dN",
                            #"LPB93_w"]
            for i in tableheaders:
                res[i]=None
             
            res["file"] = infilename                
            content = text
            #m = re.findall( '\(C\)(.*?)$', content, re.DOTALL)
            m = re.findall( 'pairwise comparison, codon.*?\n(.*?)$', content, re.DOTALL)
            m="".join(m)
            #print(m)
            #m = re.findall( '2.13\n(.*?)$', m, re.DOTALL)
            #m="".join(m)
            #print("$$$$",m,"@@@@@")
            genePairs  = re.findall( '.*?\((.*?)\n', m)
            lnL = re.findall( 'lnL(.*?)\n', m)
            #print(lnL)
            parameters = re.findall('(t=.*?)\n', m)
            parameters2 = re.findall('(t=.*?)$', m)
            #print(parameters, "parameter")
            parameters = parameters +parameters2
            #print(len(parameters))
            #print(len(genePairs), genePairs)
            #print(re.findall( '\((.*?)\)', genePairs[1]))
            #L =  re.findall( '(L.*?)\nNs', m)
            #Ns = re.findall( '(Ns.*?)\nNv', m)
            #A = re.findall( '(A.*?)\nB', m)
            #LWL85 = re.findall( '(LWL85.*?)\nLWL85m', m)
            #LWL85m = re.findall( '(LWL85m.*?)\nLPB93', m)
            #LPB93 = re.findall( '(LPB93.*?)\n', m)
            #LPB93add = re.findall( '(LPB93.*?)$', m)
            #print(len(LWL85), len(LPB93 ))
            #LPB93=LPB93+LPB93add 
            #print(LPB93[26],(LPB93[0]))
            #print(LPB93[0])
            for i in range(0,len(genePairs)):
                #print(genePairs[i])
                #print(re.findall( '(.*?)\)', genePairs[i]))
                res["GeneA"]=re.findall( '(.*?)\)', genePairs[i])[0]
                res["GeneB"]=re.findall( '(.*?)\)', genePairs[i])[1]
                #print( res["GeneA"])#,res["GeneB"])
                res["GeneB"]=re.findall( '.*?\((.*?)$', res["GeneB"])[0]
                res["codeml_lnL"] = re.findall('=(.*?)$', lnL[i])[0]
                #res["t"] = 
                #print(re.findall('t=(.*?)S')
                try:
                    a = parameters[i]
                    #print(re.findall('t=(.*?)S',parameters[i]))
                    
                except IndexError:
                    pass
                res["codeml_t"] = re.findall('t=(.*?)S', parameters[i])[0].strip()
                res["codeml_S"] =  re.findall('S=(.*?)N', parameters[i])[0].strip()
                res["codeml_N"] = re.findall('N=(.*?)dN', parameters[i])[0].strip()
                res["codeml_dN/dS"]= re.findall('dN/dS=(.*?)dN', parameters[i])[0].strip()
                res["codeml_dN"]= re.findall('dN =(.*?)dS', parameters[i])[0]
                res["codeml_dS"]= re.findall('dS =(.*?)$', parameters[i])[0]
                #print[res]
               # print(res["lnL"], lnL[i])
                #print(genePairs[i])
                
                
            #logLikelihood:
            
                #res["GeneB"]=re.findall( '\((.*?)\)', genePairs[i])[1]
               # print(LWL85[i])
                
             #   res["LWL85_dS"]=re.findall('dS =(.*?)dN', LWL85[i])[0].strip()
             #   res["LWL85_dN"]=re.findall('dN =(.*?) w', LWL85[i])[0].strip()
             #   res["LWL85_w"]=re.findall('w =(.*?)S', LWL85[i])[0].strip()
             #   res["LWL85_S"]=re.findall(' S =(.*?)N', LWL85[i])[0].strip()
             #   res["LWL85_N"]=re.findall(' N =(.*?)$', LWL85[i])[0].strip()
                #print(res["LWL85_N"])
                
                #print(LWL85m[i])
                 
              #  res["LWL85m_dS"]=re.findall('dS =(.*?)dN', LWL85m[i])[0].strip()
              #  res["LWL85m_dN"]=re.findall('dN =(.*?) w', LWL85m[i])[0].strip()
              #  res["LWL85m_w"]=re.findall('w =(.*?)S', LWL85m[i])[0].strip()
              #  res["LWL85m_S"]=re.findall(' S =(.*?)N', LWL85m[i])[0].strip()
              #  res["LWL85m_N"]=re.findall(' N =(.*?) \(', LWL85m[i])[0].strip()
              #  res["LWL85m_rho"]=re.findall(' \(rho = (.*?)\)', LWL85m[i])[0].strip()
                
                #print(res["LWL85m_N"])
                #print(LPB93[i])
                #print(LPB93[i])
               # res["LPB93_dS"]=re.findall('dS =(.*?)dN', LPB93[i])[0].strip()
               # res["LPB93_dN"]=re.findall('dN =(.*?) w',LPB93[i])[0].strip()
               # res["LPB93_w"]=re.findall('w =(.*?)$', LPB93[i])[0].strip()
               # print(res["LPB93_dN"],LPB93[i] )
                #print(res["LPB93_w"])
                #print("LWL85_S",res["LWL85_S"])
                #print(res["LPB93_w"],LPB93[i] )
                #res["LWL85_dN"]=re.findall( 'dN = (.*?)', LWL85[i])
                #res["LWL85_w"]=re.findall( 'w = (.*?)', LWL85[i])
                #res["LWL85_S"]=re.findall( 'S = (.*?)', LWL85[i])
                #res["LWL85_N"]=re.findall( 'N = (.*?)', LWL85[i])
                #print(res)
                #---------------------------------------- print(i, genePairs[i])
                #------------------------------------------------- println(L[i])
                #------------------------------------------------ println(Ns[i])
                #------------------------------------------------- println(A[i])
                #--------------------------------------------- println(LWL85[i])
                #-------------------------------------------- println(LWL85m[i])
                #--------------------------------------------- println(LPB93[i])
#------------------------------------------------------------------------------ 
            #---------------------------------------------------------- print(m)
                #print(res)
                yield(res)
            print("#########\n")
            
            
    def writeCODEML(self, text, infilename, outfilename):
        """" text is the whole paml output"""""
        written = set()
        outfile = open(outfilename,'a')
       
        #m = re.findall( '\(A\)(.*?)\(B\)', text, re.DOTALL)
        #text="".join(m)
        tableheaders = ["GeneA", "GeneB", "codeml_dN","codeml_dS",
                         "codeml_dN/dS","codeml_N","codeml_S",
                          "codeml_t", "codeml_lnL"  ]
        #                    "LWL85_dS","LWL85_dN","LWL85_w",
        #                    "LWL85_S","LWL85_N",
        #                    "LWL85m_dS","LWL85m_dN",
        #                    "LWL85m_w","LWL85m_rho",
        #                    "LPB93_dS","LPB93_dN",
         #                   "LPB93_w"]
        
         
        
           
        #outfile.write("file"+"\t"+"\t".join(tableheaders)+"\n")
        for i in self.retrievePairsCODEML(text, infilename):
           # print(i)
            tmp = [i[j] for j in tableheaders]
            #print(tmp)
            tmp = [str(i) for i in tmp]
           #print(i, )
            tableline = infilename+"\t"+"\t".join(tmp)+"\n"
            outfile.write(tableline)
            #if (i["GeneA"], i["GeneB"]) not in written:
            #    outfile.write("\t".join([infilename, i["GeneA"], i["GeneB"],]]))
            #    outfile.write("\n")
            #written.add((i["GeneA"], i["GeneB"] ))
        outfile.close()
    def parseFile(self,file, outfile):
        #fileh = open(file,'r')
        with open(file, 'r') as content_file:
            content = content_file.read()
        res ={}
        res["file"] = None
        res["GeneA"] = None
        res["GeneB"] = None
        ###################
        res["NG86_dS"] = None
        res["NG86_dN"] = None 
        res["NG86_w"] = None
        #####################
        #res["YN00_S"] = None
        #res["YN00_N"] = None
        #res["YN00_t"] = None
        #res["YN00_kappa"] = None
        #res["YN00_w"]= None
        #res["YN00_dN+SE"] = None
        #res["YN00_dN-SE"] = None
        #res["YN00_dS+SE"] = None
        #res["YN00_dS-SE"] = None
        ####################
        #res["LWL85_dS"]=None
        #res["LWL85_dN"]=None
        #res["LWL85_w"]=None
        ####################
        #res["LWL85m_dS"]=None
        #res["LWL85m_dN"]=None
        #res["LWL85m_w"]=None
        #res["LWL85m_rho"]=None
        #####################
        #res["LPB93_dS"]=None
        #res["LPB93_dN"]=None
        #res["LPB93_w"]=None
        
        
        #filename
        for i in content.split("\n"):
          
            #if i.startswith("YN00"):
            res["file"] =file.split(os.sep)[-1]
            #i.strip().split(" ")[1]
        
        #############NG:
       
                
        #NG86
        
        #m = re.findall( '\(A\)(.*?)\(B\)', content, re.DOTALL)
        #m="".join(m)
        #print(content)
        self.writeNG86(content, res["file"], outfile+"_NG86.tsv")
       
        #yn00
       # self.retrievePairsYN00(content)
       #
       
       
      
        #self.writeYN00(content, res["file"], outfile+"_YN00.tsv")
      
        #m = re.findall( '\(B\)(.*?)\(C\)', content, re.DOTALL)
        #m="".join(m)
        #print(m)
        
        #LWL85, LPB93 & LWLm
        #print[res["file"]]
        self.writeCODEML(content, res["file"], outfile+"_CODEML.tsv")
        self.retrievePairsCODEML(content, res["file"])
        #self.writeNG86(m, res["file"], outfile+"_NG86.tsv")
        
        #print(m)
        #for r in self.retrievePairsNG86(m):
        #    #print("RES ",r)
        #    #print("\n")
        #    res.update(r)
        #    print(res)
            
            #else:
                #print(i)
        #print(res) 
def println(s):
    print(s,"\n")
def writeSummary(out, filesInDir):
    fileDct = {}
    outh = open(out,'w')
    header = "\t".join(["file","GeneA", "GeneB","codeml_dN","NG86_dN",
                        "codeml_dS","NG86_dS","codeml_dN/dS","NG86_w",
                        "codeml_N","codeml_S","codeml_t","codeml_lnL"])+"\n"
    print(header)
    outh.write(header)
    outh.close()
    outh = open(out,'a')
    for f in filesInDir:
        print(f)
        with open(f, 'r') as content_file:
            content = content_file.read()
            for i in YN00Parser().retrievePairsCODEML(content,f.split(os.sep)[-1]):
                print((i["GeneA"],i["GeneB"]),i["file"])
                fileDct[(i["GeneA"],i["GeneB"])] = dict()
                fileDct[(i["GeneA"],i["GeneB"])]["codeml"] = i.copy()
                #print("QQQQQQ",fileDct[(i["GeneA"],i["GeneB"])])
                ##for k in (fileDct.keys()):
                #    println(k)
            for i in YN00Parser().retrievePairsNG86(content,f.split(os.sep)[-1]):
                
                try:
                    #print(i,"III",i["GeneA"])
                    print((i["GeneA"],i["GeneB"]),"AB", (i["GeneB"],i["GeneA"]) in fileDct,i["file"])
                    if ((i["GeneA"],i["GeneB"])in fileDct ):
                        #print("already there","AB") 
                        fileDct[(i["GeneA"],i["GeneB"])]["ng86"] = i.copy()
                        #print(fileDct[(i["GeneA"],i["GeneB"])])
                    elif ((i["GeneB"],i["GeneA"]) in fileDct):
                        #print("already there","BA") 
                        #print("RRR",fileDct[(i["GeneB"],i["GeneA"])])
                        fileDct[(i["GeneB"],i["GeneA"])]["ng86"] =i.copy()
                        #print("RRR",fileDct[(i["GeneA"],i["GeneB"])])
                        #print(fileDct[(i["GeneA"],i["GeneB"])]["ng86"])
                    else:
                        print("sth missing")
                        
                except KeyError,e:
                    print e
                    pass
                #if (i["GeneA"],i["GeneB"]) in fileDct:
                #    print("already there") 
    print(fileDct)
    for k,v in (fileDct.items()):
                #print(k,v["codeml"])
        print("\n")
                #print(k,v["ng86"])
                #print("\n")
        try: 
            s = "\t".join([v["codeml"]["file"],v["codeml"]["GeneA"],v["codeml"]["GeneB"], 
                           v["codeml"]["codeml_dN"],v["ng86"]["NG86_dN"],
                           v["codeml"]["codeml_dS"],v["ng86"]["NG86_dS"],
                           v["codeml"]["codeml_dN/dS"],v["ng86"]["NG86_w"],
                           v["codeml"]["codeml_N"],v["codeml"]["codeml_S"],
                           v["codeml"]["codeml_t"],v["codeml"]["codeml_lnL"]])+"\n"
            outh.write(s)
        except KeyError,k:
                print(v["codeml"]["file"])
                outh.write("#"+v["codeml"]["file"]+"\n")
    outh.close()
                
def main():
    out=None
    dir = None
    pattern = None
    what="all"
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "D:o:p:m:h", ["dir=","out=","pattern=","method=","help"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-D", "--dir"):
            dir = a  
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-p", "--pattern"):
            pattern = a
        elif o in ("-o", "--out"):
            out=a
        elif o in ("-m", "--method"):
            what = a
        else:
            assert False, "unhandled option"
    
    ########################################
    if dir is None:
        usage()
    if out is None:
        usage()
    if pattern is None:
        usage()
    ########################################
    filesInDir= []
    print(dir)
    #header = mkHeader(what=what)
    if dir:
        dir = dir+os.sep+pattern
        for fname in glob.glob(dir):
            #print(fname)
            filesInDir.append(fname)
    print(filesInDir)
    outhandle = open(out+"_NG86.tsv",'w')
    tableheaders = ["GeneA", "GeneB","NG86_w", "NG86_dN", "NG86_dS"]
    outhandle.write("file"+"\t"+"\t".join(tableheaders)+"\n")
    outhandle.close()
    tableheaders = ["GeneA", "GeneB","YN00_S","YN00_N", "YN00_t" ,
                       "YN00_kappa" ,"YN00_w",
                       "YN00_dN+SE", "YN00_dN-SE",
                       "YN00_dS+SE", "YN00_dS-SE"] 
    #outhandle = open(out+"_YN00.tsv",'w')
    #outhandle.write("file"+"\t"+"\t".join(tableheaders)+"\n")
    #outhandle.close()
    
    outhandle = open(out+"_CODEML.tsv",'w')
    
    
    tableheaders = ["GeneA", "GeneB", "codeml_dN","codeml_dS", "codeml_dN/dS","codeml_N","codeml_S", "codeml_t", "codeml_lnL"  ]

    #                        "LWL85_dS","LWL85_dN","LWL85_w",
    #                        "LWL85_S","LWL85_N",
    #                        "LWL85m_dS","LWL85m_dN",
    #                        "LWL85m_w","LWL85m_rho",
    #                        "LPB93_dS","LPB93_dN",
    #                        "LPB93_w"]
        
    outhandle.write("file"+"\t"+"\t".join(tableheaders)+"\n")
    outhandle.close()
    
    ynp=YN00Parser()
    writeSummary(out+"_summary.tsv",filesInDir)
    for f in filesInDir:
        #print(f)
        ynp.parseFile(f,out)
if __name__ == "__main__": 
    main()

