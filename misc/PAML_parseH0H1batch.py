import getopt
import imp
import sys
import os.path
import re
import glob
import fnmatch
from scipy.stats import chi2

BASEIND = 1
#BRANCHIND = 2
BRANCHIND = 4
#MODELIND = 3
MODELIND = 5
"""
Read H0 codeml output and H1 codeml output, write table with their information
and perform loglikelihood ratio test (write to table).
"""

def usage():
    print ("""
    #############################################
    # python PAML_parseH0H1batch -n H0.out -a H1.out -o OUT.csv
    ############################################
    files that should be compared must have a shared basename, eg:
    ORTHOMCL_ORTHOMCL1520.raxAutoMRbl.nwk.nl.M0_e1
    ORTHOMCL_ORTHOMCL1520.raxAutoMRbl.nwk.nl.M7_e1
    ORTHOMCL_ORTHOMCL1520.raxAutoMRbl.nwk.nl.M8_e1
    ORTHOMCL_ORTHOMCL1520.raxAutoMRbl.nwk.nl.M8a_e1
    ORTHOMCL_ORTHOMCL1520.raxAutoMRbl.nwk.nl.M1a_e1
    ORTHOMCL_ORTHOMCL1520.raxAutoMRbl.nwk.nl.M2a_e1
    ORTHOMCL_ORTHOMCL1520.raxAutoMRbl.nwk.19.BM_e1
    ORTHOMCL_ORTHOMCL1520.raxAutoMRbl.nwk.19.Ah0
    ORTHOMCL_ORTHOMCL1520.raxAutoMRbl.nwk.19.Ah1


    general options:
    -d, --dir=DIR            codeml output
   // -n, --H0=H0.codeml.out   codeml output for the null hypothesis
   // -a, --H1=H1.codeml.out   codeml output for the alternative
    [-o, --out=H0vsH1.csv otherwise stdout]
    -v, --verbose    print header
    -h, --help  prints this

    """)
    sys.exit(2)



class CODEMLParser(object):
    #lnL(ntime: 47  np: 51): -23761.923333      +0.000000
    def retrieveNumParameters(self, infile):
        with open(infile, 'r') as file:
            text = file.read()
        #search for np
        m = re.findall('lnL.*?np:(.*?)\).*?', text, re.DOTALL)
        if m:
            m = m[0].strip()
            m= int(m)
        else:
            m = None
        #print("np ", m)
        return(m)
    def retrieveLogLikelihood(self, infile):
        with open(infile, 'r') as file:
            text = file.read()
        m = re.findall('lnL.*?np:.*?:(.*?)\n', text)
        if m:
            m = m[0].split(" ")
            m = [n for n in m if n !=""]
            m = float(m[0])
        else:
            m = None
        return(m)
    def retrieveNtime(self, infile):
        with open(infile, 'r') as file:
            text = file.read()
        m = re.findall('lnL\(ntime:(.*?)np.*?', text, re.DOTALL)
        if m:
            m = m[0].strip()
            m = int(m)
        else:
            m = None
        #print("ntime: ",m)
        return(m)


    def getPositiveSitesNEB(self, infile):
        """Only for models that allow Positive Sites"""
        #Naive ...
        #Positive sites for foreground lineages Prob(w>1):
        res = {}
        res["significant"]=[]
        res["table"] = []
        with open(infile, 'r') as file:
            text = file.read()
        m = re.findall('Naive Empirical Bayes (.*?)Bayes Empirical Bayes.*?', text, re.DOTALL)
        if m:
            #print(m)
            m = m[0].split("\n")[3:]
            m = [x for x in m if x !=""]
            for c in m:
                c = c.split(" ")
                c = [x for x in c if x !="" and x]
                sig  = c[-1].count("*")
                c = [x.strip("*") for x in c]
                if sig > 0:
                    res["significant"].append(c)
                c.append(sig*"*")
                res["table"].append(infile.split(os.sep)[-1]+","+",".join(c))
        else:
            m = None
        res["table"]="\n".join(res["table"])
        return(res.copy())


    def getPositiveSitesBEB(self, infile):
        """Only for models That allow Positive Sites"""
        #Bayes Empirical Bayes
        #Positive
        pass
        res = {}
        res["significant"]=[]
        res["table"] = []
        with open(infile, 'r') as file:
            text = file.read()
        m = re.findall('Bayes Empirical Bayes (.*?)The.*?', text, re.DOTALL)
        if m:
            #print(m)
            m = m[0].split("\n")[3:]
            m = [x for x in m if x !=""]
            for c in m:
                c = c.split(" ")
                c = [x for x in c if x !="" and x]
                sig  = c[-1].count("*")
                c = [x.strip("*") for x in c]
                if sig > 0:
                    res["significant"].append(c)
                    c.append(sig*"*")
                if sig*"*"=="":
                    c.append("not significant")
                res["table"].append(infile.split(os.sep)[-1]+","+",".join(c))
        else:
            m = None
        res["table"]="\n".join(res["table"])
        return(res.copy())

    ####/done############################################

    def retriveTreeLength(self, infile):
        pass
    def retrieveTree(self, infile):
        tree = None
        pass
    def retrieveKappa(self, infile):
        #kappa (ts/tv) =  1.68285

        pass
    def retrieveOmegaForSiteClasses(self, infile):
        proportion = None
        background = None
        foreground = None

        pass

           #dN/dS (w) for site classes (K=4)
           #site class             0        1       2a       2b
           #proportion       0.75099  0.03808  0.20075  0.01018
           #background w     0.05475  1.00000  0.05475  1.00000
           #foreground w     0.05475  1.00000  1.00000  1.00000

         #  """
         # Note: Branch length is defined as number of nucleotide substitutions per codon (not per neucleotide site).
         # tree length =   8.67311
         # ((((((((5: 0.05225, 6: 0.05123): 0.00669, 4: 0.08560): 0.07238, 7: 0.09315): 0.06587, (8: 0.16465, 9: 0.16997): 0.05051): 0.20088, ((25: 0.30062, (24: 0.13763, 23: 0.17610): 0.17407): 0.19026, 22: 0.27528): 0.48826): 0.16448, (((21: 0.22150, 20: 0.19343): 0.11323, ((19: 0.13009, (18: 0.08349, 17: 0.04123): 0.07359): 0.09879, (16: 0.26060, 15: 0.20944): 0.05647): 0.10406): 0.25470, ((10: 0.26605, 11: 0.17738): 0.06104, (14: 0.10475, (13: 0.16188, 12: 0.06001): 0.08444): 0.11437): 0.36321): 0.31263): 1.44785, 3: 0.17434): 0.03935, 1: 0.41502, 2: 0.13029);
         # ((((((((GRMZM2G069542_T01: 0.05225, Sb04g008720.1: 0.05123): 0.00669, GRMZM2G074122_T01: 0.08560): 0.07238, Si016228m: 0.09315): 0.06587, (Bradi3g09210.1: 0.16465, LOC_Os02g14770.1: 0.16997): 0.05051): 0.20088, ((Si005789m: 0.30062, (GRMZM2G083841_T01: 0.13763, Sb10g021330.1: 0.17610): 0.17407): 0.19026, Bradi1g39167.1: 0.27528): 0.48826): 0.16448, (((Sb07g014960.1: 0.22150, LOC_Os08g27840.1: 0.19343): 0.11323, ((Si028826m: 0.13009, (GRMZM2G473001_T01: 0.08349, Sb02g021090.1: 0.04123): 0.07359): 0.09879, (LOC_Os09g14670.1: 0.26060, Bradi4g27910.1: 0.20944): 0.05647): 0.10406): 0.25470, ((Bradi2g50380.1: 0.26605, LOC_Os01g55350.1: 0.17738): 0.06104, (Si000184m: 0.10475, (GRMZM2G110714_T02: 0.16188, Sb03g035090.1: 0.06001): 0.08444): 0.11437): 0.36321): 0.31263): 1.44785, Si000160m: 0.17434): 0.03935, Bradi2g06620.1: 0.41502, LOC_Os01g11054.1: 0.13029);
         # Detailed output identifying parameters
         # """

            #tree length =   8.67311
def getBasename(string, delimiter=".", index = 1):
    return(string.split(delimiter)[index])
def calculatePvalue(h0, h1):
    res = {}
    res["pval"]= None
    res["sig"] = None
    res["stars"] = None
    res["np_H0"] = None
    res["np_H1"] = None
    res["lnL_H0"] = None
    res["lnL_H1"] = None
    res["df"] = None
    res["tablerow"]  = None
    cp = CODEMLParser()

    nph0 =    cp.retrieveNumParameters(h0)
    nph1 = cp.retrieveNumParameters(h1)
    lnLh0 = cp.retrieveLogLikelihood(h0)
    lnLh1 = cp.retrieveLogLikelihood(h1)
    print("h1", lnLh1)
    ntimeh0 = cp.retrieveNtime(h0)
    ntimeh1 = cp.retrieveNtime(h1)
    try:
        ts = 2 * (lnLh1 - lnLh0)
    except TypeError:
        print("??", h0,h1, lnLh1, lnLh0)
        sys.exit(2)
    df = nph1 -nph0
    pval =  1 - (chi2.cdf(ts,df))
    sig = 0
    if (pval < 0.05):
        sig = 1
    if (pval < 0.01):
        sig = 2
    if (pval < 0.001):
        sig = 3
    stars = sig*"*"
    if stars == "":
        stars = "not significant"
    ######## set dictionary ########
    res["np_H0"] =  nph0
    res["np_H1"] =  nph1
    res["lnL_H0"] = lnLh0
    res["lnL_H1"] = lnLh1
    res["df"] = df
    res["pval"] = pval
    res["stars"]  = stars
    res["sig"] = sig
    res["tablerow"]= "{},{},{},{},{},{},{},{}\n".format(h0.split(os.sep)[-1],h1.split(os.sep)[-1],pval, stars, lnLh0, lnLh1, nph0, nph1 )


    return(res.copy())
def main():
    out=None
    h0 = None
    h1 = None
    verbose = False
    printToStdout = True
    dir = None
    infiles = []
    bebstr =""
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "n:a:o:d:vh", ["H0=","H1=","out=","dir=","verbose","help"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-n", "--H0"):
            h0 = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-d", "--dir"):
            dir = a
        elif o in ("-a", "--H1"):
            h1 = a
        elif o in ("-o", "--out"):
            out=a
            outneb=out+".neb"
            outbeb=out+".beb"
            printToStdout = False
        elif o in ("v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    ########################################
    if h0 is None and dir is None:
        usage()
    if h1 is None and dir is None:
        usage()

    if dir:
        m0 = {}
        m0f = {}
        m1 = {}
        m2 = {}
        m7 = {}
        m8 = {}
        m8a = {}
        mAh0 = {}
        mAh1 = {}
        bm = {}
        branches = {}
        models = {}

        #fill the dictionaries
        cp = CODEMLParser()
        model = None
        branch = None
        base = None
        for fname in os.listdir(dir):
            try:
                base = getBasename(fname, ".",BASEIND)
                branch = getBasename(fname, ".",BRANCHIND)
                model = getBasename(fname, ".",MODELIND)

                #print(base,branch, model)
            except IndexError:
                pass
                #print(fname)
            if base not in branches:
                branches[base] = [branch]
                models [base] = [model]
            else:
                branches[base].append(branch)
                models[base].append(model)

            if fnmatch.fnmatch(fname, '*.M8a*'):
                if not base in m8a:
                    m8a[base]=[fname]
                else:
                    m8a[base].append(fname)
            elif fnmatch.fnmatch(fname, '*.M8*'):
                if not base in m8:
                    m8[base]=[fname]
                else:
                    m8[base].append(fname)
            elif fnmatch.fnmatch(fname, '*.M7*'):
                if not base in m7:
                    m7[base]=[fname]
                else:
                    m7[base].append(fname)
            elif fnmatch.fnmatch(fname, '*.M2*'):
                if not base in m2:
                    m2[base]=[fname]
                else:
                    m2[base].append(fname)
            elif fnmatch.fnmatch(fname, '*.M1*'):
                if not base in m1:
                    m1[base]=[fname]
                else:
                    m1[base].append(fname)

            elif fnmatch.fnmatch(fname, '*.M0_f*'):
                if not base in m0f:
                    m0f[base]=[fname]
                else:
                    m0f[base].append(fname)
            elif fnmatch.fnmatch(fname, '*.M0_e*'):
                if not base in m0:
                    m0[base]=[fname]
                else:
                    m0[base].append(fname)

            ################# BRANCH SPECIFIC #############3
            elif fnmatch.fnmatch(fname, '*.BM*'):
                #print(branch)
                if not base in bm:
                    bm[base]=[fname]
                else:
                    bm[base].append(fname)

            elif fnmatch.fnmatch(fname, '*.Ah0*'):
                print("AH0")
                print(branch)
                if not base in mAh0:
                    mAh0[base]=[fname]
                else:
                    mAh0[base].append(fname)
            elif fnmatch.fnmatch(fname, '*.Ah1*'):
                print("AH1")
                if not base in mAh1:
                    mAh1[base]=[fname]
                else:
                    mAh1[base].append(fname)
        print(mAh1, mAh0)
        try:
            for group in m0:
                for alt in m0[group]:
                    #test fixed vs estimated
                    for null in m0f[group]:
                        pval  = calculatePvalue(null, alt)
                        sys.stdout.write(pval["tablerow"])

                #for alt in bm[group]:
                #    for null in m0[group]:
                #        pval  = calculatePvalue(null, alt)
                #        sys.stdout.write(pval["tablerow"]+"\n")

                for alt in m2[group]:
                    for null in m1[group]:
                        #print(null, alt)
                        pval  = calculatePvalue(null, alt)
                        sys.stdout.write(pval["tablerow"]+"")
                for alt in m2[group]:
                    for null in m1[group]:
                        #print(null, alt)
                        pval  = calculatePvalue(null, alt)
                        sys.stdout.write(pval["tablerow"]+"")
                for alt in m8[group]:
                    for null in m7[group]:
                        #print(null, alt)
                        pval  = calculatePvalue(null, alt)
                        sys.stdout.write(pval["tablerow"]+"")
                for alt in m8[group]:
                    for null in m8a[group]:
                        #print(null, alt)
                        pval  = calculatePvalue(null, alt)
                        sys.stdout.write(pval["tablerow"]+"")
                for alt in mAh1[group]:
                    #tmpbranch = getBasename(alt, ".",4)
                    #print(alt, getBasename(alt, ".",4))
                    for null in mAh0[group]:
                        #print(null, alt,getBasename(null, ".",BRANCHIND), getBasename(alt, ".",BRANCHIND) )
                        if getBasename(null, ".",BRANCHIND) == getBasename(alt,".",BRANCHIND):
                            pval  = calculatePvalue(null, alt)
                            sys.stdout.write(pval["tablerow"]+"")
                            #SITES
                            sitesBEB = cp.getPositiveSitesBEB(alt)
                            if sitesBEB["table"] !="":
                                bebstr+=(sitesBEB["table"]+"\n")

                for alt in bm[group]:
                    for null in m0[group]:
                        pval  = calculatePvalue(null, alt)
                        sys.stdout.write(pval["tablerow"]+"")

            for group in mAh1:
                for alt in mAh1[group]:
                    #tmpbranch = getBasename(alt, ".",4)
                    print("getbasealt",alt, getBasename(alt, ".",4))
                    for null in mAh0[group]:
                        print("null",null, alt,getBasename(null, ".",BRANCHIND), getBasename(alt, ".",BRANCHIND), getBasename(null, ".",BRANCHIND).split(".")[0:-1] == getBasename(alt,".",BRANCHIND).split(".")[0:-1] )
                        if getBasename(null, ".",BRANCHIND).split(".")[0:-1]== getBasename(alt,".",BRANCHIND).split(".")[0:-1]:
                            pval  = calculatePvalue(null, alt)
                            print("pval",pval,null,alt)
                            sys.stdout.write(pval["tablerow"]+"")
                            #SITES
                            sitesBEB = cp.getPositiveSitesBEB(alt)
                            print(sitesBEB)
                            if sitesBEB["table"] !="":
                                bebstr+=(sitesBEB["table"]+"\n")
        except KeyError as e:
            with open("pamlbatchErr",'a') as errout:
                errout.write("Error w/ : "+str(e)+"\n")
            sys.stderr.write(str(e)+"\n")
        if outbeb:
            with open(outbeb, 'w') as outh:
                outh.write(bebstr)

    else:
        ########################################
        pval  = calculatePvalue(h0, h1)
        cp = CODEMLParser()
        sitesNEB = cp.getPositiveSitesNEB(h1)
        sitesBEB = cp.getPositiveSitesBEB(h1)

        if printToStdout:
            sys.stdout.write(pval["tablerow"])
            sys.stdout.write("\nNEB\n")
            sys.stdout.write(sitesNEB["table"])
            sys.stdout.write("\nBEB\n")
            sys.stdout.write(sitesBEB["table"])


        else:
            with open(out, 'w') as outh:
                outh.write(pval["tablerow"])
            #with open(outneb, 'w') as outh:
            #    outh.write(sitesNEB["table"])
            #    outh.write("\n")
            with open(outbeb, 'w') as outh:
                outh.write(sitesBEB["table"])
                outh.write("\n")

if __name__ == "__main__":
    main()

