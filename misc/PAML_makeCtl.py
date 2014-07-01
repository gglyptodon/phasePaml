#########################
# last update:
# Tue 15 Apr 2014 04:09:45 PM CEST
# [JMass]
#########################
import getopt
import imp
import sys
import os.path
import re
import glob
from ete2 import Tree
from ete2 import EvolTree
from ete2 import TreeStyle
from ete2 import NodeStyle
from ete2 import ImgFace
from ete2 import TextFace
import subprocess
import dendropy
import re

#from ete2 import random_color, COLOR_SCHEMES, SVG_COLORS

import random

#print "#%x" % x
MAJORITY=0.5
MAX_PARENT = 4
# generate MAX_PARENT random colors
NODE_COLORS = [ random.randint(0, 16777215) for i in range(0,MAX_PARENT)]

print(NODE_COLORS)
"""
Produces control files for codeml with the given trees, alignment and model.
It will automatically produce ctl for contrasting hypotheses
and different omega start values.

Example: -t tree1.nwk -a al1.paml -mM8 -n9 -A
will also create ctl files for M7 and M8a.
"""

def usage():
    print ("""
    #############################################
    # example usage:
    # python PAML_makeCtl -t tree.nwk -a alignment.paml --model="M1a,M2a,M0,M7,M8a,Ah0,Ah1" --nodes="10,11,13"
    #or:
    # python PAML_makeCtl -t tree.nwk -a alignment.paml --model="M1a,M2a,M0,M7,M8a,Ah0,Ah1" --regex='GRM.*'
    ############################################

    general options:
    -t, --tree=tree.nwk input tree file in newick format
     or
    -T, --bootstrapTree=TREE.raxmlbootstrap


    [-a, --alignment=MSA.paml]    alignment in paml format
    [-n, --nodes=ListOfNodes ]    comma-separated list of nodes
                                that should be labelled
                                if no nodes are listed,
                                the tree will be traversed
                                and node ids are shown
    [-r, --regex=REGEX]         alternative to '-n', label all nodes and their predecessors that match REGEX and their

    [-m, --model=MODEL]         MODEL can be one of ["M0","M1a","M2a",
                                "M7","M8","M8a", "Ah0", "Ah1","BM"] or a comma-separated list of those
                                for details run with "-M"
    //[-o, --out=OUT]             output prefix (default: treefile specified with -t)
    -A, --writeOther            also write ctl for complementary hypotheses
    -M, --modelhelp             prints model descriptions
    -h, --help                  prints this

    If no nodes are specified, nodes and their labels will be printed to stdout.
    If no nodes are specified but a model is specified, the program will assume that the tree is already labelled.


    """)
    sys.exit(2)

def ctlTemplate(seqfile = None, treefile = None, outfile = None, runmode = "0",\
                seqtype = "1", model = None, fix_kappa = "0",kappa= "2",\
                fix_omega= None, omega =None, NSsites=None,CodonFreq="2"):

    seqfile = "seqfile = {}\n".format(seqfile)
    treefile = "treefile = {}\n".format(treefile)
    outfile = "outfile = {}\n".format(outfile)
    runmode = "runmode = {}\n".format(runmode)
    seqtype = "seqtype = {}\n".format(str(seqtype))
    CodonFreq = "CodonFreq = {}\n".format(CodonFreq)
    model = "model = {}\n".format(model)
    NSsites = "NSsites = {}\n".format(NSsites)
    fix_kappa ="fix_kappa = {}\n".format(str(fix_kappa))
    #0 * 1: kappa fixed, 0: kappa to be estimated
    kappa = "kappa = {}\n".format(str(kappa))
    #2 * initial or fixed kappa
    fix_omega="fix_omega = {}\n".format(fix_omega) #* 1: fixed, 0: estimate
    omega = "omega = {}\n".format(str(omega)) #1 * initial or fixed omega, for codons or codon-based AAs


    runmodeExplanation = """
    *********************
    * 0: user tree;
    * 1: semi-automatic;
    * 2: automatic
    * 3: StepwiseAddition;
    * (4,5):PerturbationNNI;
    * -2: pairwise
    *********************
    """
    seqtypeExplanation = """
    * 1:codons; 2:AAs; 3:codons-->AAs
    """
    CodonFreqExplanation = """
    *0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
    """
    modelExplanation = """
    **********************
    * models for codons:
    * 0:one,
    * 1:b (free ratio),
    * 2:2 or more dN/dS ratios for branches
    * When model = 2, you have to group branches on the tree into branch groups
    * using the symbols # or $ in the
    **********************
    """
    NSsitesExplanation = """
    **********************
    * 0:one w;1:neutral;
    * 2:selection; 3:discrete;
    * 4:freqs; 5:gamma;
    * 6:2gamma;7:beta;
    * 8:beta&w;9:beta&gamma;
    * 10:beta&gamma+1; 11:beta&normal>1;
    * 12:0&2normal>1; 13:3normal>0
    **********************
    """
    stuff = """
    *icode = 0 * 0:universal code; 1:mammalian mt; 2-11: ...
    *Mgene = 0 * 0:rates, 1:separate
    ***************************************
    *fix_alpha = 1 *0: estimate gamma shape parameter; 1: fix it at alpha
    *alpha = 0.1 *initial or fixed alpha, 0:infinity (constant rate)
    *Malpha = 0 * different alphas for genes
    *ncatG = 8 *# of categories in dG of NSsites models
    *clock = 0 * 0:no clock, 1:clock; 2:local clock; 3:TipDate
getSE = 0 * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0 * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
Small_Diff = .5e-6
cleandata = 1 * remove sites with ambiguity data (1:yes, 0:no)
method = 0 * 0: simultaneous; 1: one branch at a time
    **************
    **************

    """
    template = seqfile +treefile + outfile
    template += runmode + seqtype + CodonFreq +model
    template += NSsites + fix_kappa + kappa
    template += fix_omega + omega
    template += stuff

    return(template)

def generateCtl(model, treefile, seqfile, outfile, generateOther = False):
    print("GENERATE CTL")
    if not outfile:
        outfile = treefile
    ######### Basic model ########################
    if (model.upper() == "M0"):
        outsuff0 = ".M0_e05"
        outsuff1 = ".M0_e10"
        outsuff2 = ".M0_f10"

        ctl0 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile = outfile+outsuff0, model= "0",NSsites = "0",fix_omega = "0",omega = "0.5")
        ctl1 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile = outfile+outsuff1,model= "0",NSsites = "0",fix_omega = "0",omega = "1")
        ctl2 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile = outfile+outsuff2,model= "0",NSsites = "0",fix_omega = "1",omega = "1")
        ################# write control files ###############
        with open(outfile+outsuff0+".ctl", 'w') as out:
            out.write(ctl0)
        with open(outfile+outsuff1+".ctl", 'w') as out:
            out.write(ctl1)
        with open(outfile+outsuff2+".ctl", 'w') as out:
            out.write(ctl2)
    ########## /Basic model ######################

    ##########  Branch models#####################
    elif (model.upper() =="FR"):
        outsuff0 = ".FR_e1"
        ctl0 = ctlTemplate(treefile = treefile, seqfile = seqfile, outfile = outfile+outsuff0, model= "1", NSsites = "0",fix_omega = 0, omega = 1)
        #free ratio model  - use is discouraged
        ################# write control files ###############
        with open(outfile+outsuff0+".ctl", 'w') as out:
            out.write(ctl0)

    elif (model.upper() =="BM"):
        outsuff0 = ".BM_f1"
        outsuff1 = ".BM_e1"
        ctl0 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile = outfile+outsuff0, model= "2", NSsites = "0",fix_omega = 1, omega = 1)
        ctl1 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile = outfile+outsuff1, model= "2", NSsites = "0",fix_omega = 0, omega = 1)
        ################# write control files ###############
        with open(outfile+outsuff0+".ctl", 'w') as out:
            out.write(ctl0)
        with open(outfile+outsuff1+".ctl", 'w') as out:
            out.write(ctl1)

    ############ /Branch models ######################

    ############  Site models   ######################
    # M1a
    elif (model.upper() == "M1A"):
        #M1aNearlyNeutral
        outsuff0 = ".M1a_e1"
        ctl0 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile= outfile+outsuff0, model = 0, NSsites = 1, fix_omega = 0, omega = 1)
        if generateOther:
            generateCtl("M2A", treefile, seqfile, outfile, False)
        ################# write control files ###############
        with open(outfile+outsuff0+".ctl", 'w') as out:
            out.write(ctl0)

    #M2a
    elif (model.upper() == "M2A"):
        #M2a (PositiveSelection)
        outsuff0 = ".M2a_e1"
        ctl0 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile= outfile+outsuff0, model = 0, NSsites = 2, fix_omega = 0, omega = 1)
        if generateOther:
            generateCtl("M1A", treefile, seqfile, outfile, False)
        ################# write control files ###############
        with open(outfile+outsuff0+".ctl", 'w') as out:
            out.write(ctl0)

    #M7
    elif (model.upper() == "M7"):
       #M7 (beta)
       outsuff0 = ".M7_e1"
       ctl0 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile= outfile+outsuff0, model = 0, NSsites = 7, fix_omega = 0, omega = 1)
       if generateOther:
            generateCtl("M8", treefile, seqfile, outfile, False)
       ################# write control files ###############
       with open(outfile+outsuff0+".ctl", 'w') as out:
            out.write(ctl0)
    #M8
    elif (model.upper() == "M8"):
        outsuff0 = ".M8_e1"
        ctl0 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile= outfile+outsuff0, model = 0, NSsites = 8, fix_omega = 0, omega = 1)
        if generateOther:
            generateCtl("M7", treefile, seqfile, outfile, False)
            generateCtl("M8A", treefile, seqfile, outfile, False)
        ################# write control files ###############
        with open(outfile+outsuff0+".ctl", 'w') as out:
            out.write(ctl0)
    #M8a
    elif (model.upper() == "M8A"):
        outsuff0 = ".M8a"
        ctl0 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile= outfile+outsuff0, model = 0, NSsites = 8, fix_omega = 1, omega = 1)
        if generateOther:
            generateCtl("M8", treefile, seqfile, outfile, False)
        ################# write control files ###############
        with open(outfile+outsuff0+".ctl", 'w') as out:
            out.write(ctl0)

    #BranchSite models
    #modified model A Null Hypothesis
    elif(model.upper() == "AH0"):
        outsuff0 = ".Ah0"
        #null
        ctl0 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile= outfile+outsuff0, model = 2, NSsites = 2, fix_omega = 1, omega = 1)
        if generateOther:
            generateCtl("AH1", treefile, seqfile, outfile, False)
        ################# write control files ###############
        with open(outfile+outsuff0+".ctl", 'w') as out:
            out.write(ctl0)
    #modified model A Alternative
    elif(model.upper() == "AH1"):
        outsuff0 = ".Ah1"
        ctl0 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile= outfile+outsuff0, model = 2, NSsites = 2, fix_omega = 0, omega = 1.5)
        if generateOther:
            generateCtl("AH0", treefile, seqfile, outfile, False)
        ################# write control files ###############
        with open(outfile+outsuff0+".ctl", 'w') as out:
            out.write(ctl0)

    #PAML manual:
    #To calculate the p value based on this mixture distribution,
    #you calculate the p value using X_1^2, and then divide it by 2.
    #Suppose your 2DELTAl= 2.71, you will get 0.10 from X_1^2 ,
    #the the p value for the mixture is 0.10/2 = 5%.
    #We recommend that you use X_1^2 (with critical values 3.84 and 5.99)
    #instead of the mixture to guide against violations of model assumptions.




def readTree(tree):
    t = Tree(tree)
    #print (t.get_ascii(attributes=["name", "dist", "size"]))
    #print (t.dist)
    #print(t.write(format=9))
    with open(tree+".nl","w") as tree_nolabel:
        tree_nolabel.write(t.write(format=9))
    return(t.write(format=9))

def showAlignmentWithTree(tree,alignment):
    print(tree)
    t = EvolTree(tree, alignment,alg_format="paml")
    nsFG = NodeStyle()
    nsFG["fgcolor"] = "darkgreen"
    nsFG["size"] = 15
    #print(t)
    #t.run_model ('fb.example')
   # t.show()

    t.link_to_alignment(alignment, alg_format="paml")
    for node in t.traverse():
        print(node)
        #if (node.name.split("_")[0] in GENES):
        #    print(node.name, node.node_id)
        #    if (node.name.split("_")[0] == "GRMZM2G083841"):
        #        node.add_face(ImgFace("83841.1.png", height=50, width=50), column=1, position="aligned")
        #    if (node.name.split("_")[0] == "GRMZM2G473001"):
        #        node.add_face(ImgFace("473001.png", height=50, width=50), column=1, position="aligned")
        node.add_face(TextFace(str(node.node_id)),column=0)
            #node.add_face(ImgFace("tux.png", height=50), column=1)
        #    node.set_style(nsFG)
        leaves  = node.get_leaf_names()
        #leaves = [l for l in leaves if l.split("_")[0] in GENES ]
        if leaves !=[]:
            print(node.name, node.node_id)
        #print(node.node_id)
    #t.mark_tree([8], marks=["#1"])
    #print(t.write())
    #print(alignment)
    #print(t)
    #t.show() #layout=evol_clean_layout)

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = True
    t.show(tree_style=ts)
    t.render("tree.pdf", tree_style=ts)


def showAlignmentWithTreeOld(unlabelledTreeAsString,alignment):
    t = EvolTree(unlabelledTreeAsString, alignment,alg_format="paml")
    t.link_to_alignment(alignment, alg_format="paml")
    for node in t.traverse():
        print(node)
        print(node.node_id, node.name)
    print(t.write())
    #print(t)
    t.show() #layout=evol_clean_layout)

def labelForPaml(unlabelledTreeAsString,listOfNodes, tree):
    t = EvolTree(unlabelledTreeAsString)
    marks = []
    count = 1
    for i in listOfNodes:
        marks.append("#"+str(count))
        count+=1
    t.mark_tree(listOfNodes, marks=marks)
    print(t.write())
    outfile = tree+"."+"_".join(listOfNodes)
    with open(outfile, 'w') as out:
        out.write(t.write())
    return(outfile)
def rgb(c):
    split = (c[0:2], c[2:4], c[4:6])
    return [int(x, 16) for x in split]
def labelForPamlRegex(unlabelledTreeAsString, regex, tree):
    pattern = re.compile(regex)
    t = EvolTree(unlabelledTreeAsString)
    marks = []
    count = 1
    outfiles = []
    #nsFG=TreeStyle()
    nsMatch = NodeStyle() #match
    nsMatch["fgcolor"] = "blue"
    nsMatch["size"] = 10
    nsBG = NodeStyle()
    nsBG["fgcolor"] = "black"
    nsBG["size"] = 0
    nsFG = []

    tolabelreg = []
    for i in range(0,MAX_PARENT):
        nsFG.append( NodeStyle())
        nsFG[i]["size"] = 10
        nsFG[i]["fgcolor"] = NODE_COLORS[i]

    isroot=True
    for node in t.traverse():
        node.set_style(nsBG)
        if node.is_root():
            print("root")
            node.unroot()
            node._support = None

    for node in t.get_descendants():
        node.add_face(TextFace(node.node_id), column=0)

    #traverse and match
    for node in t.traverse():
        if re.match(pattern, node.name):
            print("MATCH", node.name, node.node_id)
            node.set_style(nsMatch)
            n = node
            try:
                for i in range(0,MAX_PARENT):

                    n = n.up
                    n.set_style(nsFG[i])
                    marks.append("#"+str(count))
                    print(count)
                    t.mark_tree([str(count)], marks=marks)
                    #just label everything with #1
                    print(t.write())

                    tolabelreg.append(str(n.node_id))

                    outfile = tree+"."+"_".join([str(n.node_id)])
                    with open(outfile, 'w') as out:
                        out.write(t.write())
                    outfiles.append(outfile)

            except AttributeError:
                pass

            marks.append("#"+str(count))
            print(count)
            t.mark_tree([str(count)], marks=marks)
            #just label everything with #1
            print(t.write())

            outfile = tree+"."+"_".join([str(node.node_id)])
            with open(outfile, 'w') as out:
                out.write(t.write())
                outfiles.append(outfile)
    #t.show()
    t.render(tree+".png")
    return(outfiles, tolabelreg)
def showTreeNodesOld(unlabelledTreeAsString):
    t = EvolTree(unlabelledTreeAsString)
    for node in t.traverse():
        print(node)
        print(node.node_id, node.name)
def showTreeNodes(unlabelledTreeAsString):
    t = EvolTree(unlabelledTreeAsString)
    for node in t.traverse():
        #print(node)
        #print(node.node_id, node.name)
        #if (node.name.split("_")[0] in GENES):
        #    print(node.name, node.node_id)
        leaves  = node.get_leaf_names()
        #leaves = [l for l in leaves if l.split("_")[0] in GENES ]
        if leaves !=[]:
            print(node)
            print(leaves, node.name, node.node_id)
            print("\n")

def modelhelp():
    print ("""
    ############################################
    # MODELS
    ############################################
    modelID:    Description:             vs:
      BM     Branch model, one fg branch M0, BM with w:=0
      M0          one-ratio model;       ?
    [ M1 ]    obsolete
      M1a    branch nearly neutral       M2a
    [ M2 ]    obsolete
      M2a    branch pos. sel.            M1a
      M7     branch beta                 M8
      M8     branch beta+omega>1,alt     M7, M8a
      M8a    branch beta+omega=1,null    M8
      Ah0    branch-site                 Ah1
      Ah1    branch-site                 Ah0
    [ B ]    not yet implemented
    [ C ]    not yet implemented
    """)
    sys.exit(2)




def main():
    tree = None
    tolabel = None
    alignment = None
    model = None
    writeOther = False
    bstree =False
    regex = False

    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "t:T:a:n:m:r:hAM",
                                        ["tree=","bootstrapTree=","alignment=",
                                         "node=","model=","regex=","help", "writeOther","modelhelp"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-t", "--tree"):
            tree = a
        elif o in ("-a", "--alignment"):
            alignment = a
        elif o in ("-n", "--nodes"):
            tolabel = a.split(",")
            print("debug", tolabel)
        elif o in ("-m", "--model"):
            model = a
        #elif o in ("-A", "--writeOther"):
            #writeOther = True
        elif o in ("-M", "--modelhelp"):
            modelhelp()
        elif o in ("-T", "--bootstrapTree"):
            bstree = a
        elif o in ("-r","--regex"):
            regex = a

        elif o in ("-h", "--help"):
            usage()

        else:
            assert False, "unhandled option"

    ########################################
    if bstree:
        res = ""
        nex = False
#########################
# last update:
# Tue 15 Apr 2014 10:53:54 AM CEST
# [JMass]
#########################
        p = subprocess.Popen('sumtrees.py --unrooted --min-clade-freq='+str(MAJORITY)+" --percentages --decimals=0 "+bstree, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
             if line.startswith("#NEXUS") or nex:
                 nex = True
                 res+= line
                 retval = p.wait()
        with open(".tmptree0123.tmp","w") as tmp:
            tmp.write(res)
        test = dendropy.TreeList.get_from_path(".tmptree0123.tmp", "nexus")

        with open(bstree+".nwk", 'w') as out:
            out.write(str(test))

        tree = bstree+".nwk"
        print(tree)

    if tree is None:
        usage()
    ########################################
    unlabelledTree = readTree(tree)
    #print(unlabelledTree)
    print(model, tree, alignment)

    ##
    if regex:
        res = labelForPamlRegex(unlabelledTree, regex, tree)
        tolabel = res[1]
        print(res,"RES")
    if model and tree and alignment:
        #model is a list:
        modellist = model.split(",")
        modellist  = [m.strip() for m in modellist]
        #tag trees and return their filenames
        taggedTrees = []
        if tolabel:
            tolabel = [l.strip() for l in tolabel]
            for n in tolabel:
                print("debug",n)
                taggedTrees.append(labelForPaml(unlabelledTree,[n], tree))
        #loop through new tagged trees:
        for t in taggedTrees:
            for model in modellist:
                print(model.upper())
                if model.upper() in ["M0","M1A","M2A","M7","M8","M8A"]:
                    print(unlabelledTree)
                    if (tree.endswith("nl")):
                        generateCtl(model=model, treefile = tree+"", seqfile=alignment, outfile=None, generateOther = writeOther)
                    else:
                        generateCtl(model=model, treefile = tree+".nl", seqfile=alignment, outfile=None, generateOther = writeOther)
                else:
                    generateCtl(model=model, treefile = t, seqfile=alignment, outfile=None, generateOther = writeOther)


    elif alignment:
        #pass#
        showAlignmentWithTree(tree, alignment)
    if tolabel:
        labelForPaml(unlabelledTree,tolabel, tree)
    else:
        showTreeNodes(unlabelledTree)

if __name__ == "__main__":
    main()

