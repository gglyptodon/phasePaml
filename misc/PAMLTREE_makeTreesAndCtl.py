import getopt
import imp
import sys
import os.path
import re
import glob
from ete2 import Tree
from ete2 import EvolTree

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
    # python PAMLTREE_makeTreesAndCtl -t tree.nwk -a alignment.paml
    ############################################

    general options:
    -t, --tree=tree.nwk input tree file in newick format
    [-a, --alignment=MSA.paml]    alignment in paml format
    [-n, --nodes=ListOfNodes ]    comma-separated list of nodes
                                that should be labelled
                                if no nodes are listed,
                                the tree will be traversed
                                and node ids are shown
    [-m, --model=MODEL]         MODEL can be one of ["M0","M1a","M2a",
                                "M7","M8","M8a", "Ah0", "Ah1"]
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
        ctl0 = ctlTemplate(treefile = treefile, seqfile = seqfile, outfile = outfile+outsuff, model= "1", NSsites = "0",fix_omega = 0, omega = 1)
        #free ratio model  - use is discouraged
        ################# write control files ###############
        with open(outfile+outsuff0+".ctl", 'w') as out:
            out.write(ctl0)
        
    elif (model.upper() =="TR"):
        outsuff0 = ".TR_f1"
        outsuff1 = ".TR_e1"
        #NSsites ignored?
        ctl0 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile = outfile+outsuff, model= "2", NSsites = "0",fix_omega = 1, omega = 1)
        ctl1 = ctlTemplate(treefile = treefile, seqfile = seqfile,outfile = outfile+outsuff, model= "2", NSsites = "0",fix_omega = 0, omega = 1)
        ################# write control files ###############
        with open(outfile+outsuff0+".ctl", 'w') as out:
            out.write(ctl0)
        with open(outfile+outsuff0+".ctl", 'w') as out:
            out.write(ctl0)
    
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
        #null
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

def showAlignmentWithTree(unlabelledTreeAsString,alignment):
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

def showTreeNodes(unlabelledTreeAsString):
    t = EvolTree(unlabelledTreeAsString)
    for node in t.traverse():
        print(node)
        print(node.node_id, node.name)

def modelhelp():
    print ("""
    ############################################
    # MODELS  
    ############################################
    modelID:    Description:            vs:
      M0          one-ratio model;       ?
    [ M1 ]    not implemented
      M1a    branch nearly neutral      M2a
    [ M2 ]    not implemented
      M2a    branch pos. sel.           M1a
      M7     branch beta                M8
      M8     branch beta+omega>1,alt    M7, M8a
      M8a    branch beta+omega=1,null   M8
      Ah0    branch-site                Ah1
      Ah1    branch-site                Ah0
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
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "t:a:n:m:hAM", ["tree=","alignment=","node=","model=","help", "writeOther","modelhelp"])
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
        elif o in ("-m", "--model"):
            model = a
        elif o in ("-A", "--writeOther"):
            writeOther = True
        elif o in ("-M", "--modelhelpr"):
            modelhelp()
        else:
            assert False, "unhandled option"

    ########################################
    if tree is None:
        usage()
    ########################################
    unlabelledTree = readTree(tree)
    #print(unlabelledTree)
    print(model, tree, alignment)
    if model and tree and alignment:
        generateCtl(model=model, treefile = tree, seqfile=alignment, outfile=None, generateOther = writeOther)
    elif alignment:
        showAlignmentWithTree(unlabelledTree, alignment)
    if tolabel:
        labelForPaml(unlabelledTree,tolabel, tree)
    else:
        showTreeNodes(unlabelledTree)
    
if __name__ == "__main__":
    main()

