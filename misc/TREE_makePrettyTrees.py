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
from subprocess import call
import subprocess
#to handle nexus
import dendropy

"""
Produces control files for codeml with the given trees, alignment and model.
It will automatically produce ctl for contrasting hypotheses 
and different omega start values.

Example: -t tree1.nwk -a al1.paml -mM8 -n9 -A
will also create ctl files for M7 and M8a.
"""
Z=False
L=False
MAJORITY=0.5
PICS = []
PICDIR = ""
LOG = "./pic/Thumb_ln"
ZSC = "./pic/Thumb_Zs"
FILE = ""
INTERACTIVE = False
GENES = ["GRMZM2G039993","GRMZM2G116966",
         "GRMZM2G048371","GRMZM2G064302",
         "GRMZM2G156255","GRMZM2G310569",
         "GRMZM2G023946","GRMZM2G103847",
         "GRMZM2G155767","GRMZM5G833140",
         "GRMZM2G094452","GRMZM2G106578",
         "GRMZM2G003385","GRMZM5G833389",
         "GRMZM2G033962","GRMZM2G095727",
         "GRMZM2G007477","GRMZM2G059740",
         "GRMZM2G083016","GRMZM2G089136",
         "GRMZM2G150024","GRMZM2G384439",
         "GRMZM2G049076","GRMZM2G055331",
         "GRMZM2G050329","GRMZM2G166976",
         "GRMZM2G171354","GRMZM2G174990",
         "GRMZM2G005732","GRMZM2G095727",
         "GRMZM2G073571","GRMZM2G174990",
         "GRMZM2G005036","GRMZM2G443985",
         "GRMZM2G028535","GRMZM2G375504",
         "GRMZM2G062156","GRMZM2G064437",
         "AC233887.1_FG006","GRMZM2G175453",
         "GRMZM2G062429","GRMZM2G094452",
         "GRMZM2G089136","GRMZM2G382914",
         "GRMZM2G088196","GRMZM2G092793",
         "GRMZM2G077546","GRMZM2G098679",
         "GRMZM2G004932","GRMZM2G443985",
         "GRMZM2G055331","GRMZM5G875238",
         "GRMZM2G162529","GRMZM2G463280",
         "GRMZM2G010765","GRMZM2G166459",
         "GRMZM2G144372","GRMZM2G177912",
         "GRMZM2G080375","GRMZM2G443985",
         "GRMZM2G071119","GRMZM2G096365",
         "GRMZM2G443985","GRMZM5G879882",
         "GRMZM2G005996","GRMZM2G027891",
         "GRMZM2G092588","GRMZM2G175453",
         "GRMZM2G118737","GRMZM2G136139",
         "GRMZM2G077546","GRMZM2G160430",
         "GRMZM2G085019","GRMZM2G159724",
         "GRMZM2G055575","GRMZM2G174896",
         "GRMZM2G025854","GRMZM2G109383",
         "GRMZM2G083841","GRMZM2G473001",
         "GRMZM2G428027","GRMZM2G568636",
         "GRMZM2G062156","GRMZM2G096683",
         "GRMZM2G315848","GRMZM5G868679",
         "GRMZM2G045473","GRMZM2G326707",
         "GRMZM2G026024","GRMZM2G463280"]



def usage():
    print ("""
    #############################################
    # example usage:
    # python TREE_makePrettyTrees -t tree.nwk 
    # [-a alignment.paml] [-b] [-l] [-n "1,2,3"] [-L/-Z]
    ############################################

    general options:
    -t, --tree=tree.nwk input tree file in newick format
    or
    -T, --bstree=bootstraptrees.nwk bootstrapped trees will be converted (needs single file with bootstraps)
                                    will use 70percent majority
    [-a, --alignment=MSA.paml]  alignment in paml format
    
    [-n, --nodes="NUM,NUM,NUM"] comma-sep. list of nodes to be highlighted
  //[-o, --out=OUT]             output prefix (default: treefile specified with -t)
    [-l, --showBranchLengths]    
    [-b, --showBootstrapSupport]
    [-L, --log]                 pretty log pics
    [-Z, --zscore]              pretty zscore pics
    -h, --help                  prints this

    """)
    sys.exit(2)



def getUnlabelledTree(tree):
    t = Tree(tree)
    return(t.write(format=9))
def showTreeWithPictures(tree = None, alignment=None, branchLengths=True, bootstrapSupport=True, tolabel=None,showZScores=False,showLogs=False ):
    print(PICS)
    
    print("ShowTreeWithPictures",tree, alignment, branchLengths,bootstrapSupport, tolabel,showZScores,showLogs )
    if not alignment:
        nsFG = NodeStyle()
        nsFG["fgcolor"] = "darkgreen"
        nsFG["size"] = 8
        t = EvolTree(tree)
        #todo:label
        #
            
        for node in t.traverse():
            print(node.node_id)
            if tolabel:
                if str(node.node_id) in tolabel:
                     node.set_style(nsFG)
                #q'n'd 
            if (node.name.split("_")[0]+".png" in PICS):
                print(node.name.split("_")[0]+".png")
                node.add_face(ImgFace(PICDIR+os.sep+node.name.split("_")[0]+".png", height=50), column=1, position="aligned")
            #non GRZM identifier
            elif (node.name+".png" in PICS):
                print(node.name+".png")
                node.add_face(ImgFace(PICDIR+os.sep+node.name+".png", height=50), column=1, position="aligned")
            
            
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_length = branchLengths
        ts.show_branch_support = bootstrapSupport
        out = FILE
        if branchLengths:
            out+="_Len"
        if bootstrapSupport:
            out+="_Boot"
        if Z:
            out+="_Z"
        if L:
            out+="_L"
        t.render(out+"_tree.pdf",tree_style=ts)
        t.render(out+"_tree.png",tree_style=ts)
        if INTERACTIVE:
            t.show(tree_style=ts)
        
    else:
        t = EvolTree(tree, alignment,alg_format="paml")
        t.link_to_alignment(alignment,alg_format="paml")
        #todo label
        #todo check treestyle
        
        #ts = TreeStyle()
        #ts.show_leaf_name = True
        #ts.show_branch_length = branchLength
        #ts.show_branch_support = bootstrapSupport
        t.show()



def showTreeNodes(unlabelledTreeAsString):
    t = EvolTree(unlabelledTreeAsString)
    for node in t.traverse():
        #print(node)
        #print(node.node_id, node.name)
        #if (node.name.split("_")[0] in GENES):
        #    print(node.name, node.node_id)
        leaves  = node.get_leaf_names()
        leaves = [l for l in leaves if l.split("_")[0] in GENES ] 
        if leaves !=[]:
            print(node)
            print(leaves, node.name, node.node_id)
            print("\n")


def main():
    tree = None
    bstree = None
    tolabel = None
    alignment = None
    
    showBootstrapSupport = False
    showBranchLengths = False
    showZScores = False
    showLogs = False
    global PICS, PICDIR,FILE, INTERACTIVE,Z,L
    PICS = []
    FILE = ""
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "t:T:a:n:blZLi", ["tree=","bootstrapTree=""alignment=","node=","showBootstrapSupport","ShowBranchLengths","help","zscore","log","interactive"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-t", "--tree"):
            tree = a
            FILE = a
        elif o in ("-T", "--bootstrapTree"):
            bstree = a
            FILE =a
        elif o in ("-a", "--alignment"):
            alignment = a
        elif o in ("-n", "--nodes"):
            tolabel = a.split(",")
        elif o in ("-b", "--showBootstrapSupport"):
            showBootstrapSupport = True 
        elif o in ("-l", "--showBranchLengths"):
            showBranchLengths = True 
        elif o in ("-Z", "--zscore"):
            showZScores = True 
            Z=True
        elif o in ("-L", "--log"):
            showLogs = True 
            L=True
        elif o in ("-i", "--interactive"):
            INTERACTIVE = True
        else:
            assert False, "unhandled option"

    ########################################
    if bstree:
        res = ""
        nex = False
        p = subprocess.Popen('sumtrees.py --min-clade-freq='+str(MAJORITY)+" --percentages --decimals=0 "+bstree, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
             if line.startswith("#NEXUS") or nex:
                 nex = True
                 res+= line
                 retval = p.wait()
        with open(".tmptree0123.tmp","w") as tmp:
            tmp.write(res)
        test = dendropy.TreeList.get_from_path(".tmptree0123.tmp", "nexus")
        tree = str(test)
    if tree is None:
        usage()
    ########################################
    unlabelledTree = getUnlabelledTree(tree)
    print(FILE)
    print(unlabelledTree,tree, alignment)
    print(showBootstrapSupport,showBranchLengths )
    
    if showZScores:
        for file in os.listdir(ZSC):
            if file.endswith(".png"):
                PICS.append(file)
                PICDIR = ZSC
    elif showLogs:
        for file in os.listdir(LOG):
            if file.endswith(".png"):
                PICS.append(file)
                PICDIR = LOG
    if tree and alignment:
        showTreeWithPictures(tree=tree, alignment = alignment, bootstrapSupport = showBootstrapSupport, branchLengths = showBranchLengths, tolabel=tolabel,showZScores=showZScores,showLogs=showLogs)
    else:
        showTreeNodes(unlabelledTree)
        showTreeWithPictures(tree=tree, alignment = alignment, bootstrapSupport = showBootstrapSupport, branchLengths = showBranchLengths, tolabel=tolabel,showZScores=showZScores,showLogs=showLogs)
    
if __name__ == "__main__":
    main()

