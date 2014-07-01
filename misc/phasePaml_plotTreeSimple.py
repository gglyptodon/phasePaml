#########################
# last update:
# Wed 09 Apr 2014 09:00:03 AM CEST
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
from subprocess import call
import subprocess
#to handle nexus
import dendropy

"""
Plots a tree with its alignment.
Example: phasePaml_plotTreeSimple -t tree1.nwk -a al1.paml

"""
INTERACTIVE = False
MAJORITY=0.7
FILE = None
def usage():
    print ("""
    #############################################
    # example usage:
    # python phasePaml_plotTreeSimple -t tree.nwk
    # [-a alignment.paml] [-b] [-l] [-n "1,2,3"]
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
    [-i, --interactive]         show tree in viewer
    -h, --help                  prints this

    """)
    sys.exit(2)



def getUnlabelledTree(tree):
    t = Tree(tree)
    return(t.write(format=9))
def showTreeWithPictures(tree = None, alignment=None, branchLengths=True, bootstrapSupport=True, tolabel=None):

    print("ShowTreeWithPictures",tree, alignment, branchLengths,bootstrapSupport, tolabel)
    if alignment:
        t = EvolTree(tree, alignment,alg_format="paml")
        t.link_to_alignment(alignment,alg_format="paml")


    else:
        t = EvolTree(tree)

    nsFG = NodeStyle()
    nsFG["fgcolor"] = "darkgreen"
    nsFG["size"] = 8

    for node in t.traverse():
        print(node.node_id)
        if tolabel:
            if str(node.node_id) in tolabel:
                 node.set_style(nsFG)

    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = branchLengths
    ts.show_branch_support = bootstrapSupport
    out = FILE
    if branchLengths:
        out+="_Len"
    if bootstrapSupport:
        out+="_Boot"
    t.render(out+"_tree.pdf",tree_style=ts)
    t.render(out+"_tree.png",tree_style=ts)
    if INTERACTIVE:
        t.show(tree_style=ts)



def main():
    tree = None
    bstree = None
    tolabel = None
    alignment = None

    showBootstrapSupport = False
    showBranchLengths = False
    showZScores = False
    showLogs = False
    global INTERACTIVE
    global FILE
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "t:T:a:n:bli", ["tree=","bootstrapTree=""alignment=","node=","showBootstrapSupport","ShowBranchLengths","help","interactive"])
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

    if tree :
        showTreeWithPictures(tree=tree, alignment = alignment, bootstrapSupport = showBootstrapSupport, branchLengths = showBranchLengths, tolabel=tolabel)

if __name__ == "__main__":
    main()

