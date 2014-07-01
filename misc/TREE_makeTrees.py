import getopt
import imp
import sys
import os.path
import re
import glob
from ete2 import Tree
from ete2 import EvolTree
#from ete2.treeview.layouts import evol_clean_layout
import subprocess
import dendropy
"""
read tree in newick format
"""
MAJORITY=0.5




def usage():
    print ("""
    #############################################
    # python TREE_makeTrees -t tree.nwk -a alignment.paml
    ############################################

    general options:
    -t, --tree=tree.nwk input tree file in newick format
     or
    -T, --bootstrapTree=TREE.raxmlbootstrap

    -a, --alignment=MSA.paml    alignment in paml format
    -n, --nodes=ListOfNodes     comma-separated list of nodes
                                that should be labelled
                                if no nodes are listed,
                                the tree will be traversed
                                and node ids are shown
    -h, --help  prints this

    """)
    sys.exit(2)



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

    #print(t)
    #t.run_model ('fb.example')
   # t.show()

    t.link_to_alignment(alignment, alg_format="paml")
    for node in t.traverse():
        print(node)
        print(node.node_id, node.name)
    #t.mark_tree([8], marks=["#1"])
    #print(t.write())
    print(alignment)
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
def main():
    tree = None
    tolabel = None
    alignment = None
    bstree = None
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "t:T:a:n:h", ["tree=","bottstrapTree=","alignment=","node=","help"])
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
        elif o in ("-T", "--bootstrapTree"):
            bstree = a
        
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
        
        with open(bstree+".nwk", 'w') as out:
            out.write(str(test))
        
        tree = bstree+".nwk"
        print(tree)
    
    
    if tree is None:
        usage()
    ########################################
    unlabelledTree = readTree(tree)
    #print(unlabelledTree)
    if alignment:
        #showAlignmentWithTree(unlabelledTree, alignment)
        showAlignmentWithTree(tree, alignment)
        
    if tolabel:
        labelForPaml(unlabelledTree,tolabel, tree)
    else:
        showTreeNodes(unlabelledTree)
if __name__ == "__main__":
    main()

