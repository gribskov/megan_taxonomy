"""=================================================================================================
megan_taxonomy:ncbi_taxonomy2megan.py

07 June 2023     gribskov
================================================================================================="""
import sys
from tree.tree import Tree


class Tree_plus(Tree):
    """=============================================================================================
    Extended Tree class with additional attributes
        self.ls, self.rs        used by tree2newick to build the newick string
        self.depth              depth of a node in the tree (see also add depth
        self.rank               taxonomic rank, translated to megan level
    ============================================================================================="""

    def __init__(self, *args, **kwargs):
        super(Tree_plus, self).__init__(*args, **kwargs)

        # self.ls = ''
        self.rs = ','
        self.depth = 0
        self.rank = ''

    def createNode(self):
        """-----------------------------------------------------------------------------------------
        overrides Tree:createNode() so that childAdd() and newickLoad() functions create Tree_plus
        objects instead of Tree

        :return: Tree_plus
        -----------------------------------------------------------------------------------------"""
        return Tree_plus()

    def set_depth(self, depth_init=1):
        """-----------------------------------------------------------------------------------------
        traverse the tree and set depth of each node

        :param depth_init:int   depth of the root node
        :return: int            maximum depth
        -----------------------------------------------------------------------------------------"""
        self.depth = depth_init
        depth_max = depth_init
        n = 0
        for node in self:
            n += 1
            if node.children:
                d = node.depth + 1
                depth_max = max(depth_max, d)
                for child in node.children:
                    child.depth = d

        return depth_max


def open_safe(filename, mode):
    """---------------------------------------------------------------------------------------------
    open the file and catch exceptions. On failure, if mode is 'r' or 'rb' exit with status = 1,
    if mode is 'w', 'wb', or 'a', or 'ab', exit with status = 2.

    :param filename: string
    :param mode: string
    :return: filehandle
    ---------------------------------------------------------------------------------------------"""
    try:
        file = open(filename, mode)
    except OSError:
        sys.stderr.write(f'open_safe:Unable to open file "{filename}" in mode "{mode}"\n')
        status = 2
        if mode in 'rb':
            status = 1
        exit(status)

    return file


def read_names(name):
    """---------------------------------------------------------------------------------------------
    open and read the names.dmp file

    :param name: string     NCBI names.dmp file name
    :return: dict           key=taxid, value=name
    ---------------------------------------------------------------------------------------------"""
    names = open_safe(name, 'r')

    tax2name = {}
    for line in names:
        line = line.rstrip()
        field = line.replace('\t', '').split("|")
        if field[3] == 'scientific name':
            tax2name[field[0]] = field[1]

    names.close()
    return tax2name


def build_tree(node):
    """---------------------------------------------------------------------------------------------
    read the NCBI node.dmp and construct a tree using the Tree class

    :param node: str        name of the file containing the node data (e.g., nodes.dmp)
    :return: Tree           root node of tree
    ---------------------------------------------------------------------------------------------"""
    nodes = open_safe(node, 'r')
    tracelevel = 'order'

    ranks = {}
    taxidx = {}
    # taxid==1 is the root, it has no parent
    line = nodes.readline()
    field = line.rstrip().replace('\t', '').split("|")
    taxid, parent, rank = field[:3]
    root = Tree_plus(taxid, mode='dfs_stack')
    root.rank = rank_to_level(rank)
    taxidx[taxid] = root

    node_n = 0
    order_n = 0
    for line in nodes:
        field = line.rstrip().replace('\t', '').split("|")
        taxid, parent, rank = field[:3]
        if rank == tracelevel:
            sys.stderr.write('.')
            order_n += 1
            if not order_n % 50:
                sys.stderr.write(f'\t{order_n}\n')

        # print(f'taxon:{field[0]}\tparent:{field[1]}\trank:{field[2]}')
        if taxid not in taxidx:
            childnode = Tree_plus(taxid)
            taxidx[taxid] = childnode
            node_n += 1
        else:
            childnode = taxidx[taxid]

        childnode.rank = rank_to_level(rank)

        if parent not in taxidx:
            # parent taxon hasn't been created
            taxidx[parent] = Tree_plus(parent)
            node_n += 1

        taxidx[parent].children.append(childnode)
        if rank in ranks:
            ranks[rank] += 1
        else:
            ranks[rank] = 1

    sys.stderr.write(f'\n')
    return root


def tree_to_newick(node):
    """---------------------------------------------------------------------------------------------
    stack-based construction of newick string from a root node
    :param node: Tree       root node of tree
    :return: string         newick string
    ---------------------------------------------------------------------------------------------"""
    stack = []
    ls = ''
    rs = ';'
    stack.append([node, rs])
    count = 0
    stackmax = 0
    while stack:
        count += 1
        if not count % 2000:
            sys.stderr.write(':')
        if not count % 100000:
            sys.stderr.write(f'\t{count}\n')

        stackmax = max(len(stack), stackmax)
        node, rs = stack.pop()
        # print(f'{count:8d}     {len(stack)}     {stackmax}     {node.name}')

        if node.children:
            ls += '('
            rs = f'){node.name}{rs}'
            for n in node.children[::-1]:
                stack.append([n, rs])
                rs = ','

        else:
            ls = f'{ls}{node.name}{rs}'

    return ls


def tree_to_newick2(root, depth_max=50):
    """---------------------------------------------------------------------------------------------
    stack-based construction of newick string from a root node
    :param root: Tree       root node of tree
    :return: string         newick string
    ---------------------------------------------------------------------------------------------"""
    ls = ''

    count = 0
    for node in root:
        # count += 1
        # if not count % 2000:
        #     sys.stderr.write(':')
        # if not count % 100000:
        #     sys.stderr.write(f'\t{count}\n')

        if node.depth > depth_max:
            continue

        if node.children and node.depth < depth_max:
            ls += '('
            # node.children[0].ls = ''
            node.children[-1].rs = f'){node.name}{node.rs}'

        else:
            ls += f'{node.name}{node.rs}'

    # replace final comma with ;
    ls = f'{ls[:-1]};'
    return ls


def rank_to_level(rank):
    """---------------------------------------------------------------------------------------------
    Return the megan level (lvl) for a string representing a taxonomic rank
    
    :param rank: 
    :return: 
    ---------------------------------------------------------------------------------------------"""
    r2l = {
        'no rank': 0,
        'superkingdom': 1,
        'kingdom': 1,
        'subkingdom': 1,
        'superphylum': 2,
        'phylum': 2,
        'subphylum': 2,
        'superclass': 3,
        'class': 3,
        'subclass': 3,
        'infraclass': 3,
        'cohort': 3,
        'subcohort': 3,
        'superorder': 4,
        'order': 4,
        'suborder': 4,
        'infraorder': 4,
        'parvorder': 4,
        'superfamily': 5,
        'family': 5,
        'subfamily': 5,
        'tribe': 5,
        'subtribe': 5,
        'genus': 98,
        'subgenus': 98,
        'series': 98,
        'section': 98,
        'subsection': 98,
        'species group': 99,
        'species': 99,
        'species subgroup': 100,
        'subspecies': 100,
        'varietas': 101,
        'forma': 101,
        'forma specialis': 101,
        'pathogroup': 101,
        'morph': 101,
        'biotype': 101,
        'genotype': 101,
        'serogroup': 101,
        'clade': 101,
        'serotype': 101,
        'isolate': 101,
        'strain': 101}

    for r in r2l:
        # make the values strings
        r2l[r] = f'{r2l[r]}'

    if rank in r2l:
        return r2l[rank]
    else:
        sys.stderr.write(f'\nrank_to_level: unknown rank "{rank}". lvl set to 0\n')
        return 0


def write_map_file(mapfile, tax2name, root, depth_max=50):
    """---------------------------------------------------------------------------------------------
    write out the mapping between numeric taxid, scientific name, and taxonomic rank, tab-delimited

    1       NCBI    -1      0
    2       Bacteria        -1      0
    6       Azorhizobium    -1      98
    7       Azorhizobium caulinodans        -1      100
    9       Buchnera aphidicola     -1      100
    10      Cellvibrio      -1      98
    11      Cellulomonas gilvus     -1      100

    :param mapfile: filehandle      open for writing
    :param tax2name: dict           key=taxid, value =scientific name
    :param root: Tree object        root of tree
    :return: int                    number of taxa writen
    ---------------------------------------------------------------------------------------------"""
    tax_n = 0
    for node in root:
        if node.depth > depth_max:
            continue

        tax_n += 1
        try:
            taxon = tax2name[node.name]
        except:
            sys.stderr.write(f'write_map_file - '
                             f'unknown taxon "{node.name}" assigned as "unknown""')
            taxon = 'unknown'

        mapfile.write(f'{node.name}\t{taxon}\t-1\t{node.rank}\n')

    return tax_n


# ==================================================================================================
# Main
# ==================================================================================================
if __name__ == '__main__':
    # # tree = '((((l,m,n)h,i)d,(aa,bb)e,(cc)f)b,(j,k)c,g)a;'
    # tree = '((6,7,8,9)10239,(2,3,2759)131567,(12908)2787823,(28385)2787854)1;'
    # root = Tree_plus(newick=tree, mode='dfs_stack')
    # max1 = root.set_depth(0)
    # max0 = root.set_depth(1)
    # newick = tree_to_newick2(root)
    # print(newick)
    # newick = tree_to_newick2(root, 2)
    # print(newick)
    # max1 = root.set_depth(1)
    # max0 = root.set_depth(0)
    # print(f'max depth 0={max0}\t1={max1}')
    # exit(100)

    nodefilename = 'data/nodes.dmp'
    namefilename = 'data/names.dmp'
    newickfilename = 'new.tre'
    mapfilename = 'new.map'
    depth = 4

    tax2name = read_names(namefilename)
    # for taxid in tax2name:
    #     print(f'{taxid}\t{tax2name[taxid]}')

    root = build_tree(nodefilename)
    maxtreedepth = root.set_depth(1)
    nnodes = Tree.nnodes
    sys.stderr.write(f'\n{nnodes} taxa added to tree. maximum tree depth={maxtreedepth}\n')
    sys.stderr.write('transforming tree to Newick format\n')
    newick = tree_to_newick2(root, depth)
    sys.stderr.write(f'{newick}\n')
    sys.stderr.write(f'{newick[:100]}\n...\n{newick[-100:]}\n')
    newickfilename = 'new.tre'
    sys.stderr.write(f'writing tree to {newickfilename} in Newick format\n')
    trefile = open_safe(newickfilename, 'w')
    trefile.write(newick)
    trefile.write('\n')
    trefile.close()

    mapfile = open_safe(mapfilename, 'w')
    map_n = write_map_file(mapfile, tax2name, root, depth)
    sys.stderr.write(f'{map_n} mappings written to {mapfilename}')
    mapfile.close()

    exit(0)
