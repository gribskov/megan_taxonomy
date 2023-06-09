"""=================================================================================================
megan_taxonomy:ncbi_taxonomy2megan.py

07 June 2023     gribskov
================================================================================================="""
import sys
from tree.tree import Tree


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
    :param node:
    :return:
    ---------------------------------------------------------------------------------------------"""
    nodes = open_safe(node, 'r')

    ranks = {}
    taxidx = {}
    # taxid==1 is the root, it has no parent
    taxidx['1'] = Tree('1')
    nodes.readline()

    node_n = 0
    for line in nodes:
        field = line.rstrip().replace('\t', '').split("|")
        taxid, parent, rank = field[:3]
        # print(f'taxon:{field[0]}\tparent:{field[1]}\trank:{field[2]}')
        if taxid not in taxidx:
            childnode = Tree(taxid)
            taxidx[taxid] = childnode
            node_n += 1
        else:
            childnode = taxidx[taxid]


        if parent not in taxidx:
            # parent taxon hasn't been created
            taxidx[parent] = Tree(parent)
            node_n += 1

        taxidx[parent].children.append(childnode)
        if rank in ranks:
            ranks[rank] += 1
        else:
            ranks[rank] = 1

    nnodes = Tree.nnodes
    a=taxidx['1232737']
    print(f'{nnodes} added to tree')
    n = tree_to_newick(taxidx['1'])

    return


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
        stackmax = max(len(stack), stackmax)
        node, rs = stack.pop()
        print(f'{count:8d}     {len(stack)}     {stackmax}     {node.name}')

        if node.children:
            ls += '('
            rs = f'){node.name}{rs}'
            for n in node.children[::-1]:
                stack.append([n, rs])
                rs = ','

        else:
            ls = f'{ls}{node.name}{rs}'

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
        'pecies subgroup': 100,
        'subspecies': 100,
        'varietas': 101,
        'forma': 101,
        'forma specialis': 101,
        'pathogroup': 101,
        'morph': 101,
        'biotype': 101,
        'genotype ': 101,
        'serogroup': 101,
        'clade': 101,
        'serotype': 101,
        'isolate': 101,
        'strain': 101}

    if rank in r2l:
        return r2l[rank]
    else:
        sys.stderr.write(f'rank_to_level: unknown rank "{rank}". lvl set to 0')
        return 0


# ==================================================================================================
# Main
# ==================================================================================================
if __name__ == '__main__':
    # tree = '((d,e,f)b,c,g)a;'
    # root = Tree(newick=tree)
    # print(tree_to_newick(root))

    tax2name = read_names('data/names.dmp.test')
    for taxid in tax2name:
        print(f'{taxid}\t{tax2name[taxid]}')

    tree = build_tree('data/nodes.dmp')

    exit(0)
