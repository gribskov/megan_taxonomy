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
    read the NCBI node.dmp and construct a tree using the tree package
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

        if parent not in taxidx:
            # parent taxon hasn't been created
            taxidx[parent] = Tree(parent)
            node_n += 1

        taxidx[parent].children.append(childnode)
        if rank in ranks:
            ranks[rank] += 1
        else:
            ranks[rank] = 1

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
    while stack:
        node, rs = stack.pop()
        if node.children:
            ls += '('
            rs = f'){node.name}{rs}'
            for n in node.children[::-1]:
                stack.append([n, rs])
                rs = ','

        else:
            ls = f'{ls}{node.name}{rs}'

    return ls


# ==================================================================================================
# Main
# ==================================================================================================
if __name__ == '__main__':
    tree = '((d,e,f)b,c,g)a;'
    root = Tree(newick=tree)
    print(tree_to_newick(root))

    # tax2name = read_names('data/names.dmp.test')
    # for taxid in tax2name:
    #     print(f'{taxid}\t{tax2name[taxid]}')
    #
    # tree = build_tree('data/nodes.dmp')

    exit(0)
