import networkx
from networkx.algorithms import isomorphism
import gemmi

""" Creates a networkx graph from a gemmi object
    Taken from: https://gemmi.readthedocs.io/en/latest/analysis.html#graph-isomorphism
    Input: cc entry from gemmi
    Output: networkx graph
"""
def graph_from_chemcomp(cc):
    G = networkx.Graph()
    for atom in cc.atoms:
        G.add_node(atom.id, Z=atom.el.atomic_number)
    for bond in cc.rt.bonds:
        G.add_edge(bond.id1.atom, bond.id2.atom)
    return G

""" Converts a cif file to a networkx graph
    Input: filepath (fp) of a cif file
    Output: networkx graph of a crystal structure
"""
def cif_to_graph(fp):
    doc = gemmi.cif.read_file(fp)
    print(doc)
    block = doc.sole_block()  # mmCIF has exactly one block
    print(block.name)

    for item in block:
        print(item.line_number)
        if item.pair is not None:
            print('pair', item.pair)
        elif item.loop is not None:
            print('loop', item.loop)
        elif item.frame is not None:
            print('frame', item.frame)
    #Current problem - no atoms seem to be added
    print(list(block.find_loop('_atom_site_label')))

    cc1 = gemmi.make_chemcomp_from_block(block)
    cc1.remove_hydrogens() #Not sure why this is included - take this out later?
    print(cc1)
    G = graph_from_chemcomp(cc1)
    print(G.nodes())

def main():
    #Convert cif to networkx graph
    cc1 = cif_to_graph("Minerals/Data/CIF_Files/9014258.cif")

    #Perform MCS algorithm

    #Go from there...

if __name__ == "__main__":
    main()
