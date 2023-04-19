import sys
from ete3 import Tree


def root_at_midpoint(tree_in, tree_in_rooted):

    t = Tree(tree_in, quoted_node_names=True)
    midpoint = t.get_midpoint_outgroup()
    t.set_outgroup(midpoint)
    t.write(tree_in_rooted)



