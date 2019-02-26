from util import EucDist


class QuadTree(object):
    def __init__(self, root, level):
        self.root = root
        self.level = level
        self.nodes = []  # nodes at the last level of the quadtree
        if level > 0:
            self.build_tree([root])
        else:
            self.nodes = [root]

    def getInitialLevel(self):
        return self.nodes

    def getRoot(self):
        return self.root

    def build_tree(self, nodes):
        for n in nodes:
            if n.getLevel() == self.level or n.getSize() < 3:
                self.nodes.append(n)
            else:
                self.build_tree(n.split_node())

    def find_neighbors(self, node_search, node, max_distance, neighbor_list):
        centroid = node_search.getCentroid()
        if node.visited is False:
            node_centroid = node.getCentroid()
            node_diagonal = node.voxel.getDiagonal() / 2
            dist = EucDist(centroid.Ra, node_centroid.Ra, centroid.Dec, node_centroid.Dec)

            if dist <= max_distance + node_diagonal:

                if len(node.children) > 0:
                    for child in node.children:
                        self.find_neighbors(node_search, child, max_distance, neighbor_list)
                else:
                    if node.getCentroid().getId() != centroid.getId():

                        # for element in node.elements:
                        #    if EucDist(element.Ra, centroid.Ra, element.Dec, centroid.Dec) <= max_distance:
                        #        neighbor_list.append(element)
                        neighbor_list.extend(node.elements)

                node.mark_visited()

        return neighbor_list
