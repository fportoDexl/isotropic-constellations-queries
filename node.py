from collections import deque
from math import sqrt

from util import EucDist
from voxel import Voxel


class Node(object):
    # Structure:
    # voxel: The definition of the spatial region covered by this Node
    # centroid: the element (star) closest to the spatial centroid of the Node
    # parent: the parent node in the tree
    # level: the node level in the tree

    def __init__(self, voxel, centroid=None, parent=None, level=0):
        self.centroidStar = centroid
        self.voxel = voxel
        self.elements = deque()
        self.children = []
        self.size = 0
        self.parent = parent
        self.level = level
        self.visited = False

    def getElements(self):
        return self.elements

    def getChildren(self):
        return self.children

    def getVoxel(self):
        return self.voxel

    def getCentroid(self):
        return self.centroidStar

    def getLevel(self):
        return self.level

    def getSize(self):
        return self.size

    def addCentroid(self, centroid):
        self.centroidStar = centroid

    def addChildren(self, node):
        self.children.append(node)

    def addParent(self, node):
        self.parent = node

    def addElement(self, star):
        self.elements.append(star)

    def uniqueElement(self):
        unique = False
        if self.size == 1:
            unique = True
        return unique

    def setSize(self, size):
        self.size = size

    def isEmpty(self):
        if not self.elements:
            return True

    def mark_visited(self):
        if len(self.children) > 0:
            size = len(self.children)
            count = 0

            for n in self.children:
                if n.visited is True:
                    count += 1

            if count == size:
                self.visited = True
        else:
            pass

    def split_node(self):
        level = self.level + 1

        voxel = self.getVoxel()
        centroid_voxel = voxel.getVoxelCentroid()

        voxel11 = Voxel(voxel.x_left, centroid_voxel[0], centroid_voxel[1], voxel.y_right)
        quarter11 = Node(voxel11, None, self, level)

        voxel12 = Voxel(centroid_voxel[0], voxel.x_right, centroid_voxel[1], voxel.y_right)
        quarter12 = Node(voxel12, None, self, level)

        voxel01 = Voxel(voxel.x_left, centroid_voxel[0], voxel.y_left, centroid_voxel[1])
        quarter01 = Node(voxel01, None, self, level)

        voxel02 = Voxel(centroid_voxel[0], voxel.x_right, voxel.y_left, centroid_voxel[1])
        quarter02 = Node(voxel02, None, self, level)

        star11_cent = star12_cent = star01_cent = star02_cent = float("inf")
        star11_med_x = star11_med_y = star12_med_x = star12_med_y = star01_med_x = star01_med_y = star02_med_x = \
            star02_med_y = 0
        tot_star11 = tot_star12 = tot_star01 = tot_star02 = 0
        border_tot_star11 = border_tot_star12 = border_tot_star01 = border_tot_star02 = 0

        for _ in range(len(self.getElements())):
            star = self.elements.popleft()
            if star.getRa() < centroid_voxel[0]:
                if star.getDec() < centroid_voxel[1]:
                    quarter01.addElement(star)
                    if star.border is False:
                        star01_med_x += star.getRa()
                        star01_med_y += star.getDec()
                        tot_star01 += 1
                    else:
                        border_tot_star01 += 1
                else:
                    quarter11.addElement(star)
                    if star.border is False:
                        star11_med_x += star.getRa()
                        star11_med_y += star.getDec()
                        tot_star11 += 1
                    else:
                        border_tot_star11 += 1
            else:
                if star.getDec() > centroid_voxel[1]:
                    quarter12.addElement(star)
                    if star.border is False:
                        star12_med_x += star.getRa()
                        star12_med_y += star.getDec()
                        tot_star12 += 1
                    else:
                        border_tot_star12 += 1
                else:
                    quarter02.addElement(star)
                    if star.border is False:
                        star02_med_x += star.getRa()
                        star02_med_y += star.getDec()
                        tot_star02 += 1
                    else:
                        border_tot_star02 += 1

        node11_cent = node12_cent = node01_cent = node02_cent = None
        nodes = []

        if tot_star11 > 0:
            star11_med_x = star11_med_x / tot_star11
            star11_med_y = star11_med_y / tot_star11
            for star in quarter11.getElements():
                dist = EucDist(star.getRa(), star11_med_x, star.getDec(), star11_med_y)
                if star.border is False and dist < star11_cent:
                    node11_cent = star
                    star11_cent = dist
            quarter11.centroidStar = node11_cent
            size = tot_star11 + border_tot_star11
            quarter11.setSize(size)
            self.addChildren(quarter11)
            nodes.append(quarter11)

        if tot_star12 > 0:
            star12_med_x = star12_med_x / tot_star12
            star12_med_y = star12_med_y / tot_star12
            for star in quarter12.getElements():
                dist = EucDist(star.getRa(), star12_med_x, star.getDec(), star12_med_y)
                if star.border is False and dist < star12_cent:
                    node12_cent = star
                    star12_cent = dist
            quarter12.centroidStar = node12_cent
            size = tot_star12 + border_tot_star12
            quarter12.setSize(size)
            self.addChildren(quarter12)
            nodes.append(quarter12)

        if tot_star01 > 0:
            star01_med_x = star01_med_x / tot_star01
            star01_med_y = star01_med_y / tot_star01
            for star in quarter01.getElements():
                dist = EucDist(star.getRa(), star01_med_x, star.getDec(), star01_med_y)
                if star.border is False and dist < star01_cent:
                    node01_cent = star
                    star01_cent = dist
            quarter01.centroidStar = node01_cent
            size = tot_star01 + border_tot_star01
            quarter01.setSize(size)
            self.addChildren(quarter01)
            nodes.append(quarter01)

        if tot_star02 > 0:
            star02_med_x = star02_med_x / tot_star02
            star02_med_y = star02_med_y / tot_star02
            for star in quarter02.getElements():
                dist = EucDist(star.getRa(), star02_med_x, star.getDec(), star02_med_y)
                if star.border is False and dist < star02_cent:
                    node02_cent = star
                    star02_cent = dist
            quarter02.centroidStar = node02_cent
            size = tot_star02 + border_tot_star02
            quarter02.setSize(size)
            self.addChildren(quarter02)
            nodes.append(quarter02)

        return nodes

    def partial_match(self, node_j, distance, epsilon):
        voxel_i = self.getVoxel()
        size_i = voxel_i.getSideSize()

        centroid_i = self.getCentroid()
        centroid_j = node_j.getCentroid()
        dist = EucDist(centroid_i.getRa(), centroid_j.Ra, centroid_i.getDec(), centroid_j.getDec())

        if size_i < distance:
            if distance + epsilon + (size_i * sqrt(2) >= dist >= distance - epsilon - (size_i * sqrt(2))):
                match = True
            else:
                match = False
        else:
            if (distance + epsilon) >= dist >= (distance - epsilon):
                match = True
            else:
                match = False
        return match

    def __str__(self):
        return "%s;%s;%s;%s;%s;%s" % (self.centroidStar, self.voxel, len(self.elements), len(self.children),
                                      self.size, self.level)
