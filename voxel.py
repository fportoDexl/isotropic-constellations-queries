from math import sqrt, pow


class Voxel(object):
    def __init__(self, x_left, x_right, y_left, y_right):
        self.x_left = x_left
        self.x_right = x_right
        self.y_left = y_left
        self.y_right = y_right
        self.centroid = self.computeCentroid()

    def getBottomLeft(self):
        return [self.x_left, self.y_left]

    def getUpRight(self):
        return [self.x_right, self.y_right]

    def setCentroidPos(self):
        self.centroid = self.computeCentroid()

    def getVoxelCentroid(self):
        return self.centroid

    def getHeightSize(self):
        return self.y_right - self.y_left

    def getSideSize(self):
        return self.x_right - self.x_left

    def isPointinVoxel(self, x, y):
        return self.x_left <= x <= self.x_right and self.y_left <= y <= self.y_right

    def getDiagonal(self):
        return sqrt(pow(self.x_right - self.x_left, 2) + pow(self.y_right - self.y_left, 2))

    def computeCentroid(self):
        x_left = self.getBottomLeft()
        x_l = x_left[0]
        y_l = x_left[1]
        x_right = self.getUpRight()
        x_r = x_right[0]
        y_r = x_right[1]
        cent_ra = ((x_r - x_l) / 2) + x_l
        cent_dec = ((y_r - y_l) / 2) + y_l
        centroid = [cent_ra, cent_dec]
        return centroid

    def __str__(self):
        return "%s,%s,%s,%s,%s" % (self.x_left, self.x_right, self.y_left, self.y_right, self.centroid)

    def __eq__(self, other):
        return self.x_left == other.x_left and self.x_right == other.x_right and self.y_left == other.y_left and \
               self.y_right == other.y_right
