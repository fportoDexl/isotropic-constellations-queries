class Element(object):
    def __init__(self, pointId, ra, dec, u, u_err, g, g_err, r, r_err, i, i_err, z, z_err, distToCentroid, border):
        self.pointId = pointId
        self.Ra = ra
        self.Dec = dec
        self.u = u
        self.u_err = u_err
        self.g = g
        self.g_err = g_err
        self.r = r
        self.r_err = r_err
        self.i = i
        self.i_err = i_err
        self.z = z
        self.z_err = z_err
        self.distToCentroid = distToCentroid
        self.border = border

    def getId(self):
        return self.pointId

    def getRa(self):
        return self.Ra

    def getDec(self):
        return self.Dec

    def getDistanceTOC(self):
        return self.distToCentroid

    def __str__(self):
        return "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % (self.pointId, self.Ra, self.Dec, self.u, self.u_err,
                                                           self.g, self.g_err, self.r, self.r_err, self.i,
                                                           self.i_err, self.z, self.z_err)

    def __eq__(self, other):
        if type(other) is Element:
            return self.pointId == other.pointId
