from util import EucDist


class Tuple(object):
    def __init__(self):
        self.metadata = []
        self.stars = []
        self.scale = None

    def addTuple(self, tuple):
        self.stars.extend(tuple.getStars())
        self.metadata.extend(tuple.getMetadata())
        self.scale = tuple.getScale()

    def addStar(self, s, rel):
        self.stars.append(s)
        self.metadata.append(rel)

    def getMetadataPos(self, index):
        return self.metadata[index]

    def getStar(self, index):
        return self.stars[index]

    def getScale(self):
        return self.scale

    def getStars(self):
        return self.stars

    def getMetadata(self):
        return self.metadata


class Plan(object):
    def __init__(self, query=None):
        self.query = query

    def build_plan(self, relations, filter_operation=False):
        join_old = None

        for i in range(len(relations.getRelations()) - 1):
            if i == 0:
                join = Join(relations.getRelations()[i + 1], None, relations.getRelations()[i], self.query)
            else:
                join = Join(relations.getRelations()[i + 1], join_old, None, self.query)
            join_old = join
        if filter_operation:
            root = Filter(join_old, self.query)
        else:
            root = join_old

        return root


class Filter(object):
    def __init__(self, producer=None, query=None):
        self.producer = producer
        # self.type="filter"
        self.query = query
        self.dist = query.distances  # distG
        # self.anchorElement=self.query.getAnchor().getId()
        # self.margin=query.getMargin()
        # self.total = 0
        # self.scale_total = 0
        # self.true_tuples = 0
        # self.scale_true_tuples = 0

    def getNext(self):

        listTuple = self.producer.getNext()

        if listTuple == "end":
            return listTuple
        elif not listTuple:
            return listTuple
        else:
            listResult = []
            # self.total += len(listTuple)
            for t in listTuple:

                # if t.getScale() > 1:
                #     self.scale_total += 1
                result = self.checkScale(t, self.query.getEpsilon())
                if result:
                    # if t.getScale() > 1:
                    #     self.scale_true_tuples += 1
                    # else:
                    #     self.true_tuples += 1
                    listResult.append(t)

        if len(listResult) == 0:
            return "end"
        else:
            return listResult

    def checkScale(self, candidatesolution, epsilon):

        stars_tuple = candidatesolution.getStars()
        metadata_tuple = candidatesolution.getMetadata()
        MINMAX = float("inf")
        MAXMIN = float("-inf")
        calculated_epsilon = (candidatesolution.getScale() * epsilon)
        # print("\n\n\nENDDDDDDDDDD")

        # print("Scale = %s" % candidatesolution.getScale())
        for i in range(len(stars_tuple)):
            # id_i = stars_tuple[i].getId()
            ra_i = stars_tuple[i].getRa()
            dec_i = stars_tuple[i].getDec()
            rel_i = metadata_tuple[i]
            # print("\n(i) ID = %s\tRA = %s\tDEC = %s\tEID = %s" % (id_i,ra_i, dec_i, rel_i))
            for j in range(i + 1, len(stars_tuple)):
                # id_j = stars_tuple[j].getId()
                ra_j = stars_tuple[j].getRa()
                dec_j = stars_tuple[j].getDec()
                rel_j = metadata_tuple[j]
                dist = EucDist(ra_i, ra_j, dec_i, dec_j)
                # print("(j) ID = %s\tRA = %s\tDEC = %s\tEID = %s\t Dist = %s" % (id_j, ra_j, dec_j, rel_j, dist))

                scale_factor = dist / self.dist[rel_i][rel_j]
                minn = scale_factor - calculated_epsilon
                maxx = scale_factor + calculated_epsilon

                if minn > MAXMIN:
                    MAXMIN = minn
                if maxx < MINMAX:
                    MINMAX = maxx
        if MAXMIN > MINMAX:
            return False

        return True


class Join(object):
    def __init__(self, rightRelation, leftJoin=None, leftRelation=None, query=None, theta=0.1):
        self.rightR = rightRelation
        self.idR = self.rightR.getRelationId()
        self.type = "join"
        self.leftJ = leftJoin
        self.leftRelation = leftRelation
        self.query = query
        self.distM = query.getDistanceMatrix()
        self.scale_index = 0
        self.theta = theta
        if self.leftJ is None:
            self.stars = self.leftRelation.getStars()
            self.starIndex = 0
            self.sizeLeftleave = len(self.stars)

    def getNext(self):
        listTuple = []
        newListTuple = []
        if self.leftRelation is None:
            listTuple = self.leftJ.getNext()
            if listTuple is False:
                return listTuple
            elif listTuple == "end":
                return listTuple
        else:
            if self.starIndex == self.sizeLeftleave:
                listTuple = "end"
                return listTuple
            else:
                listTuple.append(self.stars[self.starIndex])
                self.starIndex += 1
        # match = True

        # scale = 1
        for l in range(len(listTuple)):
            genTuple = listTuple[l]
            tuple = Tuple()
            if self.leftJ is None:
                scale = genTuple[3]
                tuple.addStar(genTuple[0], self.query.maxDistGeneralPair[0])
                tuple.addStar(genTuple[1], self.query.maxDistGeneralPair[1])
                tuple.addStar(genTuple[2], self.leftRelation.getRelationId())
                tuple.scale = scale
            else:
                scale = genTuple.getScale()
                tuple.addTuple(genTuple)

            distM = self.query.getDistanceMatrix()

            dist_sm_plus = [[0 for j in range(len(self.query.getQuery()) + 1)] for i in
                            range(len(self.query.getQuery()) + 1)]

            dist_sm_minus = [[0 for j in range(len(self.query.getQuery()) + 1)] for i in
                             range(len(self.query.getQuery()) + 1)]

            calculated_epsilon = scale * self.query.getEpsilon()

            for i in range(len(self.query.getQuery()) + 1):
                for j in range(len(self.query.getQuery()) + 1):
                    dist_sm_plus[i][j] = distM[i][j] * (scale + calculated_epsilon)
                    dist_sm_minus[i][j] = distM[i][j] * (scale - calculated_epsilon)

            for t in self.rightR.getStars():
                match = True
                for i in range(len(tuple.getStars()) - 2):  # Checks pairwise distances with all elements in tuple
                    dist = EucDist(tuple.getStar(i + 2).getRa(), t[2].getRa(), tuple.getStar(i + 2).getDec(),
                                   t[2].getDec())
                    distanceMatrix_plus = dist_sm_plus[self.rightR.getRelationId()][tuple.getMetadataPos(i + 2)]
                    distanceMatrix_minus = dist_sm_minus[self.rightR.getRelationId()][tuple.getMetadataPos(i + 2)]

                    if ((distanceMatrix_minus - (2 * calculated_epsilon) <= dist) and (
                                dist <= distanceMatrix_plus + (2 * calculated_epsilon))):
                        anchor_1 = tuple.getStar(0)
                        anchor_2 = tuple.getStar(1)
                        s = tuple.getStar(i + 2)
                        if s.getId() == t[2].getId() or anchor_1.getId() == t[2].getId() or \
                                        anchor_2.getId() == t[2].getId():
                            match = False
                            break
                    else:
                        match = False
                        break

                if (match):
                    newTuple = Tuple()
                    newTuple.addTuple(tuple)
                    newTuple.addStar(t[2], self.rightR.getRelationId())
                    newListTuple.append(newTuple)

        return newListTuple


class Relation(object):
    """

    """

    def __init__(self, query_id=0, stars=[]):
        """

        :param query_id:
        :param stars:
        """
        self.query_id = query_id
        self.stars = stars

    def addStar(self, anchor_1, anchor_2, partners, scale):
        """

        :param anchor_1:
        :param anchor_2:
        :param partners:
        :param scale:
        :return:
        """
        self.stars.append([anchor_1, anchor_2, partners, scale])

    def getStars(self):
        """

        :return:
        """
        return self.stars

    def getRelationId(self):
        """

        :return:
        """
        return self.query_id


class Relations(object):
    """

    """

    def __init__(self):
        """

        """
        self.relations = []

    def addRelation(self, relation):
        """

        :param relation:
        :return:
        """
        self.relations.append(relation)

    def getRelations(self):
        """

        :return:
        """
        return self.relations
