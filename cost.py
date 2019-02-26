from util import EucDist


def produce_star_pairs(node_i, node_j, node_k, distance, dq1, dq2, epsilon, eid, scale, relations):
    """

    :param node_i:
    :param node_j:
    :param node_k:
    :param distance:
    :param dq1:
    :param dq2:
    :param epsilon:
    :param eid:
    :param scale:
    :param relations:
    :return:
    """
    elements_node_i = node_i.getElements()
    elements_node_j = node_j.getElements()
    elements_node_k = node_k.getElements()

    for element_i in elements_node_i:
        ra_i = element_i.Ra
        dec_i = element_i.Dec

        for element_j in elements_node_j:
            ra_j = element_j.Ra
            dec_j = element_j.Dec
            calculate_distance_i_j = EucDist(ra_i, ra_j, dec_i, dec_j)

            if (distance - epsilon) <= calculate_distance_i_j <= (distance + epsilon):
                for element_k in elements_node_k:
                    ra_k = element_k.Ra
                    dec_k = element_k.Dec
                    calculate_distance_i_k = EucDist(ra_i, ra_k, dec_i, dec_k)

                    if (dq1 - epsilon) <= calculate_distance_i_k <= (dq1 + epsilon):
                        calculate_distance_j_k = EucDist(ra_j, ra_k, dec_j, dec_k)

                        if (dq2 - epsilon) <= calculate_distance_j_k <= (dq2 + epsilon):
                            if relations.get((element_i.pointId, element_j.pointId)):
                                relations[element_i.pointId, element_j.pointId].append((element_i, element_j,
                                                                                        element_k, eid, scale))
                            else:
                                relations[element_i.pointId, element_j.pointId] = []
                                relations[element_i.pointId, element_j.pointId].append((element_i, element_j,
                                                                                        element_k, eid, scale))
