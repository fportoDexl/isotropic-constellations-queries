from collections import defaultdict
from os import path
from sys import argv
from time import time

from pyspark import SparkContext, SparkConf

# from cost import produce_star_pairs
from element import Element
from executionengine import Relation, Relations, Plan
from node import Node
from quadtree import QuadTree
from query import Query
from util import EucDist, compute_level
from voxel import Voxel


def build_tree(text):
    query = query_broadcast.value

    root = None
    geometric_centroid_ra = geometric_centroid_dec = None
    centroid = None
    cent_min_dist = float("inf")
    voxel = None
    for lines in text:
        for line in lines[1].split("\n"):
            split = line.split(",")
            if len(split) == 4:
                min_ra, max_ra, min_dec, max_dec = split
                voxel = Voxel(float(min_ra), float(max_ra), float(min_dec), float(max_dec))
                geometric_centroid_ra, geometric_centroid_dec = voxel.getVoxelCentroid()
                root = Node(voxel)
            elif line:
                border = False if split[13].lower() == "false" else True

                star = Element(int(split[0]), float(split[1]), float(split[2]), float(split[3]), float(split[4]),
                               float(split[5]), float(split[6]), float(split[7]), float(split[8]),
                               float(split[9]), float(split[10]), float(split[11]), float(split[12]), 0, border)
                root.addElement(star)

                if star.border is False:
                    dist = EucDist(star.getRa(), geometric_centroid_ra, star.getDec(), geometric_centroid_dec)
                    if dist < cent_min_dist:
                        centroid = star
                        cent_min_dist = dist

    root.setSize(len(root.getElements()))
    root.addCentroid(centroid)

    level = compute_level(voxel.getSideSize(), voxel.getHeightSize(), query.getMinDistGeneralQ())
    tree = QuadTree(root, level)

    print("\n**** Data Descriptions *****")
    print("Sky Voxel: %s,%s,%s,%s" % (voxel.x_left, voxel.x_right, voxel.y_left, voxel.y_right))
    print("Sky Diagonal: %s" % voxel.getDiagonal())
    print("Min Query Distance: %s" % query.getMinDistGeneralQ())
    print("Max Query Distance: %s" % query.getMaxDistance())
    print("Tree Level: %s" % level)
    print("Tree Elements: %s" % root.size)
    print("Tree Leaf nodes: %s" % len(tree.nodes))
    print("**** End Data Descriptions *****\n")

    return [tree]


"""
def produce_candidates(tree):
    query = query_broadcast.value

    query_elements = query.getQuery()
    epsilon_percentage = query.getEpsilon()
    max_dist_general_pairs = query.getMaxDistGeneralPair()
    query_matrix_distance = query.getDistanceMatrix()

    nodes = tree.nodes

    relations = {}
    for i in range(len(nodes)):
        node_i_ra, node_i_dec = nodes[i].getVoxel().getVoxelCentroid()
        if not nodes[i].isEmpty():
            for j in range(i + 1, len(nodes)):
                node_j_ra, node_j_dec = nodes[j].getVoxel().getVoxelCentroid()
                dist_anchor = EucDist(node_i_ra, node_j_ra, node_i_dec, node_j_dec)
                sc = dist_anchor
                if sc <= query.max_scale:
                    epsilon = sc * epsilon_percentage
                    matrix_distance = [[sc * query_matrix_distance[x][y] for x in range(len(query_matrix_distance))]
                                       for y in range(len(query_matrix_distance))]
                    for k in range(len(nodes)):
                        if k != i and k != j:

                            node_k_ra, node_k_dec = nodes[k].getVoxel().getVoxelCentroid()
                            d1 = EucDist(node_i_ra, node_k_ra, node_i_dec, node_k_dec)
                            d2 = EucDist(node_j_ra, node_k_ra, node_j_dec, node_k_dec)

                            for e in query_elements:
                                eid = e.getId()
                                if eid != max_dist_general_pairs[0] and eid != max_dist_general_pairs[1]:
                                    dq_1 = matrix_distance[max_dist_general_pairs[0]][eid]
                                    dq_2 = matrix_distance[max_dist_general_pairs[1]][eid]

                                    if ((sc - epsilon) * dq_1) - epsilon <= d1 <= ((sc + epsilon) * dq_1) + epsilon and \
                                                                    ((sc - epsilon) * dq_2) - epsilon <= d2 <= (
                                                        (sc + epsilon) * dq_2) + epsilon:
                                        produce_star_pairs(nodes[i], nodes[j], nodes[k], dist_anchor, dq_1, dq_2,
                                                           epsilon, eid, sc, relations)
    return relations.values()"""


def produce_candidates_color(tree):
    query = query_broadcast.value
    query_elements = query.getQuery()
    epsilon_percentage = query.getEpsilon()
    max_pairs_i, max_pairs_j = query.getMaxDistGeneralPair()
    query_matrix_distance = query.getDistanceMatrix()
    query_matrix_distance_size = len(query_matrix_distance)

    nodes = tree.nodes

    relations = defaultdict(set)
    for node in nodes:
        anchors_i = node.getElements()
        neighbors = tree.find_neighbors(node, tree.root, query.neighbor_limit, [])
        neighbors_size = len(neighbors)

        if neighbors_size > 1:
            for anchor_i in anchors_i:
                for j in range(neighbors_size):
                    anchor_j = neighbors[j]
                    if anchor_i.pointId != anchor_j.pointId:
                        u_err = anchor_i.u_err
                        g_err = anchor_i.g_err
                        r_err = anchor_i.r_err
                        i_err = anchor_i.i_err
                        z_err = anchor_i.z_err

                        if anchor_i.u - u_err <= anchor_j.u <= anchor_i.u + u_err and \
                                                        anchor_i.g - g_err <= anchor_j.g <= anchor_i.g + g_err and \
                                                        anchor_i.r - r_err <= anchor_j.r <= anchor_i.r + r_err and \
                                                        anchor_i.i - i_err <= anchor_j.i <= anchor_i.i + i_err and \
                                                        anchor_i.z - z_err <= anchor_j.z <= anchor_i.z + z_err:

                            scale = EucDist(anchor_i.getRa(), anchor_j.getRa(), anchor_i.getDec(), anchor_j.getDec())

                            epsilon = scale * epsilon_percentage
                            matrix_minus = [
                                [(scale - epsilon) * query_matrix_distance[x][y] for x in
                                 range(query_matrix_distance_size)]
                                for y in range(query_matrix_distance_size)]

                            matrix_plus = [
                                [(scale + epsilon) * query_matrix_distance[x][y] for x in
                                 range(query_matrix_distance_size)]
                                for y in range(query_matrix_distance_size)]

                            for w in range(neighbors_size):
                                if w != j:
                                    star = neighbors[w]
                                    if star.pointId != anchor_i.pointId and star.pointId != anchor_j.pointId:

                                        if anchor_i.u - u_err <= star.u <= anchor_i.u + u_err and \
                                                                        anchor_i.g - g_err <= star.g <= anchor_i.g + g_err and \
                                                                        anchor_i.r - r_err <= star.r <= anchor_i.r + r_err and \
                                                                        anchor_i.i - i_err <= star.i <= anchor_i.i + i_err and \
                                                                        anchor_i.z - z_err <= star.z <= anchor_i.z + z_err:
                                            dist_i = EucDist(anchor_i.getRa(), star.getRa(), anchor_i.getDec(),
                                                             star.getDec())
                                            dist_a = EucDist(anchor_j.getRa(), star.getRa(), anchor_j.getDec(),
                                                             star.getDec())

                                            for e in query_elements:
                                                eid = e.getId()

                                                if eid != max_pairs_i and eid != max_pairs_j:
                                                    dq_i_min = matrix_minus[max_pairs_i][eid] - (2 * epsilon)
                                                    dq_i_max = matrix_plus[max_pairs_i][eid] + (2 * epsilon)

                                                    if dq_i_min <= dist_i <= dq_i_max:
                                                        dq_a_min = matrix_minus[max_pairs_j][eid] - (2 * epsilon)
                                                        dq_a_max = matrix_plus[max_pairs_j][eid] + (2 * epsilon)

                                                        if dq_a_min <= dist_a <= dq_a_max:
                                                            relations[anchor_i.pointId, anchor_j.pointId].add(
                                                                (anchor_i, anchor_j,
                                                                 star, eid, scale))
        node.visited = True
    return relations.values()


def produce_candidates_no_color(tree):
    query = query_broadcast.value
    query_elements = query.getQuery()
    epsilon_percentage = query.getEpsilon()
    max_pairs_i, max_pairs_j = query.getMaxDistGeneralPair()
    query_matrix_distance = query.getDistanceMatrix()
    query_matrix_distance_size = len(query_matrix_distance)

    nodes = tree.nodes

    relations = defaultdict(set)
    for node in nodes:
        anchors_i = node.getElements()
        centroid_node_i = node.getCentroid()
        neighbors = tree.find_neighbors(node, tree.root, query.neighbor_limit, [])
        neighbors_size = len(neighbors)
        if neighbors_size > 1:
            for anchor_i in anchors_i:
                for j in range(neighbors_size):
                    anchor_j = neighbors[j]
                    scale = EucDist(centroid_node_i.getRa(), anchor_j.getRa(), centroid_node_i.getDec(),
                                    anchor_j.getDec())
                    epsilon = scale * epsilon_percentage
                    matrix_minus = [
                        [(scale - epsilon) * query_matrix_distance[x][y] for x in range(query_matrix_distance_size)]
                        for y in range(query_matrix_distance_size)]

                    matrix_plus = [
                        [(scale + epsilon) * query_matrix_distance[x][y] for x in range(query_matrix_distance_size)]
                        for y in range(query_matrix_distance_size)]

                    for w in range(neighbors_size):
                        if j != w:
                            star = neighbors[w]
                            dist_i = EucDist(centroid_node_i.getRa(), star.getRa(), centroid_node_i.getDec(),
                                             star.getDec())
                            dist_a = EucDist(anchor_j.getRa(), star.getRa(), anchor_j.getDec(), star.getDec())

                            for e in query_elements:
                                eid = e.getId()

                                if eid != max_pairs_i and eid != max_pairs_j:
                                    dq_i_min = matrix_minus[max_pairs_i][eid]
                                    dq_i_max = matrix_plus[max_pairs_i][eid]

                                    if dq_i_min - (epsilon * 2) <= dist_i <= dq_i_max + (epsilon * 2):
                                        dq_a_min = matrix_minus[max_pairs_j][eid]
                                        dq_a_max = matrix_plus[max_pairs_j][eid]

                                        if dq_a_min - (epsilon * 2) <= dist_a <= dq_a_max + (epsilon * 2):
                                            relations[anchor_i.pointId, anchor_j.pointId].add((anchor_i, anchor_j, star,
                                                                                               eid, scale))
        node.visited = True

    return relations.values()


def prepare_matching_pairs(data):
    query = query_broadcast.value

    dict_rel_id = dict()
    relations = Relations()
    query_l = query.getQuery()
    q_anchor_i, q_anchor_j = query.getMaxDistGeneralPair()

    for i in range(len(query_l)):
        if q_anchor_i != query_l[i].getId() and q_anchor_j != query_l[i].getId():
            relations.addRelation(Relation(query_l[i].getId(), []))
            dict_rel_id[query_l[i].getId()] = len(relations.getRelations()) - 1

    for anchor_1, anchor_2, partner, eid, scale in data:
        id_rel = dict_rel_id[eid]
        relation1 = relations.getRelations()
        relation1[id_rel].addStar(anchor_1, anchor_2, partner, scale)

    return relations


def filter_candidates(relations):
    query = query_broadcast.value
    ee = Plan(query)
    filter_operation = True
    plan = ee.build_plan(relations, filter_operation)
    solution_list = []
    if plan is not None:
        total_solutions = 0
        result = plan.getNext()
        while result != "end":
            if len(result) > 0:
                solution_list.extend(result)
                print("\n**************** SOLUTIONS ************")
                while result != "end":
                    for r in result:
                        total_solutions += 1
                        print("\nScale: %s" % r.getScale())
                        for i in range(len(r.getStars())):
                            if i == 0:
                                print("Anchor:\tId: %s\tra: %s\tDec: %s" % (r.getStar(i).getId(), r.getStar(i).getRa(),
                                                                           r.getStar(i).getDec()))
                            else:
                                print ("Elem Relation: %s\tId: %s\tra: %s\tDec: %s" % (r.getMetadataPos(i),
                                                                                       r.getStar(i).getId(),
                                                                                       r.getStar(i).getRa(),
                                                                                       r.getStar(i).getDec()))
                    result = plan.getNext()
                    if result != "end":
                        solution_list.extend(result)
                print("\n**************** SOLUTIONS ************")
                print("Solution Shapes: Total: %s\n" % total_solutions)

            else:
                result = plan.getNext()
    else:
        print("No Plan")
    print("End Approx based solution")

    return solution_list


def check_solutions(solutions):
    results = []
    query = query_broadcast.value
    matrix_distance = query.getDistanceMatrix()
    size = len(matrix_distance)

    if len(solutions) > 0:
        for r in solutions:
            valid = True
            scale_factor = r.getScale()
            calculated_epsilon = scale_factor * query.getEpsilon()
            matrix_minus = [[(scale_factor - calculated_epsilon) * matrix_distance[x][y] for x in range(size)] for y in
                            range(size)]
            matrix_plus = [[(scale_factor + calculated_epsilon) * matrix_distance[x][y] for x in range(size)] for y in
                           range(size)]

            r_size = len(r.getStars())

            for i in range(r_size):
                metadata_i = r.getMetadataPos(i)
                star_i = r.getStar(i)

                for j in range(i + 1, r_size):
                    metadata_j = r.getMetadataPos(j)
                    star_j = r.getStar(j)

                    if star_i.getId() == star_j.getId():
                        valid = False
                        break
                    else:
                        dist = EucDist(star_i.getRa(), star_j.getRa(), star_i.getDec(), star_j.getDec())

                        if matrix_minus[metadata_i][metadata_j] - calculated_epsilon <= dist <= \
                                        matrix_plus[metadata_i][metadata_j] + calculated_epsilon:
                            pass
                if valid is False:
                    break

            if valid is True:
                results.append(r)

        print("True solution = %s/%s" % (len(results), len(solutions)))
    return results


def save_result(solution_list, filename, query):

    def check_solution(query, solution):
        for e in solution.getStars():
            if e.getRa() < 0.0:
                e.Ra = 360 - e.getRa()
            else:
                if e.getRa() > 360.0:
                    e.Ra = 0 + (e.getRa() - 360)

        star_i = solution.getStar(0)
        star_j = solution.getStar(1)

        scale = EucDist(star_i.getRa(), star_j.getRa(), star_i.getDec(), star_j.getDec()) / query.getMaxDistance()
        solution.scale = scale

    if len(solution_list) > 0:
        with open(filename, "w") as f:
            count = 1
            for solution in solution_list:
                if solution:
                    for r in solution:
                        check_solution(query, r)

                        f.write("Solution = %s\n" % count)
                        f.write("Scale = %s\n" % r.getScale())

                        for element in r.getStars():
                            f.write("Star = %s\n" % element)
                        f.write("\n")
                        count += 1


if __name__ == '__main__':
    start_time = time()
    if len(argv) < 4:
        print(
            "Isotropic Constellations Queries:\n[1]InputFile\n[2]OutputFile\n[3]Epsilon Percentage\n[4]Neighbor Limit")
    else:

        spark_conf = SparkConf()
        spark_conf.setAppName("Isotropic Constellations Queries")
        spark_conf.setMaster("local[1]")  # yarn-client
        # spark_conf.set("spark.executor.instances", "120")
        # spark_conf.set("spark.executor.cores", "1")
        # spark_conf.set("spark.driver.memory","24g")
        # spark_conf.set("spark.executor.memory","4g")
        # spark_conf.set("spark.yarn.queue", "prioridade")

        sc = SparkContext(conf=spark_conf)
        dir_path = path.dirname(path.abspath(__file__))

        sc.addPyFile(dir_path + "/cost.py")
        sc.addPyFile(dir_path + "/element.py")
        sc.addPyFile(dir_path + "/node.py")
        sc.addPyFile(dir_path + "/quadtree.py")
        sc.addPyFile(dir_path + "/query.py")
        sc.addPyFile(dir_path + "/util.py")
        sc.addPyFile(dir_path + "/voxel.py")
        sc.addPyFile(dir_path + "/executionengine.py")

        # sc.addFile(path.dirname(path.abspath(__file__)), True)
        query = Query.defineQuery(float(argv[3]), neighbor_limit=float(argv[4]))

        query_broadcast = sc.broadcast(query)

        text_rdd = sc.newAPIHadoopFile(argv[1], "org.apache.hadoop.mapreduce.lib.input.TextInputFormat",
                                       "org.apache.hadoop.io.Text", "org.apache.hadoop.io.LongWritable",
                                       conf={"textinputformat.record.delimiter": "\n\n"})

        tree_rdd = text_rdd.mapPartitions(build_tree)

        candidates_rdd = tree_rdd.flatMap(produce_candidates_color)

        relations_rdd = candidates_rdd.map(prepare_matching_pairs)

        filter_rdd = relations_rdd.map(filter_candidates)

        # solution_rdd = filter_rdd.map(check_solutions)

        result = filter_rdd.collect()

        sc.stop()

        save_result(result, argv[2], query)

        end_time = time() - start_time
        print("Total Time: %s" % end_time)
