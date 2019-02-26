from math import sqrt, pow, ceil, log


def EucDist(x1, x2, y1, y2):
    """
    Distance Function: Euclidean

    This is already implemented by scipy.spatial.distance.euclidean(u, v)[source]
    http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.distance.euclidean.html
    """
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2))


def compute_level(l1, l2, distance):
    """

    :param l1:
    :param l2:
    :param distance:
    :return:
    """
    return ceil(log((sqrt((pow(l1, 2) + pow(l2, 2))) / distance), 2))


def calculate_distance(voxel_node1, voxel_node2):
    point_of_voxel_node1 = [[voxel_node1.x_left, voxel_node1.y_left],
                            [voxel_node1.x_right, voxel_node1.y_left],
                            [voxel_node1.x_right, voxel_node1.y_right],
                            [voxel_node1.x_left, voxel_node1.y_right]]

    point_of_voxel_node2 = [[voxel_node2.x_left, voxel_node2.y_left],
                            [voxel_node2.x_right, voxel_node2.y_left],
                            [voxel_node2.x_right, voxel_node2.y_right],
                            [voxel_node2.x_left, voxel_node2.y_right]]

    min_distance = float("inf")
    max_distance = float("-inf")
    for data_node1 in point_of_voxel_node1:
        for data_node2 in point_of_voxel_node2:
            calculated_distance = EucDist(data_node1[0], data_node2[0], data_node1[1], data_node2[1])
            if calculated_distance < min_distance:
                min_distance = calculated_distance
            if calculated_distance > max_distance:
                max_distance = calculated_distance

    return min_distance, max_distance
