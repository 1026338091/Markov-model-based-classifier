import random

# DM and
distanceMatrix = [[0, 12, 12, 13, 15, 15],
                  [12, 0, 2, 6, 8, 8],
                  [12, 2, 0, 6, 9, 9],
                  [13, 6, 6, 0, 8, 8],
                  [15, 8, 9, 8, 0, 4],
                  [15, 8, 9, 8, 4, 0]]

speciesList = ["M_Spacii", "T_Pain", "G_Unit", "Q_Doba", "R_Mani", "A_Finch"]

# A global dictionary for storing the height of nodes
node_heights = {}

"""
Question 1 a
"""
def findSmallest(dM):
    """
    :param dM: distance matrix
    :return: row and cloumn index of samallest value
    """
    min_dist = float('inf')
    min_positions = []
    for i in range(len(dM)):
        for j in range(i + 1, len(dM[i])):
            if dM[i][j] > 0 and dM[i][j] < min_dist:
                min_dist = dM[i][j]
                min_positions = [(i, j)]

            elif dM[i][j] == min_dist:
                min_positions.append((i, j))

    if min_positions:
        row, col = random.choice(min_positions)
        return row, col

    return 0, 1

"""
Question 1 b

Update the distance matrix: Merge the species corresponding to row and col 
Using the UPGMA algorithm: The distance from a new node to other nodes is the average of the distances between the original two nodes
"""
def updateMatrix(dM, row, col):
    n = len(dM)
    newMat = []
    if row > col:
        row, col = col, row
    newMat = []

    for i in range(n):
        if i == col:
            continue

        new_row = []

        for j in range(n):
            if j == col:
                continue

            if i == row and j == row:

                new_row.append(0)
            elif i == row:
                # new vale
                avg_dist = (dM[row][j] + dM[col][j]) / 2.0
                new_row.append(avg_dist)
            elif j == row:
                # new vale
                avg_dist = (dM[i][row] + dM[i][col]) / 2.0
                new_row.append(avg_dist)
            else:
                # other position
                new_row.append(dM[i][j])

        newMat.append(new_row)

    return newMat

def updateSpecies(sp, r, c):
    """
    update the special
    """
    sp[r] = "(" + sp[r] + "," + sp[c] + ")"
    del sp[c]
    return sp

def UPGMA(dM, sp):
    """
    UPGMA main function
    """
    global node_heights

    # 初始化所有叶节点的高度为0
    for species in sp:
        node_heights[species] = 0

    print("Starting UPGMA Algorithm")
    print("=" * 60)

    step = 1

    while len(dM) > 1:
        # find least distance
        leastRow, leastCol = findSmallest(dM)

        # obtain species name
        species1 = sp[leastRow]
        species2 = sp[leastCol]

        # obtain distance
        merge_distance = dM[leastRow][leastCol]
        node_height = merge_distance / 2.0

        # caculate branch
        branch1 = node_height - node_heights[species1]
        branch2 = node_height - node_heights[species2]

        # creat new node
        new_node = f"({species1}:{branch1:.1f},{species2}:{branch2:.1f})"

        # the height of new node
        node_heights[new_node] = node_height

        # print
        print(f"\nStep {step}:")
        print(f"  Merging: {species1} and {species2}")
        print(f"  Distance: {merge_distance:.1f}")
        print(f"  Node height: {node_height:.1f}")
        print(f"  Branch lengths: {branch1:.1f}, {branch2:.1f}")

        # new dM and sp
        dM = updateMatrix(dM, leastRow, leastCol)
        sp = updateSpecies(sp, leastRow, leastCol)
        sp[leastRow] = new_node

        print(f"  Current tree structure: {sp[0]}")

        step += 1

    print("\n" + "=" * 60)
    print("Final Tree (Newick format):")
    print(sp[0] + ";")

    return sp[0]


# UPGMA
if __name__ == "__main__":

    dM_copy = [row[:] for row in distanceMatrix]
    sp_copy = speciesList[:]

    final_tree = UPGMA(dM_copy, sp_copy)

    # save
    with open("HongyuanDeng_UPGMA.txt", "w") as f:
        f.write("UPGMA Algorithm Output\n")
        f.write("=" * 60 + "\n\n")
        f.write("Input Data:\n")
        f.write(f"Species: {speciesList}\n")
        f.write("Distance Matrix:\n")
        for row in distanceMatrix:
            f.write(str(row) + "\n")
        f.write("\n" + "=" * 60 + "\n")
        f.write("Final Tree (Newick format):\n")
        f.write(final_tree + ";\n")



















