distanceMatrix = [[0 ,12,12,13,15,15],
                  [12, 0, 2, 6, 8, 8],
                  [12, 2, 0, 6, 9, 9],
                  [13, 6, 6, 0, 8, 8],
                  [15, 8, 9, 8, 0, 4],
                  [15, 8, 9, 8, 4, 0]]
speciesList = ["M_Spacii", "T_Pain", "G_Unit", "Q_Doba", "R_Mani", "A_Finch"]

def UPGMA(dM,sp):

    while (len(dM)> 1):
        leastRow,LeastCol = findSmallest(dM)    # finds the smallest non-0 matrix coordinate
        dM = updateMatrix(dM,leastRow,LeastCol) # This is the function you are writing. See bottom of cell
        sp = updateSpecies(sp,leastRow,LeastCol)# updates the species list
        print(dM)
        print(speciesList)
    
def findSmallest(dM):# finds the smallest non-0 matrix coordinate

    #Search the matrix for the coordinate of the shortest distance between organisms
    #YOUR CODE GOES HERE
    
    return row, col # Returns row and column position of the smallest distance in the matrix

def updateSpecies(sp,r,c): # updates the species list
    sp[r]="("+sp[r]+","+sp[c]+")"
    del sp[c]
    return sp
   
def updateMatrix(dM,row,col): #the input is the previous stage of the matrix and the row and column of the minimum 
                              #distance pair in the matrix. These are the columns that should be replaced in the new matrix
                              #for the adventurous:https://en.wikipedia.org/wiki/UPGMA
    newMat = []
    
    #complete this function
    #YOUR CODE GOES HERE
    #follow the algorithm in the slide deck to create a new matrix (newMat) that has the approproate distances
    
    return newMat

UPGMA(distanceMatrix,speciesList)
