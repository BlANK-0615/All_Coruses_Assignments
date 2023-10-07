# B224911-Script of assignment for Bioinformatics Algorithms 
# 03-03-2013 version 5
# Source code from Lecture04 tutorial of Bioinformatics Algorithms
# The forlume of this script form the https://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf
# Perform Affine gaps Alignment in Python (from first principles)

# step1 create three matrix
# contains rows lists each of length cols initially set to 0
# index as my_matrix[1][2] my_matrix[R][C]
def create_matrix(rows, cols):
    matrix_m = [[0 for col in range(cols+1)] for row in range(rows+1)]
    matrix_x = [[0 for col in range(cols+1)] for row in range(rows+1)]
    matrix_y = [[0 for col in range(cols+1)] for row in range(rows+1)]
    for i in range(0,len(matrix_x[0])):
        matrix_x[0][i] = -999
        matrix_y[0][i] = -999   
    for j in range(0,len(matrix_x)):
        matrix_x[j][0] = -999
        matrix_y[j][0] = -999
    return matrix_m, matrix_x, matrix_y


# Step2 clalculate these three matrix
# x is row index, y is column index
# follows[r][c]
def calc_score(matrix_M, matrix_IX, matrix_IY, x, y):
    # print("seq1:",sequence1[y- 1]," seq2: "+sequence2[x - 1],"x:",x," y:",y)
    # x,y is the position of the seuquence.
    # match and mis match for alignment matrix
    sc = seqmatch if sequence1[y - 1] == sequence2[x - 1] else seqmismatch
    base_score_m = matrix_M[x - 1][y - 1] + sc
    insert_score_m = matrix_IX[x - 1][y - 1] + sc
    delete_score_m = matrix_IY[x - 1][y - 1] + sc
    vm=max(0, base_score_m, insert_score_m, delete_score_m) 
    # Insertion or Deletion_X
    base_score_x = matrix_M[x - 1][y] + opening + seqgap
    insert_score_x = matrix_IX[x - 1][y] + seqgap
    vx=max(base_score_x, insert_score_x)
    # Insertion or Deletion_Y
    base_score_y = matrix_M[x][y - 1] + opening + seqgap
    insert_score_y = matrix_IY[x][y - 1] + seqgap
    vy=max(base_score_y, insert_score_y)
    return vm, vx, vy


# step3 makes a single traceback step
def traceback(matrix_M, matrix_IX, matrix_IY, maxv, matrix_return):
    x=maxv[0]
    y=maxv[-1]
    if matrix_return == "M":
        val=matrix_M[x][y]
        sc = seqmatch if sequence1[y - 1] == sequence2[x - 1] else seqmismatch    
        base_score_m = matrix_M[x - 1][y - 1] + sc
        if  val == base_score_m:
            return [x-1,y-1], "M"
        insert_score_m = matrix_IX[x - 1][y - 1] + sc
        if val == insert_score_m:
            return [x-1,y-1], "X"
        delete_score_m = matrix_IY[x - 1][y - 1] + sc
        if val == delete_score_m:
            return [x-1,y-1], "Y"
    elif matrix_return == "X":
        val = matrix_IX[x][y]
        base_score_x = matrix_M[x - 1][y] + opening + seqgap
        if val == base_score_x:
            return [x-1, y], "M"
        insert_score_x = matrix_IX[x - 1][y] + seqgap
        if val == insert_score_x:
            return [x-1, y], "X"
    elif matrix_return == "Y":
        val = matrix_IY[x][y]
        base_score_y = matrix_M[x][y - 1] + opening + seqgap
        if val == base_score_y:
            return [x, y-1], "M"
        insert_score_y = matrix_IY[x][y - 1] + seqgap
        if val == insert_score_y:
            return [x, y-1], "Y" 
# It should be noted here that the coordinates are returned, the syntax is [x,y], and the list slice syntax is list[x][y]



# step4 Start calculations for three, and fill the calculated values ​​into the three matrices
# builds the initial scoring matrix used for traceback
def build_matrix(matrix_M, matrix_IX, matrix_IY):
    rows=len(matrix_M)
    cols=len(matrix_M[0])
    # because the cols in the second layer of the matrix, so you need to have len # the first factor  
    # >>> my_matrix[0]
    # [0, 0, 0, 0, 0, 0]
    for i in range(1, rows):
        for j in range(1, cols):
            matrix_M[i][j], matrix_IX[i][j], matrix_IY[i][j] = calc_score(matrix_M, matrix_IX, matrix_IY, i, j)
    return matrix_M, matrix_IX, matrix_IY

# step5 gets the max value from the built matrix
def get_max(matrix_M):
    max=matrix_M[0][0]
    mrow=0
    mcol=0
    rows = len(matrix_M)
    cols = len(matrix_M[0])
    for i in range(1, rows):
        for j in range(1, cols):
            if matrix_M[i][j]>max:
                max=matrix_M[i][j]
                mrow=i
                mcol=j
    print("max score in matrix_M: ",max)
    return [mrow,mcol] # rerurn the max position
# View the maximum score and its position

# step6 print out the best scoring path from the SW matrix
def print_matrix(matrix_M, matrix_IX, matrix_IY):
    rows = len(matrix_M)
    cols = len(matrix_M[0])
    s1="  " +sequence1
    s2=" "+sequence2
    print("Dimensions: r= %2d , c= %2d" % (rows, cols))
    # 02d means an integer, left padded with zeros up to 2 digits.
    # These are simply the way of printing the numbers. Look at this post java specifiers to understand what specifiers you are using.
    print("\n\n--------This is matrix_M for alignment matrix--------")
    for a in s1:
        print(a, end ="")
        print(" \t", end ="")
    print("\n",end="")
    for i in range(0, rows):
        print(s2[i], end ="")
        print(" \t", end ="")
        for j in range(0, cols):
            
            print("%02d\t" % (matrix_M[i][j]),end="")
        print("\n",end="") 
    print("\n\n--------This is matrix_IX for alignment insertion/deletion in X of seq2--------")
    for a in s1:
        print(a, end ="")
        print(" \t", end ="")
    print("\n",end="")
    for i in range(0, rows):
        print(s2[i], end ="")
        print(" \t", end ="")
        for j in range(0, cols):
            
            print("%02d\t" % (matrix_IX[i][j]),end="")
        print("\n",end="") 
    print("\n\n--------This is matrix_IY for insertion/deletion in Y of seq1--------")
    for a in s1:
        print(a, end ="")
        print(" \t", end ="")
    print("\n",end="")
    for i in range(0, rows):
        print(s2[i], end ="")
        print(" \t", end ="")
        for j in range(0, cols):
            
            print("%02d\t" % (matrix_IY[i][j]),end="")
        print("\n",end="") 
    print("\n\n")


#step7 rint out the traceback of the best scoring alignment
def print_traceback(matrix_M):
    #this will print as expected with internal gaps
    print("Building traceback...")
    maxv=get_max(matrix_M)
    print(maxv)
    #traverse the matrix to find the traceback elements
    #if more than one path just pick one
    topstring=""
    midstring=""
    bottomstring=""
    #pad the sequences so indexes into the sequences match the matrix indexes
    asequence1 = "#" + sequence1
    asequence2 = "#" + sequence2
    old_maxv=maxv
    #add the rest of the elements
    search=True
    matrix_return ="M"
    while(search):
        maxv, matrix_return=traceback(matrix_M, matrix_IX, matrix_IY, maxv, matrix_return)
        # print(matrix_return) check out the matrix jump
        # print(old_maxv, maxv) check out the position 
        # insertion or deletion
        if(old_maxv[-1]==maxv[-1]):
            topstring+="-"
            bottomstring += asequence2[old_maxv[0]]
            #up to cell and only the row changes 
            #mean the cols,the seq1,the top has a gap
        elif(old_maxv[0] == maxv[0]):
            topstring += asequence1[old_maxv[-1]]
            bottomstring += "-"
            #left to cell and only the col changes 
            #mean the rows,the seq2,the bottom has a gap
        else:
            # add the next element and go to previous
            topstring += asequence1[old_maxv[-1]]
            bottomstring += asequence2[old_maxv[0]]
        if (asequence1[old_maxv[-1] ] == asequence2[old_maxv[0] ]) & (old_maxv[-1]!=maxv[-1]) & (old_maxv[0] != maxv[0]):
            midstring += "|"
        else:
            midstring += " "
        if(matrix_M[maxv[0]][maxv[-1]]==0):
            search=False
            continue #this is to continue the while loop and jump out
        old_maxv = maxv
    print(topstring[::-1])
    print(midstring[::-1])
    print(bottomstring[::-1])


#step8 build the SW alignment...
def perform_smith_waterman():
#values for weights
    global  seqmatch
    global  seqmismatch
    global  seqgap
    global  opening
    global sequence1
    global sequence2

# note these are not the exact weights used in the original SW paper
seqmatch = 4
seqmismatch = -1
seqgap = -1
opening = -2

# input sequences
# test sequence1
sequence1="AGTGATAAACTAGTAATTTTT"
sequence2="TTGGGGGTAAACAGGGG"

# test sequence2
# sequence1 = "GTGTATTTTTTT"
# sequence2 = "AAAAGTGTTATT"


# another test sequence
# sequence1 = "TCGTAGACGA"
# sequence2 = "ATAGAATGCGG"

print("Sequence1: "+sequence1);
print("Sequence2: "+sequence2);

# run the main funciton
matrix_M, matrix_IX, matrix_IY=create_matrix(len(sequence2), len(sequence1))
matrix_M, matrix_IX, matrix_IY=build_matrix(matrix_M, matrix_IX, matrix_IY)
print_matrix(matrix_M, matrix_IX, matrix_IY)
print_traceback(matrix_M)
# this calls the SW algorithm when the script loads
perform_smith_waterman()

