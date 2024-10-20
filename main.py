import numpy as np

def is_complete(label_matrix): #checks if all entries in label matrix are covered.
    for element in np.nditer(label_matrix):
        if element == 0:
            return False
    return True

def operations(M): #core calculations for the rref,determinant, and inverse functions below.
    M = M.astype(float)
    A = M.copy() #A is the matrix to be converted to RREF
    # label_matrix is responsible for checking if each entry is covered or not.
    label_matrix = np.zeros_like(A)
    identity_matrix = np.identity(4)
    possible_inverse_matrix = identity_matrix #inverse starts as identity then is reduced to the inverse if det != 0


    #track changes from Type-1 and Type-2 EROS to calculate the determinant later.
    row_swap_num = 0
    determinant = 1

    # head is defined to be the index uppermost element in a "non-covered" column.
    while not is_complete(label_matrix):
        row, col = label_matrix.shape
        found_head = False

        # Find the head (the first non-covered element in a non-covered column)
        for j in range(col):
            for i in range(row):
                if label_matrix[i, j] == 0:
                    head_x, head_y = i, j
                    found_head = True
                    break
            if found_head:
                break

        # Find the first non-zero element in the column
        for tail_x in range(head_x, row):
            if A[tail_x, head_y] != 0:
                break

        if A[tail_x, head_y] == 0:
            label_matrix[:, head_y] = 1  # Mark the entire column as covered
        else:
            # Type-1 REF: switch rows head_x and tail_x
            if head_x != tail_x:
                A[[head_x, tail_x]] = A[[tail_x, head_x]]
                possible_inverse_matrix[[head_x, tail_x]] = possible_inverse_matrix[[tail_x, head_x]]
                row_swap_num += 1

            # Type-2 REF: normalize the row so that the leading coefficient is 1
            determinant = determinant * A[head_x, head_y]
            possible_inverse_matrix[head_x] = possible_inverse_matrix[head_x] / A[head_x, head_y]
            A[head_x] = A[head_x] / A[head_x, head_y]

            # Type-3 REF: eliminate all other entries in the column
            for i in range(row):
                if i != head_x:
                    possible_inverse_matrix[i] = possible_inverse_matrix[i] - A[i, head_y] * possible_inverse_matrix[head_x]
                    A[i] = A[i] - A[i, head_y] * A[head_x]

            # Update determinant based on number of row swaps.
            if(row_swap_num % 2 != 0): #odd
                determinant = -determinant

            # Mark the row and column as covered
            label_matrix[:, head_y] = 1
            label_matrix[head_x, :] = 1

    # Multiply the gathered constants in determinant calculation by the main diagonal since it is now an upper triangle
    for i in range(row):
        for j in range(col):
            if i == j:
                determinant = determinant*A[i,j]

    return A, determinant, possible_inverse_matrix

def RREF(M):
    rref,_,_ = operations(M)
    return rref
def determinant(M):
    _,determinant,_ = operations(M)
    return determinant
def inverse(M):
    if determinant(M) == 0:
        print("Matrix is Singular, therefore it is not invertible.")
    else:
        _,_,inverse_matrix = operations(M)
        return inverse_matrix

        

# Main function
matrixA = np.zeros((4,4))
print("Enter a 4x4 matrix: ")
for i in range(4):
    for j in range(4):
        matrixA[i,j] = input(f"Entry at position {i+1,j+1}:\n")


print(f"RREF is \n{RREF(matrixA)}\n")
print(f"Determinant is \n{determinant(matrixA)}\n")
print(f"Inverse is \n{inverse(matrixA)}\n")



#---------the code below is for testing purposes -----------#
1
#you can test rref using an online rref calculator or use sympy library.
# print(f"Determinant is \n{np.linalg.det(matrixA)}\n")
# print(f"Inverse is \n{np.linalg.inv(matrixA)}\n")
