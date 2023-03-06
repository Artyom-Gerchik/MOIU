# example from doc: x1 + x2 -> max
# (basis index : 3,4,5; non-basis: 1,2)


import numpy

# initial data
import lab1_help

A = numpy.array([[-1, 1, 1, 0, 0], [1, 0, 0, 1, 0], [0, 1, 0, 0, 1]])
C_with_T = numpy.array([1, 1, 0, 0, 0])
X_with_T = numpy.array([0, 0, 1, 3, 2])
B = numpy.array([2, 3, 4])  # LESS BY 1 BECAUSE ZERO-BASED ARRAYS (in task = 3,4,5)
# initial data

A_B = numpy.zeros((B.size, B.size))
A_B_INV = numpy.zeros((B.size, B.size))
C_T_B = numpy.zeros(B.size)
U_T = numpy.zeros(C_T_B.size)
TRIANGLE_T = numpy.zeros(C_with_T)
j_0 = 0  # index of first negativ on triangle_t
z = numpy.zeros(len(A))
TETHA = numpy.zeros(len(z))
TETHA_MIN_INDEX = 0
TETHA_MIN = 0


def main():
    for iteration in range(0, 100):
        print()
        print('ITERATION', iteration + 1, 'ITERATION')
        print('****************FIRST STEP****************')

        for index in range(len(A)):
            for jindex in range(len(A)):
                A_B[index][jindex] = A[index][B[jindex]]
        print('A_B =\n', A_B)

        print('****************FIRST STEP****************')

        print()

        print('****************SECOND STEP***************')

        if iteration == 0:

            A_B_INV = numpy.linalg.inv(A_B)
        else:
            tmp = numpy.zeros(B.size)

            for index in range(0, B.size):
                tmp[index] = A_B[index][TETHA_MIN_INDEX]

            A_B_INV = lab1_help.func_for_lab2(A_B_INV, tmp,
                                              TETHA_MIN_INDEX + 1,
                                              B.size)  # j_0 + 1 because in lab1 all indexes are substructed by 1

        print('A_B_INV =\n', A_B_INV)

        print('****************SECOND STEP***************')

        print()

        print('****************THIRD STEP****************')

        for index in range(B.size):
            C_T_B[index] = C_with_T[B[index]]

        print('C_T_B =\n', C_T_B)

        print('****************THIRD STEP****************')

        print()

        print('****************FOURTH STEP***************')

        U_T = numpy.dot(C_T_B, A_B_INV)

        print('U_T =\n', U_T)

        print('****************FOURTH STEP***************')

        print()

        print('****************FIFTH STEP****************')

        TRIANGLE_T = numpy.subtract(numpy.dot(U_T, A), C_with_T)

        print('TRIANGLE_T =\n', TRIANGLE_T)

        print('****************FIFTH STEP****************')

        print()

        print('****************SIXTH STEP****************')

        counter_tmp = 0

        for index in range(0, TRIANGLE_T.size):
            if TRIANGLE_T[index] < 0:
                j_0 = index
                break
            else:
                counter_tmp += 1

            if counter_tmp == len(TRIANGLE_T):
                print('NOW IT S AN OPTIMAL PLAN : ', X_with_T)
                return

        print('j_0 = ', j_0)

        print('****************SIXTH STEP****************')

        print()

        print('****************SEVENTH STEP**************')

        tmp = numpy.zeros(3)

        for index in range(0, len(A)):
            tmp[index] = A[index][j_0]

        z = numpy.dot(A_B_INV, tmp)

        print('z = ', z)

        print('****************SEVENTH STEP**************')

        print()

        print('****************EIGHTH STEP***************')

        for index in range(0, len(z)):
            if z[index] > 0:
                TETHA[index] = X_with_T[B[index]] / z[index]
                TETHA_MIN_INDEX = index
                TETHA_MIN = X_with_T[B[index]] / z[index]
            else:
                TETHA[index] = 999
        print("TETHA = ", TETHA)
        print("TETHA_MIN = ", TETHA_MIN)  # REMEMBER ABOUT ZERO-BASED

        print('****************EIGHTH STEP***************')

        print()

        print('****************NINTH STEP****************')

        B[TETHA_MIN_INDEX] = j_0
        print('NEW B IS : ', B)  # REMEMBER ABOUT ZERO-BASED

        for index in range(0, len(B)):
            if index == TETHA_MIN_INDEX:
                X_with_T[B[index]] = TETHA_MIN
            else:
                X_with_T[B[index]] = X_with_T[B[index]] - (TETHA_MIN * z[index])

        for index in range(0, len(X_with_T)):
            if index not in B:
                X_with_T[index] = 0

        print('NEW X_with_T IS: ', X_with_T)

        print('****************NINTH STEP****************')
        print()
        print('ITERATION', iteration + 1, 'ITERATION')
        print()


if __name__ == '__main__':
    main()
