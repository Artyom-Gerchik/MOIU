import numpy

import lab2_help


def main(c, A, b):
    ########### FIRST STEP ###########
    for index in range(0, b.size):
        if b[index] < 0:
            b[index] *= -1

            for jindex in range(0, A[0].size):
                A[index][jindex] *= 1

    ########### FIRST STEP ###########

    ########### SECOND STEP ###########

    n = A[0].size  # count of columns
    m = int(A.size / n)  # count of rows

    c_with_tilda = numpy.zeros(n + m)

    for index in range(n, n + m):
        c_with_tilda[index] = -1

    ########### SECOND STEP ###########

    ########### THIRD STEP ###########

    x_with_tilda = numpy.zeros(n + m)

    for index in range(n, n + m):
        x_with_tilda[index] = b[index - n]

    B = []

    for index in range(n + 1, n + m + 1):
        B.append(index)

    B_as_numpy_array = numpy.array(B)

    ########### THIRD STEP ###########

    ########### FOURTH STEP ###########

    right_A = numpy.hstack((A, numpy.eye(m)))

    result = lab2_help.help_for_lab3(right_A, c_with_tilda, x_with_tilda, B_as_numpy_array)

    optimal_plan_from_main_phase = numpy.array(result[0])
    new_B_from_main_phase = numpy.array(result[1])

    print(optimal_plan_from_main_phase)
    print(new_B_from_main_phase)

    ########### FOURTH STEP ###########

    ########### FIFTH STEP ###########

    for index in range(n, n + m):
        if optimal_plan_from_main_phase[index] != 0:
            print('Zadacha nesovmestna')
            return

    ########### FIFTH STEP ###########

    ########### SIXTH STEP ###########

    allowable_plan_x = numpy.zeros(n)

    for index in range(0, n):
        allowable_plan_x[index] = x_with_tilda[index]

    ########### SIXTH STEP ###########

    ########### SEVENTH STEP #########

    while True:
        terminate = True

        for index in range(0, len(B)):
            if B[index] > n:
                terminate = False

        if terminate:
            # print('plspls')
            return [allowable_plan_x, B, A, b]

        j_k = max(B)
        k = B.index(j_k)
        l = {}

        A_B = right_A[:, B[0] - 1]

        for index in range(1, len(B)):
            A_B = numpy.vstack((A_B, right_A[:, B[index] - 1]))

        A_B = A_B.transpose()
        A_B_INV = numpy.linalg.inv(A_B)

        for index in range(1, n + 1):
            if not index in B:
                l[index] = A_B_INV.dot(right_A[:, 1])

        zeros = True

        for index in l.keys():
            if l[index][k] != 0:
                B[k] = index
                zeros = False
                break

        if zeros:
            index = j_k - n
            A = numpy.delete(A, index, axis=0)
            b = numpy.delete(b, index - 1, axis=0)
            B.remove(j_k)
            right_A = numpy.delete(right_A, index, axis=0)

    ########### SEVENTH STEP #########


if __name__ == '__main__':
    c = numpy.array([1, 0, 0])  # from pdf vector stoimostei
    A = numpy.array([[1, 1, 1], [2, 2, 2]])  # from pdf matrix ogranicheniy
    b = numpy.array([-10, 0])  # from pdf vector of right chastei

    answer = main(c, A, b)

    print('')
    print('')
    print('ANSWER')
    print('')
    print(f'x_with_T = {answer[0]}')
    print(f'B = {answer[1]}')
    print(f'A = {answer[2]}')
    print(f'b = {answer[3]}')
    print('')
    print('')
