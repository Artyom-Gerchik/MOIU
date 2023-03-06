import copy

import numpy


def initial_step():
    matrix = numpy.array([[3, 5, 7], [6, 3, 8], [9, 5, 6]])
    matrix_inverted = numpy.linalg.inv(matrix)
    x_vector_column = numpy.array([6, 7, 8])
    index = 2
    n_size = 3
    matrix_with_overline = copy.deepcopy(matrix)
    for i in range(0, n_size):
        matrix_with_overline[i][index - 1] = x_vector_column[i]
    return matrix, matrix_inverted, x_vector_column, index, n_size, matrix_with_overline


def first_step(matrix_inverted, x_vector_column, index):
    l_vector_column = numpy.dot(matrix_inverted, x_vector_column)
    if l_vector_column[index - 1] == 0:
        print('Matrix A^- are irreversible')
        return 0
    else:
        return l_vector_column


def second_step(l, index):
    l[index - 1] = -1
    return l


def third_step(l, l_with_overline, index):
    l_with_roof = ((-1) / l[index - 1]) * l_with_overline
    return l_with_roof


def fourth_step(n_size, l_with_roof, index):
    q = numpy.eye(n_size)
    for i in range(0, n_size):
        q[i][index - 1] = l_with_roof[i]

    return q


def fifth_step(n_size, index, matrix_q, matrix_inverted):
    final_matrix = numpy.zeros((n_size, n_size))
    for i in range(0, n_size):
        for k in range(0, n_size):
            if index - 1 != i:
                final_matrix[i][k] = matrix_q[i][i] * matrix_inverted[i][k] + matrix_q[i][index - 1] * \
                                     matrix_inverted[index - 1][k]
            else:  # going here, if got index column. that means that we had 1 element.
                final_matrix[i][k] = matrix_q[i][index - 1] * matrix_inverted[index - 1][k]
    return final_matrix


def main():
    print('***************INITIAL STEP***************')
    matrix, matrix_inverted, x_vector_column, index, n_size, matrix_with_overline = initial_step()
    print('A = \n', matrix, '\n')
    print('A^-1 = \n', matrix_inverted, '\n')
    print('x =', x_vector_column, '\n')
    print('i =', index, '\n')
    print('n = ', n_size, '\n')
    print('A^- = \n', matrix_with_overline)

    print('***************INITIAL STEP***************')
    print()
    print('****************FIRST STEP****************')
    l = first_step(matrix_inverted, x_vector_column, index)
    if l.size == 0:
        return
    else:
        print('l =', l)
    print('****************FIRST STEP****************')
    print()
    print('****************SECOND STEP***************')
    l_with_overline = second_step(copy.deepcopy(l), index)
    print('l^- =', l_with_overline)
    print('****************SECOND STEP***************')
    print()
    print('****************THIRD STEP****************')
    l_with_roof = third_step(copy.deepcopy(l), copy.deepcopy(l_with_overline), index)
    print('l^^ =', l_with_roof)
    print('****************THIRD STEP****************')
    print()
    print('****************FOURTH STEP***************')
    matrix_q = fourth_step(n_size, copy.deepcopy(l_with_roof), index)
    print('q = \n', matrix_q)
    print('****************FOURTH STEP***************')
    print()
    print('****************FIFTH STEP****************')
    final_matrix = fifth_step(n_size, index, copy.deepcopy(matrix_q), copy.deepcopy(matrix_inverted))
    print('ANSWER = \n', final_matrix)
    print('****************FIFTH STEP****************')

    print('')
    print('TESTTESTTESTTESTTESTTESTTEST')
    print(numpy.linalg.inv(matrix_with_overline))
    print('TESTTESTTESTTESTTESTTESTTEST')
    print('')


if __name__ == '__main__':
    main()
