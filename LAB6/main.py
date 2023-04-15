import numpy


def quadratic_programming(x_with_T, J_with_b, J_with_b_with_star, A, D, c):
    iteration = 0

    while True:
        iteration += 1

        print(f'\n{iteration} ITERATION')

        A_B = A[:, J_with_b]
        c_x = c + numpy.dot(D, x_with_T)
        c_B = c_x[J_with_b]
        u_x_transposed = numpy.transpose(numpy.dot(-numpy.transpose(c_B), numpy.linalg.inv(A_B)))
        delta_x = numpy.transpose(numpy.dot(u_x_transposed, A) + numpy.transpose(c_x))

        optimal_plan = False
        tmp = 0

        for i in range(len(delta_x)):
            if delta_x[i] < 0:
                optimal_plan = True
            else:
                tmp += 1
        if optimal_plan and tmp == 0:
            print(f'optimal plan x is: {x_with_T}\n')
            return

        j_0 = i - 1

        l = numpy.zeros((len(x_with_T)))
        l[j_0] = 1

        D_with_star = (D[J_with_b_with_star, :])[:, J_with_b_with_star]
        A_B_with_star = A[:, J_with_b_with_star]

        H = numpy.zeros(((len(D_with_star) + len(A_B_with_star)), (len(D_with_star) + len(A_B_with_star))))

        H[:len(D_with_star), :len(D_with_star)] = D_with_star
        H[:len(A_B_with_star), len(D_with_star):] = numpy.transpose(A_B_with_star)
        H[len(D_with_star):, :len(A_B_with_star)] = A_B_with_star

        b_with_star = numpy.zeros((len(J_with_b_with_star) + len(A)))
        b_with_star[:len(J_with_b_with_star)] = D[J_with_b_with_star, j_0]
        b_with_star[len(J_with_b_with_star):] = A[:, j_0]

        x_vector = (numpy.dot(numpy.dot(-1, numpy.linalg.inv(H)), b_with_star))

        counter = 0
        for i in range(len(l)):
            if i in J_with_b_with_star:
                l[i] = x_vector[counter]
                counter += 1

        print('sdg')

        break


def main():
    x_with_T = numpy.array([2, 3, 0, 0])
    J_with_b = numpy.array([0, 1])  # original is 1 2
    J_with_b_with_star = numpy.array([0, 1])  # original is 1 2
    A = numpy.array([[1, 0, 2, 1],
                     [0, 1, -1, 2]])
    D = numpy.array([[2, 1, 1, 0],
                     [1, 1, 0, 0],
                     [1, 0, 1, 0],
                     [0, 0, 0, 0]])
    c = numpy.array([-8, -6, -4, -6])
    quadratic_programming(x_with_T, J_with_b, J_with_b_with_star, A, D, c)


if __name__ == '__main__':
    main()
