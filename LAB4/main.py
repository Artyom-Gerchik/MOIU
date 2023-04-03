import numpy


def dual_simplex_method(c_with_T, A, b, B):
    current_iteration = 0

    for current_iteration in range(0, 100):
        print(current_iteration + 1, ' ITERATION')

        print('\nFISRT STEP\n')

        A_B = numpy.zeros((len(A), len(A)))

        for index in range(0, len(B)):
            for jindex in range(0, len(B)):
                A_B[index][jindex] = A[index][B[jindex] - 1]

        print('A_B: \n', A_B, '\n')

        A_B_INV = numpy.linalg.inv(A_B)
        print('A_B_INV: \n', A_B_INV, '\n')

        print('\nFISRT STEP\n')

        print('\nSECOND STEP\n')

        C_B = numpy.zeros(len(B))

        for index in range(0, len(B)):
            C_B[index] = c_with_T[B[index] - 1]

        print('C_B: \n', C_B, '\n')

        print('\nSECOND STEP\n')

        print('\nTHIRD STEP\n')

        y_with_T = numpy.dot(C_B, A_B_INV)

        print('y_with_T: \n', y_with_T, '\n')

        print('\nTHIRD STEP\n')

        print('\nFOURTH STEP\n')

        K_B = numpy.dot(A_B_INV, b)

        print('K_B: \n', K_B, '\n')

        print('\nFOURTH STEP\n')

        print('\nFIFTH STEP\n')

        result = True
        j_k = 0

        for index in range(0, len(K_B)):
            if K_B[index] < 0:
                result = False
                j_k = index

        if result:
            result_return = numpy.zeros((1, len(c_with_T)))[0]
            for index in range(0, len(B)):
                result_return[B[index] - 1] = K_B[index]

            print('RESULT: \n', result_return, '\n')
            return

        print('\nFIFTH STEP\n')

        print('\nSIXTH STEP\n')

        print('j_k: \n', j_k, '\n')

        print('\nSIXTH STEP\n')

        print('\nSEVENTH STEP\n')

        delta_y = A_B_INV[j_k]

        print('delta_y: \n', delta_y, '\n')

        indexes = [i for i in range(0, len(c_with_T))]

        print('indexes: \n', indexes, '\n')

        for index in range(0, len(B)):
            indexes.remove(B[index] - 1)

        print('non-basis indexes: \n', indexes, '\n')

        mu = {}

        for index in indexes:
            mu[index] = delta_y.dot(A[:, index])

        print('mu: \n', mu, '\n')

        print('\nSEVENTH STEP\n')

        print('\nEIGHT STEP\n')

        stop = True

        for index in mu.keys():
            if mu[index] < 0:
                stop = False
                break
        if stop:
            print('Pryamaya zadacha nesovmestna')
            return

        print('\nEIGHT STEP\n')

        print('\nNINTH STEP\n')

        sigma = {}

        for index in mu.keys():
            sigma[index] = (c_with_T[index] - A[:, index].dot(y_with_T)) / mu[index]

        print('sigma: \n', sigma, '\n')

        print('\nNINTH STEP\n')

        print('\nTENTH STEP\n')

        sigma_0 = min(sigma.values())
        j_0 = 0

        for index in sigma.keys():
            if sigma[index] == sigma_0:
                j_0 = index
                break

        print('sigma_0: \n', sigma_0, '\n')
        print('j_0: \n', j_0, '\n')

        print('\nTENTH STEP\n')

        print('\nELEVENTH STEP\n')

        B[j_k] = j_0 + 1

        print('new B: \n', B, '\n')

        print('\nELEVENTH STEP\n')


def main():
    c_with_T = numpy.array([-4, -3, -7, 0, 0])  # vector of coefficients

    A = numpy.array([[-2, -1, -4, 1, 0],  # matrix of coefficients
                     [-2, -2, -2, 0, 1]])

    b = numpy.array = ([[-1],
                        [-3 / 2]])  # vector of right sides

    B = numpy.array = ([4, 5])  # start basis indexes

    dual_simplex_method(c_with_T, A, b, B)


if __name__ == '__main__':
    main();
