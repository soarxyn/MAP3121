import numpy as np
from typing import Tuple
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import copysign, pi
from functools import reduce


def sgn(x):
    """
    Função Sinal
    ------------

    `sgn(x)` retorna 1 se x >= 0, -1 caso contrário.
    """
    return copysign(1, x)


def qr_factorization(
    alphas: np.array, betas: np.array
) -> Tuple[np.array, np.array, np.array, np.array]:
    """
    Fatoração QR
    ------------
    Dada uma matriz tridiagonal simétrica, representada por dois vetores `alphas` e `betas`, que armazenam
    sua diagonal principal e sua sobrediagonal, retorna sua fatoração QR utilizando rotações de Givens.

    A matriz Q é retornada por meio de dois vetores `c_ks` e `s_ks` que recebem os cossenos e senos utilizados em cada etapa da
    fatoração. A matriz R, triangular superior, é retornada em dois vetores que armazenam sua diagonal principal e a
    sobrediagonal. A diagonal extra que é adicionada pela aplicação das rotações de Givens não foi calculada pois não
    será necessárias para a implementação do algoritmo QR.

    Parâmetros
    ----------

    alphas  :   np.array
        Vetor que armazena as entradas da diagonal principal da matriz.

    betas   :   np.array
        Vetor que armazena as entradas da sobrediagonal da matriz.

    Retorna
    -------
    (c_ks, s_ks, alphas, betas) : Tuple[np.array, np.array, np.array, np.array]
        Quadra que retorna os vetores `c_ks` e `s_ks`, que são os cossenos e senos utilizados nas rotações de Givens,
        e `alphas` e `betas`, que retornam a representação da matriz R.
    """
    c_ks, s_ks = [], []
    (alphas, betas) = (alphas.copy(), betas.copy())

    for k in range(len(alphas) - 1):
        if abs(alphas[k]) > abs(betas[k]):
            tau_k = -betas[k] / alphas[k]
            c_ks.append(1 / np.sqrt(1 + tau_k ** 2))
            s_ks.append(tau_k * c_ks[k])
        else:
            tau_k = -alphas[k] / betas[k]
            s_ks.append(1 / np.sqrt(1 + tau_k ** 2))
            c_ks.append(tau_k * s_ks[k])

        alphas[k] = c_ks[k] * alphas[k] - s_ks[k] * betas[k]
        betas[k] *= c_ks[k - 1] if k > 0 else 1
        (alphas[k + 1], betas[k]) = (
            s_ks[k] * betas[k] + c_ks[k] * alphas[k + 1],
            c_ks[k] * betas[k] - s_ks[k] * alphas[k + 1],
        )

    return (np.array(c_ks), np.array(s_ks), alphas, betas)


def update_matrix(
    c_ks: np.array, s_ks: np.array, alphas: np.array, betas: np.array
) -> Tuple[np.array, np.array]:
    """
    Atualização da matriz
    ------------
    Dada uma matriz triangular superior, representada por dois vetores `alphas` e `betas`, que armazenam
    sua diagonal principal e sua sobrediagonal, e as matrizes de rotação de Givens, representadas por dois
    vetores com senos e cossenos utilizados a cada rotação, retorna uma nova matriz tridiagonal simétrica.

    Os valores fora das diagonais mencionadas, mas pertencentes ao triângulo superior da matriz não são
    utilizados para o cálculo dos valores na nova matriz.

    A matriz é retornada por meio de dois vetores, `alphas` e `betas`, que recebem sua diagonal principal
    e sua sobrediagonal, respectivamente.

    Parâmetros
    ----------

    c_ks    :   np.array
        Vetor que armazena os cossenos utilizados nas rotações de Givens.

    s_ks    :   np.array
        Vetor que armazena os senos utilizados nas rotações de Givens.

    alphas  :   np.array
        Vetor que armazena as entradas da diagonal principal da matriz.

    betas   :   np.array
        Vetor que armazena as entradas da sobrediagonal da matriz.

    Retorna
    -------
    (alphas, betas) : Tuple[np.array, np.array]
        Dupla que retorna os vetores `alphas` e `betas`, que retornam a representação da matriz
        tridiagonal simétrica, através de sua diagonal principal e sua sobrediagonal.
    """
    (alphas, betas) = (alphas.copy(), betas.copy())

    for i, (c, s) in enumerate(zip(c_ks, s_ks)):
        (alphas[i], betas[i], alphas[i + 1]) = (
            c * alphas[i] - s * betas[i],
            -s * alphas[i + 1],
            c * alphas[i + 1],
        )

    return (alphas, betas)


def update_eigenvectors(V: np.array, c_ks: np.array, s_ks: np.array) -> np.array:
    """
    Atualização dos Autovetores da Matriz
    -------------------------------------
    Dada uma matriz V, que armazena os autovetores encontrados até a `(k-1)-ésima` iteração do algoritmo QR, atualiza-os
    por meio das rotações inversas de Givens, que são construídas a partir de operações nas colunas com os cossenos e
    senos obtidos anteriormente pela fatoração QR.

    Parâmetros
    ----------

    V   :   np.array
        Matriz cujas colunas são os autovetores encontrados até a `(k-1)ésima` iteração do algoritmo QR.

    c_ks    :   np.array
        Vetor que armazena os cossenos utilizados nas rotações de Givens.

    s_ks    :   np.array
        Vetor que armazena os senos utilizados nas rotações de Givens.

    Retorna
    -------

    V_k :   np.array
        Matriz cujas colunas são os autovetores atualizados até a `k=ésima` iteração do algoritmo QR.
    """
    V_k = V.copy()
    for i, (c, s) in enumerate(zip(c_ks, s_ks)):
        (V_k[:, i], V_k[:, i + 1]) = (
            c * V_k[:, i] - s * V_k[:, i + 1],
            s * V_k[:, i] + c * V_k[:, i + 1],
        )
    return V_k


def wilkinson_h(alphas: np.array, betas: np.array) -> float:
    """
    Coeficientes de Deslocamento Espectral
    --------------------------------------
    Dada a matriz A, representada por seus valores da diagonal principal e sobrediagonal com os vetores `alphas` e `betas`,
    calcula o valor do coeficiente de deslocamento espectral da `k-ésima` iteração por meio da heurística de Wilkinson.

    Parâmetros
    ----------

    alphas  :   np.array
        Vetor da diagonal principal da matriz A.

    betas   :   np.array
        Vetor da sobrediagonal da matriz A.

    Retorna
    -------

    mu_k   :   float
        Valor do coeficiente de deslocamento espectral para a `k-ésima` iteração.
    """
    d_k = (alphas[len(alphas) - 2] - alphas[len(alphas) - 1]) / 2
    return (
        alphas[len(alphas) - 1]
        + d_k
        - sgn(d_k) * np.sqrt(d_k ** 2 + betas[len(alphas) - 2] ** 2)
    )


def qr_algorithm(
    alphas: np.array,
    betas: np.array,
    V0: np.array,
    spectralShift: bool = True,
    epsilon: float = 1e-6,
) -> Tuple[np.array, np.array, np.array, int]:
    """
    Algoritmo QR
    --------------------------------------
    Dada uma matriz tridiagonal simétrica, representada por dois vetores `alphas` e `betas`, que armazenam
    sua diagonal principal e sua sobrediagonal, efetua o Algoritmo QR, com ou sem deslocamento espectral, até
    atingir um determinado erro epsilon, e retorna a matriz com os auto-valores calculados, a matriz com os
    auto-vetores e o número total de iterações executadas pelo algoritmo

    Parâmetros
    ----------

    alphas  :   np.array
        Vetor da diagonal principal da matriz A.

    betas   :   np.array
        Vetor da sobrediagonal da matriz A.

    spectralShift : bool
        Se for True, a função executa o algoritmo com deslocamento espectral. Se for False, executa o algoritmo
        sem deslocamento espectral.

    epsilon : float
        Valor utilizado para determinar convergência dos valores calculados. Quanto menor for, menor será o erro
        do valor final calculado em relação ao ideal.

    Retorna
    -------

    (alphas, betas)   :   Tuple[np.array, np.array]
        Dupla que retorna os vetores `alphas` e `betas`, que retornam a representação da matriz
        tridiagonal simétrica dos auto-valores, através de sua diagonal principal e sua sobrediagonal.

    V : np.array
        Matriz com os auto-vetores da matriz A.

    iterations : int
        Número de iterações executadas pelo algoritmo.
    """
    alphas_k = alphas.copy()
    betas_k = betas.copy()
    V = V0.copy()
    mu = 0
    iterations = 0
    for m in reversed(range(1, len(alphas))):
        while abs(betas_k[m - 1]) >= epsilon:
            (c_ks, s_ks, alphas_sub, betas_sub) = qr_factorization(
                alphas_k[: m + 1] - mu * np.ones(m + 1), betas_k[: m + 1]
            )
            (alphas_k[: m + 1], betas_k[: m + 1]) = update_matrix(
                c_ks, s_ks, alphas_sub, betas_sub
            )

            alphas_k[: m + 1] += mu * np.ones(m + 1)

            V = update_eigenvectors(V, c_ks, s_ks)

            mu = (
                wilkinson_h(alphas_k[: m + 1], betas_k[: m + 1]) if spectralShift else 0
            )

            iterations += 1

    return (alphas_k, betas_k, V, iterations)


# def Householder(x: np.array, w: np.array) -> np.array:
#     return x - 2 * np.dot(w, x) / np.dot(w, w) * w


def tridiagonalization(A: np.array) -> Tuple[np.array, np.array, np.array]:
    A = A.copy()
    alphas = []
    betas = []

    H = np.identity(np.size(A, 0))

    for m in reversed(range(2, np.size(A, 0))):
        w_i = A[1:, 0]

        alphas.append(A[0, 0])
        betas.append(np.sqrt(np.dot(w_i, w_i)))

        w_i[0] -= betas[-1]
        w_i2 = np.dot(w_i, w_i)

        A = A[1:, 1:]

        for col in np.transpose(A):
            col -= 2 * np.dot(w_i, col) / w_i2 * w_i

        for row in A:
            row -= 2 * np.dot(w_i, row) / w_i2 * w_i

        for row in H[:, -m:]:
            row -= 2 * np.dot(w_i, row) / w_i2 * w_i

    alphas.extend(np.diag(A))
    betas.append(A[1, 0])

    return (np.array(alphas), np.array(betas), H)


def teste_1():
    A = np.array(
        [
            [2.0, 4.0, 1.0, 1.0],
            [4.0, 2.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 2.0],
            [1.0, 1.0, 2.0, 1.0],
        ]
    )
    Lambda, _, V, _ = qr_algorithm(*tridiagonalization(A))
    np.set_printoptions(precision=4, suppress=True)
    print(
        f"Autovalores:\n\tEsperados: {np.array([7.0, 2.0, -1.0, -2.0])}\n\tObtidos:   {np.array(sorted(Lambda, reverse = True))}\nAutovetores:\n{V}\n"
    )
    for i in range(4):
        result = np.matmul(A, V[:, i])

        def ratio(a, b):
            return np.array(
                [Lambda[i] if a[j] < 1e-6 else a[j] / b[j] for j in range(len(a))]
            )

        print(
            f"AV[{i}]: {result}\nProporção:\n\tEsperada: {Lambda[i]:.4f}\n\tObtida:   {ratio(result, V[:, i])}\n"
        )

    print(f"Teste de ortogonalidade: VVt =\n{np.matmul(V, np.transpose(V))}\n")


def teste_2():
    n = 20
    A = np.diag([float(n - i) for i in range(n)])
    for i in range(1, n):
        A += np.diag([float(n - j) for j in range(i, n)], i)
        A += np.diag([float(n - j) for j in range(i, n)], -i)

    Lambda, _, V, _ = qr_algorithm(*tridiagonalization(A))
    eigenvalues = np.array(
        [1.0 / (2 * (1 - np.cos((2 * i + 1) * pi / (2 * n + 1)))) for i in range(n)]
    )
    np.set_printoptions(precision=3, suppress=True, linewidth=500)
    print(
        f"Autovalores:\n\tEsperados: {eigenvalues}\n\tObtidos:   {np.array(sorted(Lambda, reverse = True))}\nAutovetores:\n{V}\n"
    )
    for i in range(n):
        result = np.matmul(A, V[:, i])

        def ratio(a, b):
            return np.array(
                [Lambda[i] if a[j] < 1e-6 else a[j] / b[j] for j in range(len(a))]
            )

        print(
            f"AV[{i}]: {result}\nProporção:\n\tEsperada: {Lambda[i]:.3f}\n\tObtida:   {ratio(result, V[:, i])}\n"
        )

    print(f"Teste de ortogonalidade: VVt =\n{np.matmul(V, np.transpose(V))}\n")


def addBar(
    i: int,
    j: int,
    L: float,
    c: float,
    s: float,
    E: float,
    p: float,
    A: float,
    M: np.array,
    K: np.array,
):
    mass_contribution = 0.5 * p * A * L
    M[i] += mass_contribution

    local_stiffness = (A * E) / L * np.array([[c ** 2, c * s], [c * s, s ** 2]])
    K[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] += local_stiffness

    if j in range(len(M)):
        M[j] += mass_contribution
        K[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] += -local_stiffness
        K[2 * j : 2 * j + 2, 2 * i : 2 * i + 2] += -local_stiffness
        K[2 * j : 2 * j + 2, 2 * j : 2 * j + 2] += local_stiffness


def aplicacao():
    K = np.zeros((24, 24), float)
    M = np.zeros(12, float)

    L = 10
    E = 200e9
    p = 7.8e3
    A = 1e-1

    bars = [
        (1, 2, 0, 10),
        (1, 4, 90, 10),
        (1, 5, 45, 14.14213562373095),
        (2, 5, 90, 10),
        (2, 4, 135, 14.14213562373095),
        (3, 4, 0, 10),
        (3, 7, 90, 10),
        (3, 8, 45, 14.14213562373095),
        (4, 5, 0, 10),
        (4, 8, 90, 10),
        (4, 9, 45, 14.14213562373095),
        (5, 6, 0, 10),
        (5, 9, 90, 10),
        (5, 8, 135, 14.14213562373095),
        (6, 9, 135, 14.14213562373095),
        (6, 10, 90, 10),
        (7, 8, 0, 10),
        (8, 9, 0, 10),
        (8, 11, 90, 10),
        (8, 12, 45, 14.14213562373095),
        (9, 10, 0, 10),
        (9, 11, 135, 14.14213562373095),
        (9, 12, 90, 10),
        (11, 12, 0, 10),
        (11, 13, 116.5650511770780, 11.18033988749895),
        (11, 14, 33.69006752597979, 18.02775637731995),
        (12, 13, 146.3099324740202, 18.02775637731995),
        (12, 14, 63.43494882292201, 11.18033988749895),
    ]

    for bar in bars:
        (i, j, theta, L) = bar
        theta = theta * pi / 180
        addBar(i - 1, j - 1, L, np.cos(theta), np.sin(theta), E, p, A, M, K)

    M = 1.0 / np.sqrt(M)

    for i in range(2 * len(M)):
        for j in range(2 * len(M)):
            K[i, j] *= M[i // 2] * M[j // 2]

    Lambda, _, V, _ = qr_algorithm(*tridiagonalization(K))

    eigenval = np.amin(Lambda)
    mode = 100 * V[:][np.argmin(Lambda)]

    for i in range(len(mode)):
        mode[i] *= M[i // 2]

    freq = np.sqrt(eigenval)

    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    fig.suptitle("Evolução do Sistema no Tempo")

    X0 = [15, 5, 25, 15, 5, -5, 25, 15, 5, -5, 15, 5, 20, 0]
    Y0 = [40, 40, 30, 30, 30, 30, 20, 20, 20, 20, 10, 10, 0, 0]

    bar_lines = [0] * len(bars)

    ax = fig.add_subplot(111)

    mat = ax.plot(X0, Y0, "o")

    for k, (i, j, _, _) in enumerate(bars):
        bar_lines[k] = ax.plot([X0[i - 1], X0[j - 1]], [Y0[i - 1], Y0[j - 1]], c="k")

    patch = reduce(lambda a, b: a + b, bar_lines) + list(mat)

    def animate(index):
        t = index * 0.025

        solution = [mode[j] * np.cos(freq * t) for j in range(24)]
        solution += (0, 0, 0, 0)

        X = np.array([X0[i // 2] + solution[i] for i in range(0, 28, 2)])
        Y = np.array([Y0[i // 2] + solution[i] for i in range(1, 28, 2)])

        for k, (i, j, _, _) in enumerate(bars):
            bar_lines[k][0].set_data([X[i - 1], X[j - 1]], [Y[i - 1], Y[j - 1]])

        mat[0].set_data(X, Y)

        return patch

    anim = FuncAnimation(fig, animate, frames=600, interval=1, blit=True)
    plt.show()


if __name__ == "__main__":
    np.set_printoptions(precision=5, suppress=True, linewidth=500)
    aplicacao()
