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

def lines_to_rows(lines):
    line_to_list = map(lambda l: l.split(), lines)
    str_to_float = map(lambda l: list(map(float, l)), line_to_list)

    return enumerate(str_to_float)

def matrix_from_file(filename):
    with open(filename, encoding = 'utf-8') as file:
        matrix_size : int = int(file.readline())
        matrix = np.zeros((matrix_size, matrix_size))

        for i, line in lines_to_rows(file.readlines()):
            matrix[i,:] = line

    return matrix

def teste_1():
    print("""
      >> Você selecionou o teste A proposto pelo relatório.
      
      - Realizando leitura de arquivo: 'input-a'\n""")

    A = matrix_from_file('input-a')

    print("""      Matriz de entrada:\n""")
    print("     ", np.array2string(A, prefix = "      "))

    alphas, betas, H = tridiagonalization(A)

    print("""\n      Matriz tridiagonalizada:\n""")
    print("     ", np.array2string(np.diag(betas, k = 1) + np.diag(betas, k = -1) + np.diag(alphas), prefix = "      "))

    Lambda, _, V, _ = qr_algorithm(alphas, betas, H)

    print(f"\n      Autovalores:\n\t- Encontrados:\t{np.array(sorted(Lambda, reverse = True))}\n\t- Esperados:\t{np.array([7.0, 2.0, -1.0, -2.0])}")
    print(f"\n      Matriz de Autovetores:\n")
    print("\t    v1\t      v2\tv3\t   v4")
    print("     ", np.array2string(V, prefix = "      "))

    print("\n      OBS.: Caso λ não seja renderizado corretamente em seu terminal, este \n      char é um Lambda, indicando o autovalor da matriz A.")

    for i in range(len(Lambda)):
        result = np.matmul(A, V[:, i])

        print(f"\n      A * v{i+1} = {result}")
        print(f"      λ * v{i+1} = {Lambda[i] * V[:, i]}")

        def ratio(a, b):
            return np.array(
                [Lambda[i] if a[j] < 1e-6 else a[j] / b[j] for j in range(len(a))]
            )

        print(
            f"\n      Proporção Entrada Depois/Antes Transformação:\n\tEsperada: {Lambda[i]:.4f}\n\tObtida:   {ratio(result, V[:, i])}"
        )

        input("\n     Pressione [ENTER] para continuar.")
        sys.stdout.write('\x1b[1A')
        sys.stdout.write('\x1b[2K')
        sys.stdout.write('\x1b[1A')
        sys.stdout.write('\x1b[2K')
        print("      -")

    print(f"\n      Teste de ortogonalidade:\n")
    print("      VVt =", np.array2string(np.matmul(V, np.transpose(V)), prefix = "            "))

    print("\n      Rotina de teste concluída! Obrigado pela execução!")

def teste_2():
    print("""
      >> Você selecionou o teste B proposto pelo relatório.
      
      - Realizando leitura de arquivo: 'input-b'\n""")

    A = matrix_from_file('input-b')

    print("""      Matriz de entrada:\n""")
    print("     ", np.array2string(A, prefix = "      "))

    alphas, betas, H = tridiagonalization(A)

    n, _ = A.shape
    expectedEigenvalues = np.array( [1.0 / (2 * (1 - np.cos((2 * i + 1) * pi / (2 * n + 1)))) for i in range(n)] )
    print("""\n      Matriz tridiagonalizada:\n""")
    print("     ", np.array2string(np.diag(betas, k = 1) + np.diag(betas, k = -1) + np.diag(alphas), prefix = "      "))

    Lambda, _, V, _ = qr_algorithm(alphas, betas, H)

    print(f"\n      Autovalores:\n\t- Encontrados:\t{np.array(sorted(Lambda, reverse = True))}\n\t- Esperados:\t{expectedEigenvalues}")
    print(f"\n      Matriz de Autovetores:\n")
    print("\t    v1\t      v2\tv3\t   v4\t     v5\t  ...\t  v16\t    v17\t      v18\tv19\t  v20")
    print("     ", np.array2string(V, prefix = "      "))

    print("\n      OBS.: Caso λ não seja renderizado corretamente em seu terminal, este \n      char é um Lambda, indicando o autovalor da matriz A.")

    for i in range(len(Lambda)):
        result = np.matmul(A, V[:, i])

        print(f"\n      A * v{i+1} = {result}")
        print(f"      λ * v{i+1} = {Lambda[i] * V[:, i]}")

        def ratio(a, b):
            return np.array(
                [Lambda[i] if a[j] < 1e-6 else a[j] / b[j] for j in range(len(a))]
            )

        print(
            f"\n      Proporção Entrada Depois/Antes Transformação:\n\tEsperada: {Lambda[i]:.6f}\n\tObtida:   {ratio(result, V[:, i])}"
        )

        input("\n     Pressione [ENTER] para continuar.")
        sys.stdout.write('\x1b[1A')
        sys.stdout.write('\x1b[2K')
        sys.stdout.write('\x1b[1A')
        sys.stdout.write('\x1b[2K')
        print("      -")

    print(f"\n      Teste de ortogonalidade:\n")
    print("      VVt =", np.array2string(np.matmul(V, np.transpose(V)), prefix = "            "))

    print("\n      >> Para melhor apreciação do usuário, os resultados sem truncamento da exibição das grandezas matriciais foram\n      colocados no arquivo `result-b.txt`!")

    print("\n      Rotina de teste concluída! Obrigado pela execução!")

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

def truss_from_file(filename):
    with open(filename, encoding = 'utf-8') as file:
        total_nodes, free_nodes, _ = map(int, file.readline().split())
        p, A, E = map(float, file.readline().split())
        E *= 1e9

        treatline = lambda line: tuple( map(int, line.split()[:2]) ) + tuple( map(float, line.split()[2:]) )
        bars = list(filter(lambda line: len(line) > 0, map(treatline, file.readlines())))

        K = np.zeros((2 * free_nodes, 2 * free_nodes), float)
        M = np.zeros(free_nodes, float)

        for bar in bars:
            (i, j, theta, L) = bar
            theta = np.deg2rad(theta)
            addBar(i - 1, j - 1, L, np.cos(theta), np.sin(theta), E, p, A, M, K)
    
    return M, K, total_nodes, free_nodes, bars

def teste_3():
    print("""
      >> Você selecionou a Aplicação em Treliças Planas.
      
      - Realizando leitura de arquivo: 'input-c'\n""")

    M, K, total_nodes, free_nodes, bars = truss_from_file('input-c')

    print("""\n      Matriz de Rigidez (K):\n""")
    print("     K = ", np.array2string(K, prefix = "          "))


    print("""\n      Matriz de Massa (M):\n""")
    print("     M = ", np.array2string(np.diag(M), prefix = "          "))

    M = 1.0 / np.sqrt(M)

    for i in range(2 * len(M)):
        for j in range(2 * len(M)):
            K[i, j] *= M[i // 2] * M[j // 2]

    print("""\n      Matriz da Equação Diferencial (K~):\n""")
    print("      K~ = ", np.array2string(K, prefix = "            "))

    alphas, betas, H = tridiagonalization(K)
    Lambda, _, V, _ = qr_algorithm(alphas, betas, H)

    print("""\n      K~ Tridiagonalizado:\n""")
    print("      HK~Ht = ", np.array2string(np.diag(alphas) + np.diag(betas, -1) + np.diag(betas, 1), prefix = "               "))
    print(f"\n      Autovalores Encontrados: {np.array(sorted(Lambda))}")

    frequencies = list(map(np.sqrt, sorted(Lambda)))
    eigenvalues = sorted(Lambda)[:5]
    modes = []

    for i in map(lambda eig: list(Lambda).index(eig), eigenvalues):
        modes.append(V[:, i])

    for mode in modes:
        for i in range(len(mode)):
            mode[i] *= M[i // 2]

    print(f"\n      5 menores Frequências Encontrados: {frequencies[:5]}")
    print("""\n      Modos de vibração associados às 5 menores frequências:\n""")

    if str(input("\n      Deseja visualizar uma animação das treliças vibrando conforme cada modo de vibração? (S/n): ")).lower() != "n":
        titles = ["primeira", "segunda", "terceira", "quarta", "quinta"]
        scale = [100, 150, 100, 100, 150]

        for i in range(5):
            print(f"\n      Exibindo a {titles[i]} animação!")
            
            freq = frequencies[i]
            mode = modes[i]

            for j in range(len(mode)):
                mode[j] *= scale[i]

            print(f"\n      Frequência natural de oscilação: {freq:.4f} rad/s")
            print(f"      Modo de vibração:")
            print("      z = ", np.array2string(mode, prefix = "         "))

            print(f"\n      Iremos excitar a treliça com uma condição inicial igual a {scale[i]} vezes o modo de oscilação associado a esta frequência.")

            print(f"\n      Animação aberta. Por favor, feche a janela para prosseguir.")

            fig = plt.figure()
            fig.set_size_inches(18.5, 10.5)
            fig.suptitle(f"Evolução do Sistema no Tempo - Frequência: {freq:.4f}")

            X0 = [15, 5, 25, 15, 5, -5, 25, 15, 5, -5, 15, 5, 20, 0]
            Y0 = [40, 40, 30, 30, 30, 30, 20, 20, 20, 20, 10, 10, 0, 0]

            bar_lines = [0] * len(bars)

            ax = fig.add_subplot(111)

            mat = ax.plot(X0, Y0, "o")

            for k, (i, j, _, _) in enumerate(bars):
                bar_lines[k] = ax.plot([X0[i - 1], X0[j - 1]], [Y0[i - 1], Y0[j - 1]], c="k")

            patch = reduce(lambda a, b: a + b, bar_lines) + list(mat)

            def animate(index):
                t = index / (2 * pi * freq)
                solution = [mode[j] * np.cos(freq * t) for j in range(24)]
                solution += (0, 0, 0, 0)

                X = np.array([X0[i // 2] + solution[i] for i in range(0, 28, 2)])

                Y = np.array([Y0[i // 2] + solution[i] for i in range(1, 28, 2)])

                for k, (i, j, _, _) in enumerate(bars):
                    bar_lines[k][0].set_data([X[i - 1], X[j - 1]], [Y[i - 1], Y[j - 1]])

                mat[0].set_data(X, Y)

                return patch

            anim = FuncAnimation(fig, animate, frames=600000, interval=1, blit=True)
            plt.show()
    print("\n      Rotina de teste concluída! Obrigado pela execução!")

def teste_4():
    print("""
      >> Você selecionou o Teste com uma Matriz Arbitrária.
      
      Para este teste, há duas opções: 
        (1) Ler uma matriz de um arquivo
        (2) Entrar com uma matriz, manualmente.\n""")
    
    choice = int(input("      Escolha: "))
    sys.stdout.write('\x1b[1A')
    sys.stdout.write('\x1b[2K')
    sys.stdout.write('\x1b[1A')
    sys.stdout.write('\x1b[2K')

    if choice == 1:
        print("""
      Você escolheu o modo de leitura de uma matriz de um arquivo. Este arquivo DEVE
      ser formatado conforme as entradas para o teste A e B. Isto é, a primeira entrada
      do arquivo deve conter o tamanho da matriz e as linhas subsequentes contêm as
      entradas, linha a linha. IMPORTANTE: a última linha não pode estar em branco.\n""")

        filename = str(input("      Entre com o nome e extensão do arquivo: "))
        matrix = matrix_from_file(filename)
    
    else:
        n = int(input("""
      Você escolheu o modo de entrada manual. Entre, primeiro, com o tamanho da matriz: """))
        sys.stdout.write('\x1b[1A')
        sys.stdout.write('\x1b[2K')
        print(f"""
      Você escolheu o modo de leitura de uma matriz de um arquivo. Agora, você deve inserir as {n}
      entradas de sua matriz. Para isto, digite a entrada e pressione ENTER, até que a matriz esteja
      completa em sua exibição.\n""")

        matrix = np.zeros((n,n))
        for i in range(n):
            if i != 0:
                print('')
            for j in range(n):
                if i == 0:
                    print("      A = [", end = "")
                else:
                    print("           ", end = "")
                print("[", *matrix[i, :j], end = ' ')
                matrix[i,j] = float((input("")))
                sys.stdout.write('\x1b[1A')
            if i == 0:
                print("      A = [", end = "")
            else:
                print("           ", end = "")
            print("[", *matrix[i, :], "]", end = '')   
        print("]")     

    if choice != 2:
        print("""      Matriz de entrada:\n""")
        print("     ", np.array2string(matrix, prefix = "      "))

    alphas, betas, H = tridiagonalization(matrix)

    print("""\n      Matriz tridiagonalizada:\n""")
    print("     ", np.array2string(np.diag(betas, k = 1) + np.diag(betas, k = -1) + np.diag(alphas), prefix = "      "))

    Lambda, _, V, _ = qr_algorithm(alphas, betas, H)

    print(f"\n      Autovalores Encontrados:\t{np.array(sorted(Lambda, reverse = True))}")
    print(f"\n      Matriz de Autovetores:\n")
    print("     ", np.array2string(V, prefix = "      "))

    print("\n      OBS.: Caso λ não seja renderizado corretamente em seu terminal, este \n      char é um Lambda, indicando o autovalor da matriz A.")

    for i in range(len(Lambda)):
        result = np.matmul(matrix, V[:, i])

        print(f"\n      A * v{i+1} = {result}")
        print(f"      λ * v{i+1} = {Lambda[i] * V[:, i]}")

        def ratio(a, b):
            return np.array(
                [Lambda[i] if a[j] < 1e-6 else a[j] / b[j] for j in range(len(a))]
            )

        print(
            f"\n      Proporção Entrada Depois/Antes Transformação:\n\tEsperada: {Lambda[i]:.4f}\n\tObtida:   {ratio(result, V[:, i])}"
        )

        input("\n     Pressione [ENTER] para continuar.")
        sys.stdout.write('\x1b[1A')
        sys.stdout.write('\x1b[2K')
        sys.stdout.write('\x1b[1A')
        sys.stdout.write('\x1b[2K')
        print("      -")

    print(f"\n      Teste de ortogonalidade:\n")
    print("      VVt =", np.array2string(np.matmul(V, np.transpose(V)), prefix = "            "))

    print("\n      Rotina de teste concluída! Obrigado pela execução!")


import sys

if __name__ == "__main__":
    np.set_printoptions(precision = 6, linewidth = 250, suppress = True, sign = ' ', threshold = 10, edgeitems = 5)

    teste = int(input("""
           _____ ____ ____       __  __    _    ____ _____ _ ____  _ 
          | ____|  _ \___ \     |  \/  |  / \  |  _ \___ // |___ \/ |
          |  _| | |_) |__) |____| |\/| | / _ \ | |_) ||_ \| | __) | |
          | |___|  __// __/_____| |  | |/ ___ \|  __/___) | |/ __/| |
          |_____|_|  |_____|    |_|  |_/_/   \_\_|  |____/|_|_____|_|
                                                            
                [ Exercício Programa # 2 - Métodos  Numéricos ]

      Autovalores e Autovetores de Matrizes Reais Simétricas & Aplicações
      ===================================================================

         Gabriel Macias de Oliveira - NUSP: 11260811   - Eng. Elétrica
         Rodrigo Ryuji Ikegami      - NUSP: 10297265   - Eng. Elétrica

      Por favor, selecione um dos testes para ser executado:

      (1) Teste A.
      (2) Teste B.
      (3) Aplicação: Treliças Planas.
      (4) Matriz arbitrária.

      Digite um número (1 - 4): """))

    if teste == 1:
        teste_1()
    elif teste == 2:
        teste_2()
    elif teste == 3:
        teste_3()
    elif teste == 4:
        teste_4()
    else:
        print("\n      Inválido. Por favor, recomece, obedecendo as instruções.")    
