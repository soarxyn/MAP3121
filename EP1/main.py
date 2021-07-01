from contextlib import suppress
import numpy as np
from typing import Tuple
from math import copysign, cos, sin, pi
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def sgn(x):
    """
        Função Sinal
        ------------

        `sgn(x)` retorna 1 se x >= 0, -1 caso contrário.
    """
    return copysign(1, x)

def qr_factorization(alphas : np.array, betas : np.array) -> Tuple[np.array, np.array, np.array, np.array]:
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
            tau_k = - betas[k] / alphas[k]
            c_ks.append(1 / np.sqrt(1 + tau_k**2))
            s_ks.append(tau_k * c_ks[k])
        else:
            tau_k = - alphas[k] / betas[k]
            s_ks.append(1 / np.sqrt(1 + tau_k**2))
            c_ks.append(tau_k * s_ks[k])

        alphas[k] = c_ks[k] * alphas[k] - s_ks[k] * betas[k]
        betas[k] *= c_ks[k - 1] if k > 0 else 1
        (alphas[k + 1], betas[k]) = (s_ks[k] * betas[k] + c_ks[k] * alphas[k + 1], c_ks[k] * betas[k] - s_ks[k] * alphas[k + 1])

    return (c_ks, s_ks, alphas, betas)

def update_matrix(c_ks : np.array, s_ks : np.array, alphas : np.array, betas : np.array) -> Tuple[np.array, np.array]:
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
        (alphas[i], betas[i], alphas[i + 1]) = (c * alphas[i] - s * betas[i], -s * alphas[i + 1], c * alphas[i + 1])

    return (alphas, betas)

def update_eigenvectors(V : np.array, c_ks : np.array, s_ks : np.array) -> np.array:
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
        (V_k[:, i], V_k[:, i + 1]) = (c * V_k[:, i] - s * V_k[:, i + 1], s * V_k[:, i] + c * V_k[:, i + 1])
    return V_k

def wilkinson_h(alphas : np.array, betas : np.array) -> float:
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

        \mu_k   :   float
            Valor do coeficiente de deslocamento espectral para a `k-ésima` iteração.
    """
    d_k = (alphas[len(alphas) - 2]  - alphas[len(alphas) - 1]) / 2
    return alphas[len(alphas) - 1] + d_k - sgn(d_k) * np.sqrt(d_k**2 + betas[len(alphas) - 2]**2)

def qr_algorithm(alphas : np.array, betas : np.array, spectralShift : bool = True, epsilon : float = 1e-6) -> Tuple[np.array, np.array, np.array, int]:
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
    V = np.identity(len(alphas_k))
    mu = 0
    iterations = 0
    for m in reversed(range(1, len(alphas))):
        while abs(betas_k[m - 1]) >= epsilon:
            (c_ks, s_ks, alphas_sub, betas_sub) = qr_factorization(alphas_k[: m + 1] - mu * np.ones(m + 1), betas_k[: m + 1])
            (alphas_k[: m + 1], betas_k[: m + 1]) = update_matrix(c_ks, s_ks, alphas_sub, betas_sub)

            alphas_k[: m + 1] += mu * np.ones(m + 1)

            V = update_eigenvectors(V, c_ks, s_ks)

            mu = wilkinson_h(alphas_k[: m + 1], betas_k[: m + 1]) if spectralShift else 0

            iterations += 1

    return (alphas_k, betas_k, V, iterations)

def qr_1(alphas : np.array, betas : np.array, shift : bool = True, eps : float = 1e-6) -> Tuple[np.array, np.array, np.array, np.array, int]:
    alphas_k = alphas.copy()
    betas_k = betas.copy()
    V = np.identity(len(alphas_k))
    mu = 0
    iterations = 0
    eigenvalues = [2 * (1 - cos(i * pi / (len(alphas_k) + 1))) for i in range(1, (len(alphas_k) + 1))][::-1]
    E_max = []
    E_avg = []
    for m in reversed(range(1, len(alphas))):
        while abs(betas_k[m - 1]) >= eps:
            (c_ks, s_ks, alphas_sub, betas_sub) = qr_factorization(alphas_k[: m + 1] - mu * np.ones(m + 1), betas_k[: m + 1])
            (alphas_k[: m + 1], betas_k[: m + 1]) = update_matrix(c_ks, s_ks, alphas_sub, betas_sub)

            alphas_k[: m + 1] += mu * np.ones(m + 1)

            V = update_eigenvectors(V, c_ks, s_ks)

            mu = wilkinson_h(alphas_k[: m + 1], betas_k[: m + 1]) if shift else 0

            E_avg.append(np.mean(np.array([abs(sorted(alphas_k, reverse = True)[i] - eigenvalues[i]) for i in range(len(alphas_k))])))
            E_max.append(max(abs(sorted(alphas_k, reverse = True)[i] - eigenvalues[i]) for i in range(len(alphas_k))))

            iterations += 1

    return (alphas_k, betas_k, V, np.array([np.array(E_avg), np.array(E_max)]), iterations)

def teste_1():
    """
        Rotinas de Testes para o Problema 1.
        Exibe autovalores e autovetores das matrizes propostas, bem como o número de iterações para convergência com e sem deslocamento espectral.
    """
    print("""
      Você selecionou o teste: Matriz com diagonal principal e subdiagonal constantes.""")

    text = ["Primeira", "Segunda", "Terceira", "Quarta"]
    for i, n in enumerate([4, 8, 16, 32]):
        print(f"""
    [=== {text[i]} Rotina: n = {n} ===]
    """)

        alphas = np.array(n * [2.0])
        betas = np.array((n - 1) * [-1.0])

        print("""      Matriz original:""")
        print("     ", np.array2string(np.diag(betas, k = 1) + np.diag(betas, k = -1) + np.diag(alphas), prefix = "      "))

        (alphas_k, _, V, iterations_sem) = qr_algorithm(alphas, betas, spectralShift = False)

        print("""\n      > Procedimentos sem deslocamento espectral <
        """)
        print(f"""      Concluído em {iterations_sem} iterações.
        """)
        print(f"""      Autovalores Encontrados: {alphas_k}\n""")

        print("""      Matriz dos Autovetores:
        """)
        print("     ", np.array2string(V, prefix = "      "))

        (alphas_k, _, V, iterations_com) = qr_algorithm(alphas, betas)

        print("""\n      > Procedimentos com deslocamento espectral <
        """)
        print(f"""      Concluído em {iterations_com} iterações. Diferença com/sem deslocamento: {iterations_sem - iterations_com} iterações.
        """)
        print(f"""      Autovalores Encontrados: {alphas_k}\n""")

        print("""      Matriz dos Autovetores:
        """)
        print("     ", np.array2string(V, prefix = "      "))

        eigenvectors = np.array([[sin(i * j * pi/ (n + 1)) for j in range(1, (n + 1))][::-1] for i in range(1, (n + 1))])

        (_, _, V, _) = qr_algorithm(alphas, betas)
        print(f"      Razão de proporcionalidade: {np.divide(eigenvectors, V)[0,0]}")

        input("\n     Pressione [ENTER] para continuar para a próxima rotina.")
        sys.stdout.write('\x1b[1A')
        sys.stdout.write('\x1b[2K')
        print("\n")

def teste_2():
    """
        Rotinas de Testes para o Problema 2.
        Exibe parâmetros das massas e molas, frequências e modos de oscilação, bem como as funções que descrevem o sistema.
    """
    print("""
      Você selecionou o teste: Sistema massa-mola com 5 molas.""")
    print("""
      =====================
      Massa das Molas: 2 kg.
      Constantes elásticas:""")

    k = [40 + 2 * i for i in range(1, 7)]

    for i, j in enumerate(k):
        print(f"        k{i + 1} = {j} N/m.")
    print("      =====================")

    alphas = np.array([(a + b)/2 for (a, b) in zip(k, k[1:])])
    betas = np.array([-b/2 for b in k[1:-1]])

    print("""
      Matriz A dos Coeficientes da EDO:
      """)

    print("     ", np.array2string(np.diag(alphas) + np.diag(betas, k = 1) + np.diag(betas, k = -1), prefix = "      "))

    (alphas_k, _, V, iterations_w) = qr_algorithm(alphas, betas)

    print("\n      > Solução da EDO <")
    print(f"""
      Número de iterações necessárias: {iterations_w}
      
      Autovalores: {alphas_k}\n
      Frequências de vibração das massas: {np.sqrt(alphas_k)}

      Modos de vibração:
      """)
    print("     ", np.array2string(V, prefix = "      "))

    for X0 in [np.array([-2.0, -3.0, -1.0, -3.0, -1.0]), np.array([1.0, 10.0, -4.0, 3.0, -2.0]), V[:, 0]]:
        print(f"\n      # Condições Iniciais: X(0) = {X0}")

        Y0 = np.matmul(np.transpose(V), X0)
        # Y(0) = Q^T x X(0)

        print(f"      Y(0) = {Y0}")

        print(f"      X(t) = [")
        print(f"\t {V[0][0] * Y0[0]:9.6f} cos({np.sqrt(alphas_k[0]):8.6f} t)", end = "")
        for j in range(1, 5):
            a = V[0][j] * Y0[j]
            if abs(a) > 1e-6:
                if a > 0:
                    print(" + ", end = "")
                elif a < 0:
                    print(" - ", end = "")

                print(f"{abs(a):8.6f} cos({np.sqrt(alphas_k[j]):8.6f} t)", end = "")

        for i in range(1, 5):
            print(",")
            print(f"\t {V[i][0] * Y0[0]:9.6f} cos({np.sqrt(alphas_k[0]):8.6f} t)", end = "")
            for j in range(1, 5):
                a = V[i][j] * Y0[j]
                if abs(a) > 1e-6:
                    if a > 0:
                        print(" + ", end = "")
                    elif a < 0:
                        print(" - ", end = "")

                    print(f"{abs(a):8.6f} cos({np.sqrt(alphas_k[j]):8.6f} t)", end = "")

        print("\n      ]\n")
        # C(t) = [cos(np.sqrt(alphas_k[i]) t) for i in range(5)]^T
        # X(t) = Q x Y(t) = Q x (Y(0) \cdot C(t))
        
        for t in [0, 5, 10]:
            print(f"      X(t = {t:2}) = {np.array([sum([V[i][j] * Y0[j] * cos(np.sqrt(alphas_k[j]) * t) for j in range(5)]) for i in range(5)])}")

def teste_3():
    """
        Rotinas de Testes para o Problema 3.
        Exibe parâmetros das massas e molas, frequências e modos de oscilação, bem como as funções que descrevem o sistema.
    """
    print("""
      Você selecionou o teste: Sistema massa-mola com 10 molas.""")
    print("""
      =====================
      Massa das Molas: 2 kg.
      Constantes elásticas:""")

    k = [40 + 2 * (-1) ** i for i in range(1, 12)]

    for i, j in enumerate(k):
            print(f"        k{i + 1} = {j} N/m.")
    print("      =====================")

    alphas = np.array([(a + b)/2 for (a, b) in zip(k, k[1:])])
    betas = np.array([-b/2 for b in k[1:-1]])

    print("""
      Matriz A dos Coeficientes da EDO:
      """)

    print("     ", np.array2string(np.diag(alphas) + np.diag(betas, k = 1) + np.diag(betas, k = -1), prefix = "      "))

    (alphas_k, _, V, iterations_w) = qr_algorithm(alphas, betas)

    print("\n      > Solução da EDO <")
    print(f"""
      Número de iterações necessárias: {iterations_w}
      
      Autovalores: {alphas_k}\n
      Frequências de vibração das massas: {np.sqrt(alphas_k)}

      Modos de vibração:
      """)
    print("     ", np.array2string(V, prefix = "      "))

    for X0 in [np.array([-2.0, -3.0, -1.0, -3.0, -1.0, -2.0, -3.0, -1.0, -3.0, -1.0]), np.array([1.0, 10.0, -4.0, 3.0, -2.0, 1.0, 10.0, -4.0, 3.0, -2.0]), V[:, 0]]:
        print(f"\n      # Condições Iniciais: X(0) = {X0}")

        Y0 = np.matmul(np.transpose(V), X0)
        # Y(0) = Q^T x X(0)

        print(f"      Y(0) = {Y0}")

        print(f"      X(t) = [")
        print(f"\t {V[0][0] * Y0[0]:9.6f} cos({np.sqrt(alphas_k[0]):8.6f} t)", end = "")
        for j in range(1, 10):
            a = V[0][j] * Y0[j]
            if abs(a) > 1e-6:
                if a > 0:
                    print(" + ", end = "")
                elif a < 0:
                    print(" - ", end = "")

                print(f"{abs(a):8.6f} cos({np.sqrt(alphas_k[j]):8.6f} t)", end = "")

        for i in range(1, 10):
            print(",")
            print(f"\t {V[i][0] * Y0[0]:9.6f} cos({np.sqrt(alphas_k[0]):8.6f} t)", end = "")
            for j in range(1, 10):
                a = V[i][j] * Y0[j]
                if abs(a) > 1e-6:
                    if a > 0:
                        print(" + ", end = "")
                    elif a < 0:
                        print(" - ", end = "")

                    print(f"{abs(a):8.6f} cos({np.sqrt(alphas_k[j]):8.6f} t)", end = "")

        print("\n      ]\n")
        # C(t) = [cos(np.sqrt(alphas_k[i]) t) for i in range(5)]^T
        # X(t) = Q x Y(t) = Q x (Y(0) \cdot C(t))
        
        for t in [0, 5, 10]:
            print(f"      X(t = {t:2}) = {np.array([sum([V[i][j] * Y0[j] * cos(np.sqrt(alphas_k[j]) * t) for j in range(10)]) for i in range(10)])}")

def plot_1():
    """
        Exibe os gráficos associados ao problema 1.
    """

    n = int(input("      Forneça o tamanho máximo da matriz (maior que 3) para calcular o número de iterações: "))
    sys.stdout.write('\x1b[1A')
    sys.stdout.write('\x1b[2K')
    print(f"      n máximo: {n}                                     \n")

    dots = 0
    print("      Preparando os gráficos")
    print(f"      Calculando o número de iterações com e sem deslocamento. (n = 3/{n})")

    iters_sem = []
    iters_com = []
    for i in range(3, n + 1):
        alphas = np.array(i * [2.0])
        betas = np.array((i-1) * [-1.0])
    
        (_, _, _, iterations) = qr_algorithm(alphas, betas, spectralShift = False)
        iters_sem.append(iterations)

        (_, _, _, iterations) = qr_algorithm(alphas, betas)
        iters_com.append(iterations)

        sys.stdout.write('\x1b[1A')
        sys.stdout.write('\x1b[1A')

        if dots == 0: 
            print("      Preparando os gráficos.")
            dots += 1
        elif dots == 1: 
            print("      Preparando os gráficos..")
            dots += 1
        elif dots == 2: 
            print("      Preparando os gráficos...")
            dots += 1
        else: 
            print("      Preparando os gráficos      ")
            dots = 0
        
        print(f"      Calculando o número de iterações com e sem deslocamento. (n = {i}/{n})")

    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    (axis1, axis2) = fig.subplots(1, 2)
    fig.suptitle("Número de Iterações por Tamanho da Matriz")

    n_space = np.linspace(3, len(iters_com), len(iters_com))

    axis1.plot(n_space, iters_sem, color = "red", lw = 1.5)
    axis1.set_title("Sem Deslocamento Espectral")
    axis1.set_xlabel("Tamanho da Matriz")
    axis1.set_ylabel("Número de Iterações")

    axis2.plot(n_space, iters_com, color = "green", lw = 1.5)
    axis2.set_title("Com Deslocamento Espectral")
    axis2.set_xlabel("Tamanho da Matriz")
    axis2.set_ylabel("Número de Iterações")

    plt.show()

    print("\n      O gráfico a seguir exibirá a evolução do erro médio dos autovalores por iteração. Calcularemos o Algoritmo QR para uma matriz tridiagonal simétrica n = 512, o que pode demorar um pouco!")

    alphas = np.array(512 * [2.0])
    betas = np.array(511 * [-1.0])

    (_, _, _, E, iterations) = qr_1(alphas, betas)
    E = E[1]

    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)

    ax = fig.add_subplot(111)
    
    ax.set_yscale('log')

    n_space = np.linspace(0, iterations, iterations)

    ax.plot(n_space, E, color = "red", lw = 1.5)
    ax.set_title("Evolução do Erro por Iteração")
    ax.set_xlabel("Iteração")
    ax.set_ylabel("Erro Médio dos Autovalores")

    plt.show()

def plot_2():
    """
        Exibe os gráficos e animações associados ao problema 2.
    """
    print("")
    k = [40 + 2 * i for i in range(1, 7)]

    alphas = np.array([(a + b)/2 for (a, b) in zip(k, k[1:])])
    betas = np.array([-b/2 for b in k[1:-1]])

    (alphas_k, _, V, _) = qr_algorithm(alphas, betas)

    initial_conditions = [np.array([-2.0, -3.0, -1.0, -3.0, -1.0]), np.array([1.0, 10.0, -4.0, 3.0, -2.0]), V[:, 0]]
    mass_positions = [[10 * i for i in range(0, 6)], [20 * i for i in range(0, 6)], [10 * i for i in range(0, 6)]]
    xlims = [50, 90, 50]
    ylims = [4, 12, 2]

    for (X0, P0, xlim, ylim) in zip(initial_conditions, mass_positions, xlims, ylims):
        print(f"      Exibindo gráfico para X(0) = {X0}")
        Y0 = np.matmul(np.transpose(V), X0)
        
        time = np.linspace(0, 10, int(10 / 0.025))
        solution = np.array([[sum([V[i][j] * Y0[j] * np.cos(np.sqrt(alphas_k[j]) * t) for j in range(5)]) for t in time] for i in range(5)])

        fig, axes = plt.subplots(2, 3)
        fig.set_size_inches(18.5, 10.5)
        fig.suptitle('Evolução do Sistema no Tempo')
        fig.subplots_adjust(hspace = 0.4)
        fig.delaxes(axes[1, 2])

        plt.figtext(0.5, .94, 'Deslocamento da Mola em relação ao Equilíbrio', fontsize = 10, ha = "center")

        axes[0, 0].plot(time, solution[0], color = "blue", lw = 1.5)
        axes[0, 1].plot(time, solution[1], color = "red", lw = 1.5)
        axes[0, 2].plot(time, solution[2], color = "lime", lw = 1.5)
        axes[1, 0].plot(time, solution[3], color = "purple", lw = 1.5)
        axes[1, 1].plot(time, solution[4], color = "coral", lw = 1.5)

        for i, axis in enumerate(axes.flat):
            axis.set_ylabel("Deslocamento (cm)")
            axis.set_xlabel("Tempo (s)")
            axis.set_title(f"Massa {i+1}")

        plt.show()  

        option = input("      Deseja visualizar uma animação para este caso? (S/n): ")

        if option != 'n' and option != 'N':
            print("")

            time = []
        
            solution = [[], [], [], [], []]

            fig, axes = plt.subplots(2, 3)
            plt.figtext(0.5, .94, 'Deslocamento da Mola em relação ao Equilíbrio', fontsize = 10, ha = "center")
            fig.set_size_inches(18.5, 10.5)
            fig.suptitle('Evolução do Sistema no Tempo')
            fig.subplots_adjust(hspace = 0.4)

            mat0 = axes[0, 0].plot([], [], color = "blue", lw = 1.5)
            mat1 = axes[0, 1].plot([], [], color = "red", lw = 1.5)
            mat2 = axes[0, 2].plot([], [], color = "lime", lw = 1.5)
            mat3 = axes[1, 0].plot([], [], color = "purple", lw = 1.5)
            mat4 = axes[1, 1].plot([], [], color = "coral", lw = 1.5)
            axes[1, 2].set_xlim(-10, xlim)
            mat5 = axes[1, 2].plot([], [], 'o')

            patches = list(mat0) + list(mat1) + list(mat2) + list(mat3) + list(mat4) + list(mat5)

            axes[0, 0].set_ylim(-ylim, ylim)
            axes[0, 1].set_ylim(-ylim, ylim)
            axes[0, 2].set_ylim(-ylim, ylim)
            axes[1, 0].set_ylim(-ylim, ylim)
            axes[1, 1].set_ylim(-ylim, ylim)

            def init():
                global time
                global solution
                time = []
                solution = [[], [], [], [], []]
                mat0[0].set_data([], [])
                mat1[0].set_data([], [])
                mat2[0].set_data([], [])
                mat3[0].set_data([], [])
                mat4[0].set_data([], [])
                mat5[0].set_data([], [])
                return patches

            def animate(index):
                t = index * 0.025
                time.append(t)

                for i in range(5):
                    solution[i].append(sum([V[i][j] * Y0[j] * np.cos(np.sqrt(alphas_k[j]) * t) for j in range(5)]))

                mat0[0].set_data(time, solution[0])
                mat1[0].set_data(time, solution[1])
                mat2[0].set_data(time, solution[2])
                mat3[0].set_data(time, solution[3])
                mat4[0].set_data(time, solution[4])
                mat0[0].axes.set_xlim(t - 1, t)
                axes[0, 1].set_xlim(t - 1, t)
                axes[0, 2].set_xlim(t - 1, t)
                axes[1, 0].set_xlim(t - 1, t)
                axes[1, 1].set_xlim(t - 1, t)

                X = np.array([P0[i] + solution[i][index] for i in range(5)])
                Y = np.zeros(np.shape(X))

                mat5[0].set_data(X,Y)
                return patches
            
            for i, axis in enumerate(axes.flat[:-1]):
                axis.set_ylabel("Deslocamento (cm)")
                axis.set_title(f"Massa {i+1}")
                axis.set_xticks([])

            anim = FuncAnimation(fig, animate, init_func=init, frames=1000000, interval=10, blit=True)

            for i, axis in enumerate(axes.flat[:-1]):
                axis.set_ylabel("Deslocamento (cm)")
            axes[1, 2].set_xlabel("Posição (cm)")
            axes[1, 2].set_yticks([])
            axes[1, 2].set_title("Sistema")

            plt.show()

def plot_3():
    """
        Exibe os gráficos e animações associados ao problema 3.
    """
    print("")
    k = [40 + 2 * i for i in range(1, 7)]

    alphas = np.array([(a + b)/2 for (a, b) in zip(k, k[1:])])
    betas = np.array([-b/2 for b in k[1:-1]])

    (alphas_k, _, V, _) = qr_algorithm(alphas, betas)

    initial_conditions = [np.array([-2.0, -3.0, -1.0, -3.0, -1.0]), np.array([1.0, 10.0, -4.0, 3.0, -2.0]), V[:, 0]]
    mass_positions = [[10 * i for i in range(0, 6)], [20 * i for i in range(0, 6)], [10 * i for i in range(0, 6)]]
    xlims = [50, 90, 50]
    ylims = [4, 12, 2]

    for (X0, P0, xlim, ylim) in zip(initial_conditions, mass_positions, xlims, ylims):
        print(f"      Exibindo gráfico para X(0) = {X0}")
        Y0 = np.matmul(np.transpose(V), X0)
        
        time = np.linspace(0, 10, int(10 / 0.025))
        solution = np.array([[sum([V[i][j] * Y0[j] * np.cos(np.sqrt(alphas_k[j]) * t) for j in range(5)]) for t in time] for i in range(5)])

        fig, axes = plt.subplots(2, 3)
        fig.set_size_inches(18.5, 10.5)
        fig.suptitle('Evolução do Sistema no Tempo')
        fig.subplots_adjust(hspace = 0.4)
        fig.delaxes(axes[1, 2])

        plt.figtext(0.5, .94, 'Deslocamento da Mola em relação ao Equilíbrio', fontsize = 10, ha = "center")

        axes[0, 0].plot(time, solution[0], color = "blue", lw = 1.5)
        axes[0, 1].plot(time, solution[1], color = "red", lw = 1.5)
        axes[0, 2].plot(time, solution[2], color = "lime", lw = 1.5)
        axes[1, 0].plot(time, solution[3], color = "purple", lw = 1.5)
        axes[1, 1].plot(time, solution[4], color = "coral", lw = 1.5)

        for i, axis in enumerate(axes.flat):
            axis.set_ylabel("Deslocamento (cm)")
            axis.set_xlabel("Tempo (s)")
            axis.set_title(f"Massa {i+1}")

        plt.show()  

        option = input("      Deseja visualizar uma animação para este caso? (S/n): ")

        if option != 'n' and option != 'N':
            print("")

            time = []
        
            solution = [[], [], [], [], []]

            fig, axes = plt.subplots(2, 3)
            plt.figtext(0.5, .94, 'Deslocamento da Mola em relação ao Equilíbrio', fontsize = 10, ha = "center")
            fig.set_size_inches(18.5, 10.5)
            fig.suptitle('Evolução do Sistema no Tempo')
            fig.subplots_adjust(hspace = 0.4)

            mat0 = axes[0, 0].plot([], [], color = "blue", lw = 1.5)
            mat1 = axes[0, 1].plot([], [], color = "red", lw = 1.5)
            mat2 = axes[0, 2].plot([], [], color = "lime", lw = 1.5)
            mat3 = axes[1, 0].plot([], [], color = "purple", lw = 1.5)
            mat4 = axes[1, 1].plot([], [], color = "coral", lw = 1.5)
            axes[1, 2].set_xlim(-10, xlim)
            mat5 = axes[1, 2].plot([], [], 'o')

            patches = list(mat0) + list(mat1) + list(mat2) + list(mat3) + list(mat4) + list(mat5)

            axes[0, 0].set_ylim(-ylim, ylim)
            axes[0, 1].set_ylim(-ylim, ylim)
            axes[0, 2].set_ylim(-ylim, ylim)
            axes[1, 0].set_ylim(-ylim, ylim)
            axes[1, 1].set_ylim(-ylim, ylim)

            def init():
                global time
                global solution
                time = []
                solution = [[], [], [], [], []]
                mat0[0].set_data([], [])
                mat1[0].set_data([], [])
                mat2[0].set_data([], [])
                mat3[0].set_data([], [])
                mat4[0].set_data([], [])
                mat5[0].set_data([], [])
                return patches

            def animate(index):
                t = index * 0.025
                time.append(t)

                for i in range(5):
                    solution[i].append(sum([V[i][j] * Y0[j] * np.cos(np.sqrt(alphas_k[j]) * t) for j in range(5)]))

                mat0[0].set_data(time, solution[0])
                mat1[0].set_data(time, solution[1])
                mat2[0].set_data(time, solution[2])
                mat3[0].set_data(time, solution[3])
                mat4[0].set_data(time, solution[4])
                mat0[0].axes.set_xlim(t - 1, t)
                axes[0, 1].set_xlim(t - 1, t)
                axes[0, 2].set_xlim(t - 1, t)
                axes[1, 0].set_xlim(t - 1, t)
                axes[1, 1].set_xlim(t - 1, t)

                X = np.array([P0[i] + solution[i][index] for i in range(5)])
                Y = np.zeros(np.shape(X))

                mat5[0].set_data(X,Y)
                return patches
            
            for i, axis in enumerate(axes.flat[:-1]):
                axis.set_ylabel("Deslocamento (cm)")
                axis.set_title(f"Massa {i+1}")
                axis.set_xticks([])

            anim = FuncAnimation(fig, animate, init_func=init, frames=1000000, interval=10, blit=True)

            for i, axis in enumerate(axes.flat[:-1]):
                axis.set_ylabel("Deslocamento (cm)")
            axes[1, 2].set_xlabel("Posição (cm)")
            axes[1, 2].set_yticks([])
            axes[1, 2].set_title("Sistema")

            plt.show()

import sys

if __name__ == "__main__":
    np.set_printoptions(precision = 12, linewidth = 200, suppress = True, sign = ' ')

    teste = int(input("""
         _____ ____  _       __  __    _    ____ _____ _ ____  _
        | ____|  _ \/ |     |  \/  |  / \  |  _ \___ // |___ \/ |
        |  _| | |_) | |_____| |\/| | / _ \ | |_) ||_ \| | __) | |
        | |___|  __/| |_____| |  | |/ ___ \|  __/___) | |/ __/| |
        |_____|_|   |_|     |_|  |_/_/   \_\_|  |____/|_|_____|_|
            [ Exercício Programa # 1 - Métodos Numéricos ]

      Autovalores e Autovetores de Matrizes Tridiagonais Simétricas
      =============================================================

      Gabriel Macias de Oliveira - NUSP: 11260811
      Rodrigo Ryuji Ikegami      - NUSP: 10297265

      Por favor, escolha uma das seguintes rotinas de teste para proseguir:

      (1) Matriz com diagonal principal e subdiagonal constantes.
      (2) Sistema massa-mola com 5 molas.
      (3) Sistema massa-mola com 10 molas.
      (4) Matriz arbitrária.
      (5) Gráficos

      Digite um número (1 - 5): """))

    if teste == 1:
        teste_1()
    elif teste == 2:
        teste_2()
    elif teste == 3:
        teste_3()
    elif teste == 4:
        print("""
      Você selecionou o teste: Matriz arbitrária.""")

        n = int(input("""
      Entre com o tamanho da matriz tridiagonal simétrica a ser diagonalizado: """))

        alphas = []
        betas = []

        print("""
      OBSERVAÇÃO: Ao inserir dados nos campos a seguir, entre com um por vez, pressionando [ENTER] entre cada entrada.""")
        print("""
      Insira as entradas da diagonal principal da matriz: """, end = '')

        alphas.append(int(input("""[""")))
        for i in range(1, n):
            sys.stdout.write('\x1b[1A')
            print("""      Insira as entradas da diagonal principal da matriz: [""", end = "")
            for elem in alphas:
                print(f"{elem}, ", end = "")
            alphas.append(int(input("")))
        sys.stdout.write('\x1b[1A')
        print(f"""      Diagonal principal: {alphas}                                      """)

        print("""
      Insira as entradas da sobrediagonal da matriz: """, end = '')

        betas.append(int(input("""[""")))
        for i in range(1, n - 1):
            sys.stdout.write('\x1b[1A')
            print("""      Insira as entradas da sobrediagonal da matriz: [""", end = "")
            for elem in betas:
                print(f"{elem}, ", end = "")
            betas.append(int(input("")))
        sys.stdout.write('\x1b[1A')
        print(f"""      Sobrediagonal: {betas}                                      """)

        spectralShift = True

        if input("""
      Utilizar deslocamento espectral? (S/n): """) == 'n':
            spectralShift = False

        (alphas_k, betas_k, V, iterations_w) = qr_algorithm(alphas, betas, spectralShift)

        print("""
      Matriz a ser diagonalizada:
        """)
        print("     ", np.array2string(np.diag(alphas) + np.diag(betas, k = -1) + np.diag(betas, k = 1), prefix = "      "))

        print(f"""
      Concluído em {iterations_w} iterações.
        """)

        print(f"""
      Autovalores: {alphas_k}
        """)

        print("""
      Autovetores:
        """)
        print("     ", np.array2string(V, prefix = "      "))

        print("\n\n")

    elif teste == 5:
        grafico = int(input("""
      Selecione o caso para o qual devem ser gerados os gráficos:

      (1) Matriz com diagonal principal e subdiagonal constantes.
      (2) Sistema massa-mola com 5 molas.
      (3) Sistema massa-mola com 10 molas.

      Digite um número (1 - 3): """))

        if grafico == 1:
            plot_1()
        elif grafico == 2:
            plot_2()
        elif grafico == 3:
            plot_3()
        else:
            print("\nInválido\n\n")
    else:
        print("\nInválido!\n\n")


    print("      Rotinas de teste concluídas! Obrigado pela execução!")