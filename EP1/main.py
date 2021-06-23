from turtle import update
import numpy as np
from typing import Tuple
from functools import partial
from math import copysign

def sgn(x):
    """
        Função Sinal
        ------------

        `sgn(x)` retorna 1 se x >= 0, -1 caso contrário.
    """
    return partial(copysign, 1)(x)

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
    d_k = (alphas[len(alphas) - 1]  - alphas[len(alphas) - 2]) / 2
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

import sys

if __name__ == "__main__":
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

      Por favor, escolha uma das seguintes rotinas de teste para proseeguir:

      (1) Matriz com diagonal principal e subdiagonal constantes.
      (2) Sistema massa-mola com 5 molas.
      (3) Sistema massa-mola com 10 molas.
      (4) Matriz arbitrária.

      Digite um número (1 - 4): """))

    np.set_printoptions(precision = 6)
    if teste == 1:
        print("""
      Você selecionou o teste: Matriz com diagonal principal e subdiagonal constantes.""")

        text = ["Primeira", "Segunda", "Terceira", "Quarta"]
        for i, n in enumerate([4, 8, 16, 32]):
            print(f"""
      [=== {text[i]} Rotina: n = {n} ===]
      """)

            alphas = np.array(n * [2.0])
            betas = np.array((n - 1) * [-1.0])

            print("""      Matriz original:
            """)
            print(np.diag(betas, k = 1) + np.diag(betas, k = -1) + np.diag(alphas))

            (alphas_k, betas_k, V, iterations_sem) = qr_algorithm(alphas, betas, spectralShift = False)

            print("""\n      > Procedimentos sem deslocamento espectral <
            """)
            print(f"""      Concluído em {iterations_sem} iterações.
            """)
            print(f"""      Autovalores Encontrados: {alphas_k}\n""")

            print("""      Matriz dos Autovetores:
            """)
            print(V)

            (alphas_k, betas_k, V, iterations_com) = qr_algorithm(alphas, betas, spectralShift = True)

            print("""\n      > Procedimentos com deslocamento espectral <
            """)
            print(f"""      Concluído em {iterations_com} iterações. Diferença com/sem deslocamento: {iterations_sem - iterations_com} iterações.
            """)
            print(f"""      Autovalores Encontrados: {alphas_k}\n""")

            print("""      Matriz dos Autovetores:
            """)
            print(V)

            input("\n     Pressione [ENTER] para continuar para a próxima rotina.")
            sys.stdout.write('\x1b[1A')
            sys.stdout.write('\x1b[2K')
            print("\n")

    elif teste == 2:
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

        print(np.diag(alphas) + np.diag(betas, k = 1) + np.diag(betas, k = -1))

        (alphas_k, betas_k, V, iterations_w) = qr_algorithm(alphas, betas)

        print(f"""
      Número de iterações necessárias: {iterations_w}

      Frequências de vibração das massas: {np.sqrt(alphas_k)}

      Modos de vibração:
      """)
        print(V)
        print("\n\n")

    elif teste == 3:
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

        print(np.diag(alphas) + np.diag(betas, k = 1) + np.diag(betas, k = -1))

        (alphas_k, betas_k, V, iterations_w) = qr_algorithm(alphas, betas)

        print(f"""
      Número de iterações necessárias: {iterations_w}

      Frequências de vibração das massas: {np.sqrt(alphas_k)}

      Modos de vibração:
      """)
        print(V)
        print("\n\n")

    elif teste == 4:
        print("""
      Você selecionou o teste: Matriz arbitrária.""")

        n = int(input("""
      Entre com o tamanho da matriz tridiagonal simétrica a ser diagonalizado: """))

        alphas = []
        betas = []

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
        print(np.diag(alphas) + np.diag(betas, k = -1) + np.diag(betas, k = 1))

        print(f"""
      Concluído em {iterations_w} iterações.
        """)

        print("""
      Autovalores:
        """)
        print(alphas_k)

        print("""
      Autovetores:
        """)
        print(V)

        print("\n\n")

    else:
        print("\nInválido!\n\n")


    print("Rotinas de teste concluídas! Obrigado pela execução!")