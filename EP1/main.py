import numpy as np 
from typing import Tuple
from functools import partial
from math import copysign

"""
    Função Sinal
    ------------

    sgn(x) retorna 1 se x > 0, -1 caso contrário.
"""
sgn = partial(copysign, 1)

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
    n : int = len(alphas)
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
    n : int = len(alphas)
    (alphas, betas) = (alphas.copy(), betas.copy())

    for k in range(n - 1):
        alphas[k] = c_ks[k] * alphas[k] - s_ks[k] * betas[k]
        betas[k] = -s_ks[k] * alphas[k + 1]

        alphas[k + 1] = c_ks[k] * alphas[k + 1]

    return (alphas, betas)