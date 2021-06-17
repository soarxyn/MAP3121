import numpy as np 
from typing import Tuple

def qr_factorization(alphas : np.array, betas : np.array) -> Tuple[np.array, np.array, np.array, np.array]:
    """
        Fatoração QR
        ------------
        Dada uma matriz tridiagonal simétrica, representada por dois vetores `alphas` e `betas`, que armazenam
        sua diagonal principal e sua sobrediagonal, retorna sua fatoração QR utilizando rotações de Givens.

        A matriz Q é retornada por meio de dois vetores `c_ks` e `s_ks` que recebem os cossenos e senos utilizados em cada etapa da
        fatoração. A matriz R, triangular superior, é retornada em dois vetores que armazenam sua diagonal principal e a 
        sobrediagonal. A diagonal extra que foi adicionada pela aplicação das rotações de Givens não foram calculadas pois não
        são necessárias para a implementação do algoritmo QR.

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

    (alphas, betas, lower_betas) = (alphas.copy(), betas.copy(), betas.copy())
    
    for k in range(n - 1):
        if abs(alphas[k]) > abs(lower_betas[k]):
            tau_k = - lower_betas[k] / alphas[k]
            c_k = 1 / np.sqrt(1 + tau_k**2)
            s_k = tau_k * c_k
        else:
            tau_k = - alphas[k] / lower_betas[k]
            s_k = 1 / np.sqrt(1 + tau_k**2)
            c_k = tau_k * s_k

        c_ks.append(c_k)
        s_ks.append(s_k)

        alpha_k = c_k * alphas[k] - s_k * betas[k]
        beta_k = c_k * betas[k] - s_k * alphas[k + 1]

        alpha_k1 = s_k * betas[k] + c_k * alphas[k + 1]

        if k < n - 2:
            betas[k + 1] = betas[k + 1] * c_k

        alphas[k] = alpha_k
        alphas[k + 1] = alpha_k1
        betas[k] = beta_k

    return (c_ks, s_ks, alphas, betas)

