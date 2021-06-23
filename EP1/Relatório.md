---
geometry:
    - left=3cm
    - top=3cm
    - bottom=2cm
    - right=2cm
numbersections: true
papersize: a4
table-of-contents: true
bibliography: biblio.bib
citation-style: abnt
monofont: JetBrainsMono-Regular.ttf
colorlinks: true
header-includes: |
    \usepackage[utf8]{inputenc}
    \usepackage{amsmath,amsthm}
    \usepackage{graphicx,cite,enumerate}
    \usepackage[brazil]{babel}
---

\begin{titlepage}
   \begin{center}
        \hspace{0pt}
        \large
        MAP3121 - Métodos Numéricos e Aplicações \\
        \textit{Escola Politécnica da Universidade de São Paulo}
        \vfill

        \Huge
        \textbf{Exercício Programa 1} \\
        \huge
        Autovalores e Autovetores de Matrizes Tridiagonais Simétricas

        \vspace{1cm}
        \Large
        Gabriel Macias de Oliveira, NUSP 11260811 \\
        Rodrigo Ryuji Ikegami, NUSP 10297265

        \vfill

        \normalsize
        São Paulo,\\        
        2021.           
   \end{center}
\end{titlepage}

\pagebreak

\tableofcontents

\pagebreak

# Introdução

## Descrição do Problema

## Ferramentas Utilizadas

## Execução dos _Scripts_

\pagebreak
# Implementação

## As Rotações de Givens

Conforme [@MAT3121], as _Rotações de Givens_ são transformações lineares ortogonais $Q:\symbb{R}^n\rightarrow\mathbb{R}^n$ da forma $Q(i, j, \theta)$ que operam sobre o espaço de matrizes como uma rotação no plano gerado pelas coordenadas $i$ e $j$. Dada uma matriz de interesse $A$, se $Y=Q(i,j,\theta)A$, então:

\begin{equation}
y_{k,l}=
    \begin{cases} 
        a_{k,l} & k\neq i,j \\
        ca_{k,l}-sa_{j,l} & k=i \\
        sa_{k,l}+ca_{j,l} & k=j
    \end{cases}\label{eq:1}
\end{equation}

sendo $c=\cos{\theta}$ e $s=\sin{\theta}$. Se $A$ é tridiagonal simétrica, podemos, particularmente, explorar essa operação para anular as entradas abaixo da diagonal principal. Se $A_n$ é definido pelos vetores $\alpha = (\alpha_1, \alpha_2, \cdots, \alpha_n)$ e $\beta = (\beta_1, \beta_2, \cdots, \beta_{n-1})$, ou seja,

$$A = 
    \begin{bmatrix}
        \alpha_1 & \beta_1 \\
        \beta_1 & \alpha_2 & \beta_2 \\
        & \ddots & \ddots & \ddots & \\
        & & \beta_{n-2} & \alpha_{n-1} & \beta_{n-1} \\
        & &  & \beta_{n-1} & \alpha_n
    \end{bmatrix}
$$

faz-se, iterativamente, $n-1$ rotações com o objetivo de transformar $A$ em uma matriz triangular superior $R$. Desta forma, para a $k$-ésima iteração, definimos $Q_k=Q(k,k+1,\theta_k)$ onde $\theta_k$ é o ângulo que permite anular a entrada $\beta_k$\ através da rotação. Podemos encontrar tal $\theta_k$ de modo numericamente estável fazendo:

$$
\begin{matrix}
    \tau_k = -\beta_k/\alpha_k &
    c_k = 1/\sqrt{1 + \tau_k^2} &
    s_k = c_k\tau_k
\end{matrix}
$$

se $\left|\alpha_k\right| >\left|\beta_k\right|$ e

$$
\begin{matrix}
    \tau_k=-\alpha_k/\beta_k \ & s_k=1/\sqrt{1+\tau_k^2} & c_k = s_k\tau_k
\end{matrix}
$$

caso contrário.

Observamos que a aplicação de sucessivas rotações geram novas entradas acima da sobrediagonal, entretanto, embora tais valores existam, não serão calculados, dado que não são necessários para a execução do algoritmo nem para encontrar autovalores e autovetores da matriz $A$.

Sumarizando o que se discutiu, as sucessivas iterações produzem $R=Q_{n-1}\cdots Q_2Q_1A$. Da inversibilidade das rotações de Givens, podemos fazer:

$$A=(Q_1^{-1}Q_2^{-1}\cdots Q_{n-1}^{-1})R$$

Da ortogonalidade da transformação, se escreve: $(Q_1^{-1}Q_2^{-1}\cdots Q_{n-1}^{-1})=(Q_1^TQ_2^T\cdots Q_{n-1}^T)=Q$. Assim,

$$A=QR$$

### Representação das Matrizes em Código

Nesta implementação, lidamos em todas as etapas com Matrizes Tridiagonais Simétricas, conforme já descrito. Uma implementação inicial pode partir de operações utilizando representações puramente matriciais, porém se trabalharmos com matrizes $n\times n$, teremos $n^2$ posições de memória ocupadas, das quais a maior parte vale $0$. Desta forma, estamos desperdiçando memória e tempo (quando consideramos que também iteramos sobre tais zeros).

Destarte, todas as matrizes nesta implementação, à exceção da matriz de autovetores, serão representadas por dois vetores, `alphas`, que armazena as entradas da diagonal principal, e `betas`, que armazena as entradas da sobrediagonal (e, portanto, da subdiagonal também, quando simétrica). Desta seção em diante, nos referenciaremos livremente a tais vetores.

### Função de Implementação da Fatoração QR

Com as considerações acima, criou-se a função `qr_factorization`. Dada uma matriz de entrada $A$, representada por seus vetores `alphas` e `betas`, retorna sua fatoração QR. 

As matrizes $Q$ e $R$ da fatoração são representadas por uma 4-tupla ordenada. As posições $0$ e $1$ da 4-tupla contêm os valores de $c_k$ e $s_k$ das rotações de Givens. Igualmente, as posições $2$ e $3$ armazenam a diagonal principal e a sobrediagonal da matriz $R$ na forma de `alphas` e `betas`. 

A implementação está no Código 1, abaixo, do qual se retiraram os comentários, mantidos no arquivo original do _script_.

\footnotesize
~~~~ {#qrfactor .python .numberLines}
def qr_factorization(alphas : np.array, betas : np.array) -> Tuple[np.array, np.array, np.array, np.array]:
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
\normalsize
**Código 1:** Função que implementa a Fatoração QR de uma matriz dada por seus vetores `alphas` e `betas`.

O código segue a descrição formal apresentada anteriormente. As linhas 6 a 13 calculam os valores de $\tau_k$, $c_k$ e $s_k$. As linhas 15 a 17 correspondem à atualização das linhas da matriz a partir dos valores da diagonal principal e sobrediagonal. Nota-se que na linha 17 utilizou-se a atribuição simultânea da linguagem para permitir a atualização dos valores de $\alpha_{k+1}$ e $\beta_k$ sem a introdução de uma variável adicional, uma vez que a atualização delas é interdependente e, se feita em sequência, não apresentaria o valor correto.

Uma vez construída a fatoração QR, implementa-se o Algoritmo QR com deslocamento espectral, o que se fará em sequência.

## A Heurística de Wilkinson

### Função `sgn`

### Função `wilkinson_h`

## O Algoritmo QR

### Função `update_matrix`

### Função `update_eigenvectors`

### Função `qr_algorithm`

\pagebreak
# Construção dos Testes

## Teste 1: Verificação do Algoritmo

### Implementação do Teste

## Teste 2: Sistema Massa-Mola com 5 Massas

### Implementação do Teste

## Teste 3: Sistema Massa-Mola com 10 Massas

### Implementação do Teste

\pagebreak
# Resultados e Discussão

\pagebreak
# Referências {-}
