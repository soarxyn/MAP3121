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
monofont: FiraCode-Regular.ttf
colorlinks: true
header-includes: |
    \usepackage[bb=ams]{mathalpha}
    \usepackage{algorithm}
    \usepackage{algpseudocode}
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

Conforme [@MAT3121], as _Rotações de Givens_ são transformações lineares ortogonais $Q:\mathbb{R}^n\rightarrow\mathbb{R}^n$ da forma $Q(i, j, \theta)$ que operam sobre o espaço de matrizes como uma rotação no plano gerado pelas coordenadas $i$ e $j$. Dada uma matriz de interesse $A$, se $Y=Q(i,j,\theta)A$, então:

\begin{equation}
y_{k,l}=
    \begin{cases}
        a_{k,l} & k\neq i,j \\
        ca_{i,l}-sa_{j,l} & k=i \\
        sa_{i,l}+ca_{j,l} & k=j
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

A implementação está no Código \ref{code:qr_fact}, abaixo, do qual se retiraram os comentários, mantidos no arquivo original do _script_.

\scriptsize

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
**\label{code:qr_fact}Código \ref{code:qr_fact}:** Função que implementa a Fatoração QR de uma matriz dada por seus vetores `alphas` e `betas`.

O código segue a descrição formal apresentada anteriormente. As linhas 6 a 13 calculam os valores de $\tau_k$, $c_k$ e $s_k$. As linhas 15 a 17 correspondem à atualização das linhas da matriz a partir dos valores da diagonal principal e sobrediagonal. Nota-se que na linha 17 utilizou-se a atribuição simultânea da linguagem para permitir a atualização dos valores de $\alpha_{k+1}$ e $\beta_k$ sem a introdução de uma variável adicional, uma vez que a atualização delas é interdependente e, se feita em sequência, não apresentaria o valor correto.

Uma vez construída a fatoração QR, implementa-se o Algoritmo QR com deslocamento espectral, o que se fará em sequência.

## O Algoritmo QR

Em conformidade com a discussão de [@MAT3121], o _Algoritmo QR_ determina os autovalores e autovetores de uma matriz $A\in\mathbb{R}^{n\times n}$ e é construído a partir do pseudocódigo abaixo, cujos retornos são a matriz $A$ em sua forma diagonalizada $\Lambda$ e sua matriz de autovetores $V$, em que se armazena os autovetores nas respectivas colunas.

\begin{algorithm}
    \caption{Algoritmo QR} 
    \begin{algorithmic}[1]
        \State $A^{(0)}=A$
        \State $V^{(0)}=I_n$
		\For {$k=1,2,\ldots$}
            \State $A^{(k)} \rightarrow Q^{(k)} R^{(k)}$
            \State $A^{(k+1)}=R^{(k)}Q^{(k)}$
            \State $V^{(k+1)}=V^{(k)}Q^{(k)}$
		\EndFor
	\end{algorithmic} 
\end{algorithm}

Vale notar que $A^{(k+1)}=R^{(k)}Q^{(k)}=(Q^{(k)})^TQ^{(k)}R^{(k)}Q^{(k)}=(Q^{(k)})^TA^{(k)}Q^{(k)}$. Disso podemos dizer que $A^{(k+1)}$ e $A^{(k)}$ são (ortogonalmente) semelhantes. Isso significa que, se $A^{(k)}$ é tridiagonal simétrica, então $A^{(k+1)}$ também o é.

Consequentemente podemos fatorar a matriz de entrada, tridiagonal simétrica, e restaurar uma matriz, também tridiagonal simétrica, de maneira iterativa, de modo a convergir a uma matriz diagonal com os autovalores da matriz de entrada.

### Função de Atualização da Matriz

Para se implementar o Algoritmo QR, é necessária, além da fatoração QR da matriz $A$, sua reconstrução, ao fazer o produto $RQ$, de modo a criar convergência à matriz diagonal dos autovalores da matriz. Iremos, nos próximos parágrafos, explorar essa operação e como implementá-la[^1].

[^1]: **Observação:** por brevidade, nos trechos desta descrição, omitiremos os índices de iteração das operações matriciais, que ficarão subentendidos.

Recuperando a definição de $Q$, temos $Q=(Q_1^TQ_2^T\cdots Q_{n-1}^T)$, assim $RQ=R(Q_1^TQ_2^T\cdots Q_{n-1}^T)$ que podemos reescrever como $RQ=\Big\{(Q_1^TQ_2^T\cdots Q_{n-1}^T)^TR^T\Big\}^T=\Big\{Q_{n-1}\cdots Q_2Q_1R^T\Big\}^T$. Logo, a reconstrução da matriz $A$ é a aplicação das rotações de Givens às colunas de $R$.

Portanto, recuperando a definição de \ref{eq:1}, podemos analisar a operação sob mesma perspectiva, resultando, portanto nas considerações da relação \ref{eq:2} abaixo. Se $Y=AQ(i,j,\theta)^T$, vale:

\begin{equation}
y_{k,l}=
    \begin{cases}
        a_{k,l} & k\neq i,j \\
        ca_{k,i}-sa_{k,j} & l=i \\
        sa_{k,i}+ca_{k,j} & l=j
    \end{cases}\label{eq:2}
\end{equation}

Logo, computamos $RQ$ iterativamente aplicando a rotação $Q_k^T = Q^T(k, k+1, \theta_k)$ a cada passo em que $\theta_k$ é originado do par `(c_k, s_k)` de cossenos e senos definido previamente pela representação da matriz Q da fatoração de A. 

É válido notar que as parcelas em $R_{k,k+2}$ produzidas pelas rotações de Givens originais não influenciam nos cálculos em \ref{eq:2}, haja vista que as entradas da sobrediagonal $\beta'=(\beta_1', \beta_2', \cdots, \beta_{n-1}')$ são gerados pela subdiagonal, por simetria da matriz, justificando a ausência de um vetor $\gamma=(\gamma_1,\cdots,\gamma_{n-2})$ para as armazenar.

A partir dessa descrição, criou-se a função `update_matrix`. Dadas as matrizes $Q$ e $R$ oriundas da fatoração da matriz $A$, representadas pela 4-tupla ordenada, com a mesma ordem da saída de \ref{code:qr_fact}, retorna a matriz $A$, reconstruída pela aplicação das rotações reversas a $R$.

A saída é representada por uma dupla de vetores, `alphas` e `betas`, que armazenam sua diagonal principal e sua sobrediagonal, respectivamente.

A implementação está no Código \ref{code:updt_matrix}, abaixo, do qual se retiraram os comentários, mantidos no arquivo original do _script_.

\footnotesize

~~~~ {#updatematrix .python .numberLines}
def update_matrix(c_ks : np.array, s_ks : np.array, alphas : np.array, betas : np.array) -> Tuple[np.array, np.array]:
    (alphas, betas) = (alphas.copy(), betas.copy())

    for i, (c, s) in enumerate(zip(c_ks, s_ks)):
        (alphas[i], betas[i], alphas[i + 1]) = (c * alphas[i] - s * betas[i], -s * alphas[i + 1], c * alphas[i + 1])

    return (alphas, betas)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
\normalsize
**\label{code:updt_matrix}Código \ref{code:updt_matrix}:** Função que implementa a reconstrução de uma matriz fatorada, dada pelos vetores `alphas` e `betas` da matriz R e pelos vetores `c_ks` e `s_ks`, que representam os senos e cossenos que compõem a matriz Q.

O código atualiza, na linha 5, os valores das colunas da matriz a partir dos valores dos senos e cossenos das rotações de Givens. Foi utilizada atribuição simultânea das variáveis, na mesma linha, mas a atribuição sequencial dos valores, da esquerda para a direita, também funciona.

Na linha 4, utilizam-se duas funções da linguagem, `enumerate` e `zip`. A função `zip` recebe dois ou mais vetores e retorna outro vetor com os elementos dos vetores pareados em uma n-tupla, cujo comprimento é equivalente ao menor dos vetores. A função `enumerate` recebe um vetor e retorna outro vetor cujos elementos são duplas `(índice, elemento)` do vetor passado.

### Função de Atualização dos Autovetores

O mesmo desenvolvimento a respeito da operação $A^{(k+1)}=R^{(k)}Q^{(k)}$ na seção anterior vale para a operação $V^{(k+1)}=V^{(k)}Q^{(k)}$ e, por causa disso, a implementação desta operação é praticamente análoga à primeira. A matriz $Q^{(k)}=(Q_1^TQ_2^T\cdots Q_{n-1}^T)^{(k)}$ aplicada à direita atua como uma série de rotações de Givens sobre as colunas de $V^{(k)}$. Apesar disso, devemos ter em mente que a matriz de autovalores _não é tridiagonal simétrica_. Por conseguinte, devemos operar sobre a matriz $V$ pura, conforme \ref{eq:2}, sem considerar a abstração em `alphas` e `betas` como feito para as demais matrizes.

É com essa constatação que se desenvolveu a função `update_eigenvectors`, cuja implementação está no Código \ref{code:updt_eigen} abaixo, do qual omitimos os comentários que estão no _script_ original. 

\footnotesize
~~~~ {#updateeigen .python .numberLines}
def update_eigenvectors(V : np.array, c_ks : np.array, s_ks : np.array) -> np.array:
    V_k = V.copy()

    for i, (c, s) in enumerate(zip(c_ks, s_ks)):
        (V_k[:, i], V_k[:, i + 1]) = (c * V_k[:, i] - s * V_k[:, i + 1], s * V_k[:, i] + c * V_k[:, i + 1])

    return V_k
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
\normalsize
**\label{code:updt_eigen}Código \ref{code:updt_eigen}:** Função que implementa a atualização dos autovetores de $A$, armazenando-os nas colunas de $V$.

As entradas da função são as matrizes $V$ e $Q$ da `k`-ésima iteração, esta última representada por seus vetores `c_ks` e `s_ks` de cossenos e senos. Na linha 4, executamos um laço sobre as colunas de $V$, associando a cada coluna $i$ o par $(c_i, s_i)$ da rotação $Q(i, i+1, \theta_i)^T$. A rotação de Givens é executada na linha 5. A função retorna a matriz $V^{(k+1)}$, contendo os autovetores atualizados até a iteração atual do algoritmo QR.

## A Heurística de Wilkinson

### Função Sinal

### Função de Cálculo dos Coeficientes de Deslocamento

## O Algoritmo QR com Deslocamento Espectral

Para o _Algoritmo QR_, podemos acelerar sua convergência ao utilizar os Coeficientes de Deslocamento, pela Heurística de Wilkinson, para alterar os valores da diagonal principal, de tal modo que, ao efetuar uma iteração do algoritmo QR sobre ela, a nova matriz será numericamente mais próxima da matriz $\Lambda$ desejada.

[PSEUDOCÓDIGO]

Todas as deduções feitas para o _Algoritmo QR_ sem deslocamento espectral valem. Ou seja, ao introduzir o deslocamento espectral, não há grandes modificações a serem feitas à implementação do algoritmo para possibilitar o deslocamento.

### Função de Implementação do Algoritmo

Finalmente, foi criada a função `qr_algorithm`. Dada uma matriz de entrada, representada pelos vetores `alphas` e `betas`, retorna uma 4-tupla com a matriz calculada após se atingir um determinado erro `epsilon`, além da matriz $V$, com seus autovetores, e o número de iterações necessárias para se atingir o critério de parada. Além disso, a função recebe um booleano `spectralShift`, que indica se o algoritmo deve ser efetuado com ou sem deslocamento espectral.

A implementação está no Código \ref{code:qr_fact}, abaixo, do qual se retiraram os comentários, mantidos no arquivo original do _script_.

\footnotesize

~~~~ {#qralgo .python .numberLines}
def qr_algorithm(alphas : np.array, betas : np.array, spectralShift : bool = True, epsilon : float = 1e-6) -> Tuple[np.array, np.array, np.array, int]:
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
\normalsize
**\label{code:qr_algo}Código \ref{code:qr_algo}:** Função que implementa o Algoritmo QR, com ou sem deslocamento espectral, a partir de uma matriz dada por seus vetores `alphas` e `betas`, até atingir um certo erro `epsilon`.

O código segue a descrição formal apresentada anteriormente. Na linha 9 é executada a fatoração QR da submatriz cujos valores ainda não convergiram. Na linha 10, são atualizados os valores da submatriz, de acordo com a função `update_matrix`. Na linha 12 é desfeito o deslocamento espectral. Na linha 14 é atualizada a matriz dos autovalores. Na linha 16 é calculado $\mu_{k+1}$, o coeficiente da próxima iteração.

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
