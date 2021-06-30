---
geometry:
    - left=3cm
    - top=3cm
    - bottom=2cm
    - right=2cm
numbersections: true
papersize: a4
linestretch: 1.5
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

Matrizes reais simétricas surgem comumente no estudo de aplicações de métodos em engenharia. Além disso, seus autovalores e autovetores carregam informações sobre diversos modelos e descrições de muito interesse na análise e projeto de sistemas.

Neste Exercício Programa, implementamos o _Algoritmo QR_ aplicado a _Matrizes Tridiagonais Simétricas_, embora seu uso seja relevante para matrizes quaisquer. Abordaremos aspectos mais formais da implementação, como a descrição em código, bem como características de desempenho e acertividade. Por fim, o método será aplicado à solução de um Sistema de EDOs Lineares de Segunda Ordem.

Nosso objetivo é, dada uma matriz $A \in \mathbb{R}^{n\times n}$ tridiagonal simétrica, encontrar seus autovalores $\{\lambda_1, \lambda_2, \cdots, \lambda_n\}$ e seus respectivos autovetores $\{\symbf{v}_1, \symbf{v}_2, \cdots, \symbf{v}_n\}$, de forma eficiente e atendendo limites de erro e convergência. Recordamos que, segundo [@Algelin], o Teorema Espectral garante que uma matriz real simétrica é ortogonalmente diagonalizável, todos os seus autovalores são reais e podemos escolher os respectivos autovetores de modo a formar uma base ortonormal do $\mathbb{R}^n$.

Na Seção \ref{sec:impl}, detalha-se a implementação do Algoritmo QR em funções. Primeiro, há uma descrição da base matemática que orienta a construção de cada função, acompanhada pelo código e alguns comentários. A execução dos testes propostos é detalhada na Seção \ref{sec:test}, em que se retrata especificamente como foram implementados e o que se espera observar nas variáveis de retorno. Também é descrita a interface de comando (CLI) que acompanha o programa. Ao final, a Seção \ref{sec:results} contém a exibição, análise e discussão dos resultados dos testes, bem como comentários gerais sobre o desempenho do algoritmo.

## Ferramentas Utilizadas

Foram utilizadas as seguintes ferramentas para construção do código:

* Linguagem de Programação: _Python 3.7.9+_
* Bibliotecas Externas:
    * `numpy`, para trabalhar com aritmética de vetores
    * `matplotlib`, para produção de gráficos e animações
* _IDE_: _Visual Studio Code_
* Desenvolvimento Paralelo: _Git_

Além das bibliotecas externas, utilizaram-se as bibliotecas nativas `math`, para funções matemáticas básicas, `typing`, para utilizar tipos estáticos em _Python_, e `sys` para personalização da CLI.

Todos os testes em que são envolvidas métricas de tempo / número de iterações foram executados com base em um AMD Ryzen 5 3600X @ 4.2 GHz, portanto sendo suscetíveis a variações.

Todo o código está concentrado no arquivo `main.py`, cujos detalhes de execução se encontram em sequência e no arquivo `LEIA-ME.txt`.

Este relatório foi escrito em \LaTeX\.

## Execução dos _Scripts_

Estando o _Python_ atualizado para uma versão compatível, isto é, 3.7.9 ou mais recente, deve-se certificar que ambas bibliotecas `numpy` e `matplotlib` estejam instaladas. Caso contrário, basta executar `pip install -r requirements.txt` em um terminal, para recebê-las.

O arquivo principal deve ser executado no mesmo diretório em que foi descompactado, utilizando o comando `python main.py`. A exibição do terminal deve ser da CLI que acompanha o programa, conforme a Figura \ref{fig:1}.

\begin{figure}[h]
    \includegraphics[width = \linewidth]{fig_term1.png}
    \centering
    \caption{Exibição inicial da Command Line Interface (CLI) do programa.}
    \label{fig:1}
\end{figure}

Na visualização principal do programa, deve-se escolher uma entre cinco rotinas de exibição. As 3 primeiras rotinas exibem autovalores, autovetores e outras grandezas relevantes para os 3 testes propostos no enunciado em [@MAT3121]. A quarta rotina permite a diagonalização de uma matriz tridiagonal simétrica arbitrária. Por último, a quinta rotina exibe os gráficos relevantes para cada um dos testes, conforme se esclarecerá em sequência.

Ao escolher a primeira rotina, são exibidos, em sequência, os resultados para o primeiro teste considerando cada especificação do tamanho da matriz. Entre uma execução e outra, deve-se pressionar `[ENTER]` para prosseguir para a próxima entrada, como se observa na Figura \ref{fig:2} abaixo. Desejamos que o usuário possa analisar os autovalores e autovetores das matrizes propostas, bem como comparar a (_grande!_) diferença do número de iterações quando utilizamos ou não o deslocamento espectral.

\begin{figure}[h]
    \includegraphics[height = 10cm]{fig_term2.png}
    \centering
    \caption{Resultado da escolha da primeira rotina de testes.}
    \label{fig:2}
\end{figure}

\pagebreak

A segunda rotina, associada ao teste do sistema de 5 massas e 6 molas, exibe os valores dos coeficientes elásticos de cada mola, a matriz $A$ dos coeficientes do sistema de EDOs $X''(t)+AX(t)=0$ solucionado, bem como a frequência de vibração das massas e seus modos naturais de vibração.

A terceira rotina exibe as mesmas informações que a segunda, porém considerando o sistema de 10 massas e 11 molas. O que varia, de uma execução para outra, é o tamanho das matrizes exibidas, bem como a matriz a ser diagonalizada.

A quarta rotina permite ao usuário a inserção de uma matriz tridiagonal simétrica, fornecendo sua diagonal principal e sua sobrediagonal.

**Nota:** Ao fornecer as entradas, é essencial que o usuário pressione `[ENTER]` entre uma entrada e outra. Logo, o padrão de digitação deve ser, por exemplo `1 [ENTER] 2 [ENTER]` etc, para garantir que todas as entradas sejam lidas corretamente. Um exemplo de execução está na Figura \ref{fig:4} abaixo.

\begin{figure}[h]
    \includegraphics[height = 10cm]{fig_term4.png}
    \centering
    \caption{Exemplo de execução ao selecionar a quarta rotina.}
    \label{fig:4}
\end{figure}

\newpage

# Implementação {#sec:impl}

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

Em conformidade com a discussão de [@MAT3121], o _Algoritmo QR_ determina os autovalores e autovetores de uma matriz $A\in\mathbb{R}^{n\times n}$ e é construído a partir do pseudocódigo \ref{alg:qr} abaixo, cujos retornos são a matriz $A$ em sua forma diagonalizada $\Lambda$ e sua matriz de autovetores $V$, em que se armazena os autovetores nas respectivas colunas.

\begin{algorithm}
    \label{alg:qr}
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

Na linha 4, utilizam-se duas funções da linguagem, `enumerate` e `zip`. A função `zip` recebe dois ou mais vetores e retorna outro vetor com os elementos dos vetores pareados em uma n-tupla, cujo comprimento é equivalente ao menor dos vetores. A função `enumerate` recebe um vetor e retorna outro vetor cujos elementos são duplas `(índice, elemento)` do vetor de entrada.

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

## O Algoritmo QR com Deslocamento Espectral

### A Heurística de Wilkinson

De acordo com [@MAT3121], a taxa de convergência do Algoritmo QR depende da razão $|\lambda_{j+1}/\lambda_j|$, o que a torna muito lenta caso a razão entre os módulos de autovalores consecutivos esteja próxima de 1.

Visando acelerar a convergência do método, podemos subtrair da matriz $A^{(k)}$, em cada iteração, a matriz identidade múltiplicada por uma constante escalar $\mu_k$, denominada _constante de deslocamento espectral_, que esteja próxima a um autovalor.

É importante perceber que $\mu_k$ deve ser recalculado a cada iteração. Seguindo o argumento de [@MAT3121], ao alterarmos o pseudocódigo \ref{alg:qr}, fazendo iterações da forma
$$A^{(k)}-\mu_kI_n\rightarrow Q^{(k)}R^{(k)}$$
$$A^{(k+1)}=R^{(k)}Q^{(k)}+\mu_kI_n$$
a taxa de convergência de uma entrada da diagonal principal $\alpha_j^{(k)}$ para o respectivo autovalor $\lambda_j$ será proporcional a $|(\lambda_{j}-\mu_k)/(\lambda_{j-1}-\mu_k)|$. Se $\mu_k$ está próximo de $\lambda_j$, esperamos que $\beta_{j-1}^{(k)}\to 0$ mais rápido que $\beta_{i}^{(k)}, \forall i<j$, o que corresponde a $\alpha_j^{(k)}$ convergir rapidamente para $\lambda_j$.

Para calcular os coeficientes de deslocamento $\mu_k$ de cada iteração, utilizamos a heurística de Wilkinson.

**Heurística de Wilkinson:** Seja $d_k=(\alpha_{n-1}^{(k)}-\alpha_n^{(k)})/2$, definimos $\mu_k=\alpha_n^{(k)}+d_k-sgn(d_k)\sqrt{d_k^2+\big(\beta_{n-1}^{(k)} \big)^2}$, onde $sgn(d)$ é a função _sinal_, isto é:

\begin{equation}
sgn(d)=
\begin{cases}
    1,  & d \geq 0 \\
    0,  & d < 0
\end{cases}\label{eq:3}
\end{equation}

### Função Sinal

A implementação em código da função `sgn` é imediata, como se observa no código \ref{code:sgn} abaixo, do qual se retiraram os comentários, mantidos no _script_ original.

\footnotesize
~~~~ {#sgn .python .numberLines}
def sgn(x):
    return copysign(1, x)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
\normalsize
**\label{code:sgn}Código \ref{code:sgn}:** Função _sinal_ da heurística de Wilkinson.

A função `copysign(x, y)`, da biblioteca `math` retorna um número construído pelo módulo de `x` e o sinal de `y`, assim a função `sgn` recebe um $x \in \mathbb{R}$ e copia o sinal de $x$ para o número $1$, resultando no mesmo efeito descrito pela equação \ref{eq:3}.

### Função de Cálculo dos Coeficientes de Deslocamento

A função exibida no Código \ref{code:wilk} implementa o cálculo dos coeficientes de deslocamento $\mu_k$ a partir da heurística de Wilkinson.

\footnotesize
~~~~ {#wilk .python .numberLines}
def wilkinson_h(alphas : np.array, betas : np.array) -> float:
    d_k = (alphas[len(alphas) - 1]  - alphas[len(alphas) - 2]) / 2
    return alphas[len(alphas) - 1] + d_k - sgn(d_k) * np.sqrt(d_k**2 + betas[len(alphas) - 2]**2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
\normalsize
**\label{code:wilk}Código \ref{code:wilk}:** Função de cálculo dos coeficientes de deslocamento pela heurística de Wilkinson.

Todas as deduções feitas para o _Algoritmo QR_ sem deslocamento espectral valem. Isto é, todas as matrizes $A^{(k)}$ são ortogonalmente semelhantes a $A$ e, ao introduzir o deslocamento espectral, as hipóteses realizadas anteriormente para particularização do algoritmo a Matrizes Tridiagonais Simétricas ainda são verificadas.

Fixado um $\epsilon$ de tolerância, consideramos a convergência de $\alpha_n^{(k)}$ para $\lambda_n$ quando $|\beta_{n-1}^{(k)}|<\epsilon$. Já determinado o autovalor associado a `n`-ésima posição, continua-se a execução do algoritmo com a submatriz tridiagonal $n-1\times n-1$ obtida pelo _slice_ de $A$. Esta rotina se repete até obtermos todos os autovalores de $A$ mediante a tolerância fornecida.

### Função de Implementação do Algoritmo

Finalmente, foi criada a função `qr_algorithm`. Dada uma matriz de entrada, representada pelos vetores `alphas` e `betas`, retorna uma 4-tupla com a matriz calculada após se atingir um determinado erro `epsilon`, além da matriz $V$, com seus autovetores, e o número de iterações necessárias para se atingir o critério de parada. Além disso, a função recebe um booleano `spectralShift`, que indica se o algoritmo deve ser efetuado com ou sem deslocamento espectral.

A implementação está no Código \ref{code:qr_fact}, abaixo, do qual se retiraram os comentários, mantidos no arquivo original do _script_.

\tiny

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

O código segue a descrição formal apresentada anteriormente. Na linha 7 é criada uma variável `m`, que vai de $N - 1$ até $1$, utilizada para determinar o tamanho da submatriz na iteração atual. A cada iteração, a matriz sobre a qual os cálculos são feitos é menor que a anterior, de tal forma que, quando um determinado $\beta_j$ converge, seu valor não é mais calculado nas próximas iterações. A verificação da condição de convergência é feita na linha 8, ao se comparar o $\beta_j$ mais à direita, na submatriz atual, com um dado $\epsilon$. Na linha 9 é executada a fatoração QR da submatriz cujos valores ainda não convergiram. Na linha 10, são atualizados os valores da submatriz, de acordo com a função `update_matrix`. Na linha 12 é desfeito o deslocamento espectral. Na linha 14 é atualizada a matriz dos autovalores. Na linha 16 é calculado $\mu_{k+1}$, o coeficiente da próxima iteração.

\pagebreak
# Construção dos Testes {#sec:test}

## Teste 1: Verificação do Algoritmo

### Implementação do Teste

Na seção 2.3 a) de [@MAT3121], é apresentada uma família de matrizes cujos autovalores e autovetores são conhecidos. É pedido que se execute o Algoritmo QR sobre algumas dessas matrizes para testar o funcionamento da implementação feita. Além disso, é pedido para que se compare o número de iterações necessárias para a convergência

A implementação está no Código \ref{code:teste_1}, abaixo, do qual se retiraram os comentários, mantidos no arquivo original do _script_.

\footnotesize

~~~~ {#teste1 .python .numberLines}
def teste_1():
    iters_com = []
    iters_sem = []
    for n in [4, 8, 16, 32]:
        print(f"n = {n}")

        alphas = np.array(n * [2.0])
        betas = np.array((n - 1) * [-1.0])

        print("Com deslocamento espectral")
        (alphas_k, betas_k, V, E, iterations) = qr_1(alphas, betas)
        iters_com.append(iterations)

        print(f"{iterations} iterações. Autovalores: {alphas_k}\n Autovetores: \n{V}\n")
        print(f"Erro médio por iteração: {E[0]}\n Erro máximo por iteração: {E[1]}\n")

        print("Sem deslocamento espectral")
        (alphas_k, betas_k, V, E, iterations) = qr_1(alphas, betas, shift = False)
        iters_sem.append(iterations)

        print(f"{iterations} iterações. Autovalores: {alphas_k}\n Autovetores: \n{V}\n")
        print(f"Erro médio por iteração: {E[0]}\n Erro máximo por iteração: {E[1]}\n")

        eigenvalues = [2 * (1 - cos(i * pi / (n + 1))) for i in range(1, n + 1)][::-1]
        eigenvectors = np.array([[sin(i * j * pi/ (n + 1)) for j in range(1, (n + 1))][::-1] for i in range(1, (n + 1))])

        print(f"Valores esperados: {eigenvalues}. \nVetores esperados: \n{eigenvectors}\n")
        (alphas_k, betas_k, V, E, iterations) = qr_1(alphas, betas)
        print(f"Razão de proporcionalidade: \n{np.divide(eigenvectors, V)}\n")

    print(f"Iterações por n (com deslocamento): {iters_com}\n Iterações por n (sem deslocamento): {iters_sem}")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
\normalsize
**\label{code:teste_1}Código \ref{code:teste_1}:** Função que implementa o teste de verificação do Algoritmo QR, tanto com deslocamento espectral quanto sem, e imprime os dados relevantes no terminal.

O código apresentado itera sobre os casos desejados (matriz de dimensões 4x4, 8x8, 16x16 e 32x32). Nas linhas 2 e 3 são criados dois vetores, para armazenar as iterações em função do tamanho da matriz. Nas linhas 6 e 7 são criadas a diagonal principal e a sobrediagonal de acordo com o especificado em [@MAT3121], com os comprimentos adequados ao tamanho da matriz desejado. Na linha 10 é feita a execução do Algoritmo QR com deslocamento espectral e, na linha 17, é feita a execução sem deslocamento espectral. Nas linhas 11 e 18 são adicionados os valores do número de iterações para os vetores apropriados. Nas linhas 23 e 24 são calculados os valores teóricos esperados, cujas fórmulas são definidas em [@MAT3121].

A função auxiliar `qr_1` utilizada nesse código executa o algoritmo QR de mesma maneira que a função em \ref{code:qr_algo}, mas calcula os erros máximo e médio absolutos a cada iteração e os retornam em um `np.array`. Os erros foram definidos como

$E^{(k)}_{máx}=máx|\alpha^{(k)}_j-\lambda_j|$

Com $j\in[1,2,3...n]$ e

$E^{(k)}_{avg}=\frac{1}{n}\sum\limits^n_{j=1}|\alpha^{(k)}_j-\lambda_j|$

Sendo $k$ o número da iteração.

\scriptsize

~~~~ {#qr1 .python .numberLines}
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
\normalsize
**\label{code:qr_1}Código \ref{code:qr_1}:** Função auxiliar para cálculo dos erros médio e máximo absolutos por iteração do algoritmo.

## Teste 2: Sistema Massa-Mola com 5 Massas

### Implementação do Teste

A implementação está no Código \ref{code:teste_2}, abaixo, do qual se retiraram os comentários e a seção de geração dos gráficos, mantidos no arquivo original do _script_.

\scriptsize

~~~~ {#teste2 .python .numberLines}
def teste_2():
    k = [40 + 2 * i for i in range(1, 7)]

    alphas = np.array([(a + b) / 2 for (a, b) in zip(k, k[1:])])
    betas = np.array([-b/2 for b in k[1:-1]])

    print("Sem deslocamento espectral")
    (alphas_k, betas_k, V, iterations) = qr_algorithm(alphas, betas, spectralShift = False)

    print(f"{iterations} iterações. Autovalores: {alphas_k}\n Autovetores: \n{V}\n")

    print("Com deslocamento espectral")
    (alphas_k, betas_k, V, iterations) = qr_algorithm(alphas, betas)

    print(f"{iterations} iterações. Autovalores: {alphas_k}\n Autovetores: \n{V}\n")

    for X0 in enumerate([np.array([-2.0, -3.0, -1.0, -3.0, -1.0]), np.array([1.0, 10.0, -4.0, 3.0, -2.0]), V[:, 0]]):
        Y0 = np.matmul(np.transpose(V), X0)

        print(f"X(0) = {X0}\nY(0) = {Y0}")

        print(f"X(t) = [")
        print(f"\t{V[0][0] * Y0[0]:9.6f} cos({np.sqrt(alphas_k[0]):8.6f} t)", end = "")
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
            print(f"\t{V[i][0] * Y0[0]:9.6f} cos({np.sqrt(alphas_k[0]):8.6f} t)", end = "")
            for j in range(1, 5):
                a = V[i][j] * Y0[j]
                if abs(a) > 1e-6:
                    if a > 0:
                        print(" + ", end = "")
                    elif a < 0:
                        print(" - ", end = "")

                    print(f"{abs(a):8.6f} cos({np.sqrt(alphas_k[j]):8.6f} t)", end = "")

        print("\n]")

        t = 0
        print(f"X(t = {t}) = {np.array([sum([V[i][j] * Y0[j] * cos(np.sqrt(alphas_k[j]) * t) for j in range(5)]) for i in range(5)])}\n")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
\normalsize
**\label{code:teste_2}Código \ref{code:teste_2}:** Função que resolve a EDO equivalente ao sistema de 5 massas e 6 molas, para os três casos pedidos.

O código apresentado utiliza o Algoritmo QR para encontrar os autovalores e autovetores da matriz A da EDO do sistema. Após isso, imprime na tela os valores encontrados e o número de iterações necessárias para a convergência. Por fim, para cada $X(0)$ dado apresenta na tela a solução $X(t)$ encontrada para o sistema e verifica o caso $X(t = 0)$, para fins de verificação.

Na linha 2 são calculados os $k_n$ das molas e, nas linhas 4 e 5 esses valores são utilizados para construir a diagonal principal e a sobrediagonal da matriz $A$. Na linha 8 é executado o Algoritmo QR sem deslocamento espectral e, na linha 13, o algoritmo é executado com deslocamento espectral. Para o cálculo da solução do sistema, é utilizado o resultado da execução do algoritmo com deslocamento espectral.

## Teste 3: Sistema Massa-Mola com 10 Massas

### Implementação do Teste

A implementação está no Código \ref{code:teste_3}, abaixo, do qual se retiraram os comentários e a seção de geração dos gráficos, mantidos no arquivo original do _script_.

\tiny

~~~~ {#teste3 .python .numberLines}
def teste_3():
    k = [40 + 2 * (-1) ** i for i in range(1, 12)]

    alphas = np.array([(a + b) / 2 for (a, b) in zip(k, k[1:])])
    betas = np.array([-b/2 for b in k[1:-1]])

    print("Sem deslocamento espectral")
    (alphas_k, betas_k, V, iterations) = qr_algorithm(alphas, betas, spectralShift = False)

    print(f"{iterations} iterações. Autovalores: {alphas_k}\n Autovetores: \n{V}\n")

    print("Com deslocamento espectral")
    (alphas_k, betas_k, V, iterations) = qr_algorithm(alphas, betas)

    print(f"{iterations} iterações. Autovalores: {alphas_k}\n Autovetores: \n{V}\n")

    for X0 in enumerate([np.array([-2.0, -3.0, -1.0, -3.0, -1.0, -2.0, -3.0, -1.0, -3.0, -1.0]), np.array([1.0, 10.0, -4.0, 3.0, -2.0, 1.0, 10.0, -4.0, 3.0, -2.0]), V[:, 0]]):
        Y0 = np.matmul(np.transpose(V), X0)

        print(f"X(0) = {X0}\nY(0) = {Y0}")

        print(f"X(t) = [")
        print(f"\t{V[0][0] * Y0[0]:9.6f} cos({np.sqrt(alphas_k[0]):8.6f} t)", end = "")
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
            print(f"\t{V[i][0] * Y0[0]:9.6f} cos({np.sqrt(alphas_k[0]):8.6f} t)", end = "")
            for j in range(1, 10):
                a = V[i][j] * Y0[j]
                if abs(a) > 1e-6:
                    if a > 0:
                        print(" + ", end = "")
                    elif a < 0:
                        print(" - ", end = "")

                    print(f"{abs(a):8.6f} cos({np.sqrt(alphas_k[j]):8.6f} t)", end = "")

        print("\n]")

        t = 0
        print(f"X(t = {t}) = {np.array([sum([V[i][j] * Y0[j] * cos(np.sqrt(alphas_k[j]) * t) for j in range(10)]) for i in range(10)])}\n")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
\normalsize
**\label{code:teste_3}Código \ref{code:teste_3}:** Função que resolve a EDO equivalente ao sistema de 10 massas e 11 molas, para os três casos pedidos.

O código apresentado utiliza o Algoritmo QR para encontrar os autovalores e autovetores da matriz A da EDO do sistema. Após isso, imprime na tela os valores encontrados e o número de iterações necessárias para a convergência. Por fim, para cada $X(0)$ dado apresenta na tela a solução $X(t)$ encontrada para o sistema e verifica o caso $X(t = 0)$, para fins de verificação.

Sua implementação é análoga à do teste 2.

\pagebreak
# Resultados e Discussão {#sec:results}

\pagebreak
# Referências {-}
