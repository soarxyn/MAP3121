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
        \textbf{Exercício Programa 2} \\
        \huge
        Autovalores e Autovetores de Matrizes Reais Simétricas

        \vspace{1cm}
        \Large
        Gabriel Macias de Oliveira, NUSP 11260811, Eng. Elétrica \\
        Rodrigo Ryuji Ikegami, NUSP 10297265, Eng. Elétrica

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

<!-- Matrizes reais simétricas surgem comumente no estudo de aplicações de métodos em engenharia. Além disso, seus autovalores e autovetores carregam informações sobre diversos modelos e descrições de muito interesse na análise e projeto de sistemas.

Neste Exercício Programa, implementamos o _Algoritmo QR_ aplicado a _Matrizes Tridiagonais Simétricas_, embora seu uso seja relevante para matrizes quaisquer. Abordaremos aspectos mais formais da implementação, como a descrição em código, bem como características de desempenho e acertividade. Por fim, o método será aplicado à solução de um Sistema de EDOs Lineares de Segunda Ordem.

Nosso objetivo é, dada uma matriz $A \in \mathbb{R}^{n\times n}$ tridiagonal simétrica, encontrar seus autovalores $\{\lambda_1, \lambda_2, \cdots, \lambda_n\}$ e seus respectivos autovetores $\{\symbf{v}_1, \symbf{v}_2, \cdots, \symbf{v}_n\}$, de forma eficiente e atendendo limites de erro e convergência. Recordamos que, segundo [@Algelin], o Teorema Espectral garante que uma matriz real simétrica é ortogonalmente diagonalizável, todos os seus autovalores são reais e podemos escolher os respectivos autovetores de modo a formar uma base ortonormal do $\mathbb{R}^n$.

Na Seção \ref{sec:impl}, detalha-se a implementação do Algoritmo QR em funções. Primeiro, há uma descrição da base matemática que orienta a construção de cada função, acompanhada pelo código e alguns comentários. A execução dos testes propostos é detalhada na Seção \ref{sec:test}, em que se retrata especificamente como foram implementados e o que se espera observar nas variáveis de retorno. Também é descrita a interface de comando (CLI) que acompanha o programa. Ao final, a Seção \ref{sec:results} contém a exibição, análise e discussão dos resultados dos testes, bem como comentários gerais sobre o desempenho do algoritmo. -->

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

Este relatório foi tipografado em \LaTeX\.

## Execução dos _Scripts_

Estando o _Python_ atualizado para uma versão compatível, isto é, 3.7.9 ou mais recente, deve-se certificar que ambas bibliotecas `numpy` e `matplotlib` estejam instaladas. Caso contrário, basta executar `pip install -r requirements.txt` em um terminal, para recebê-las.

O arquivo principal deve ser executado no mesmo diretório em que foi descompactado, utilizando o comando `python main.py`. A exibição do terminal deve ser da CLI que acompanha o programa, conforme a Figura \ref{fig:1}.

<!-- \begin{figure}[h]
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

\begin{figure}[H]
    \includegraphics[height = 10cm]{fig_term4.png}
    \centering
    \caption{Exemplo de execução ao selecionar a quarta rotina.}
    \label{fig:4}
\end{figure}

A quinta rotina gera gráficos relevantes aos testes propostos no enunciado em [@MAT3121]. Ao selecioná-la, será pedido um número ao usuário, correspondente ao teste para o qual se deseja gerar os gráficos.

Para o teste 1, são mostrados gráficos do número de iterações até a convergência, para a execução do Algoritmo QR sem e com deslocamento espectral, e da evolução do erro médio, todos em função de $n$, o comprimento da diagonal de $A$.

Para os testes 2 e 3, são mostrados gráficos dos deslocamentos durante os primeiros 10 segundos para todas as massas. É possível também visualizar uma animação das massas e molas, com o desenvolvimento do sistema físico ao longo do tempo.

Apenas para o teste 2, pode-se visualizar também as formas de onda dos deslocamentos avançando no tempo. A fim de evitar poluição visual, essa opção não foi inclusa para o teste 3, dada a grande quantidade de gráficos.

\begin{figure}[H]
    \includegraphics[height = 10cm]{fig_term5.png}
    \centering
    \caption{Exemplo de execução ao selecionar a quinta rotina.}
    \label{fig:5}
\end{figure} -->

\newpage

# Implementação {#sec:impl}

Para o Algoritmo QR, foi utilizada a mesma implementação do EP anterior. A única modificação feita sobre a função `qr_algorithm` foi a adição de uma parâmetro de entrada, `V0`, que é utilizado ao invés da identidade para o cálculo dos autovetores.

## As Transformações de Householder

Conforme [@MAP3121], as _Transformações de Householder_ são transformações lineares ortogonais $H_w:\mathbb{R}^n\rightarrow\mathbb{R}^n$ da forma $H_w=I-\frac{2ww^T}{w\cdot w}$ que operam sobre o espaço de vetores como uma reflexão em relação ao espaço $w^\perp$. Dado um vetor de interesse $x$, se $y=H_wx$, então:

$$
    y=x-2\frac{w\cdot x}{w\cdot w}w\,.
$$

Dados dois vetores $x$ e $y$, não nulos em $\mathbb{R}^n$, é possível definir uma transformação de Householder tal que $H_wx=\lambda y$, com $\lambda\in\mathbb{R}$. Para tanto basta se definir $w=x\pm\frac{||x||}{||y||}y$. (1)

Esta propriedade se torna extremamente útil para a tridiagonalização de matrizes reais simétricas. Para cada coluna $i$ de uma matriz dada $A$, podemos definir uma transformação de Householder com:

$$
    w_i=\tilde{a}_i+\delta\frac{||\tilde{a}_i||}{||e_{i+1}||}e_{i+1}=\tilde{a}_i+\delta||\tilde{a}_i||e_{i+1}\,.
$$

sendo $\delta=\pm1$, $e_{i+1}$ o i+1-ésimo versor da base canônica de $\mathbb{R}^n$ e $\tilde{a}_i$ composta pelos elementos da coluna $i$ de $A$, exceto os pertencentes à diagonal principal e acima dela. Isto é,

$$
    \tilde{a}_i=(0,0,...,0,A_{i+1,i},A_{i+2,i},...,A_{n-1,i},A_{n,i})^T\,.
$$

De acordo com o proposto em [@MAP3121], utilizamos $\delta$ com sinal igual ao de $A_{i+1,i}$ para cada $w_i$.

Utilizando a propriedade em (1), podemos provar que
$$
    H_{w_i}\tilde{a}_i=\tilde{a}_i-\delta w_i
    H_{w_i}\tilde{a}_i=(0,0,...,0,-\delta||\tilde{a}_i||,0,...,0)^T\,.
$$

Ou seja, a coluna $i$, após a transformação $H_{w_i}$, possui como único elemento não nulo o módulo de $\tilde{a}_i$ na posição $i$, com sinal oposto ao de $A_{i+1,i}$.

Assim, temos que

$$H_{w_1}A =
    \begin{bmatrix}
        x & x & x & ... & x \\
        x & x & x & ... & x \\
        0 & x & x & ... & x \\
        ... & ... & ... & ... & ... \\
        0 & x & x & ... & x \\
    \end{bmatrix}
$$

onde os $x$ representam valores quaisquer.

E, como $H_{w_1}$ e $A$ são simétricas, podemos fazer

$$H_{w_1}AH_{w_1} =
    \begin{bmatrix}
        x & x & 0 & ... & 0 \\
        x & x & x & ... & x \\
        0 & x & x & ... & x \\
        ... & ... & ... & ... & ... \\
        0 & x & x & ... & x \\
    \end{bmatrix}
$$

para zerar os elementos à direita da sobrediagonal na primeira coluna.

A matriz resultante, pelo mesmo motivo, também é simétrica. Assim, podemos aplicar sucessivas transformações de Householder à direita e à esquerda de $A$ de forma a obter uma matriz semelhante a $A$, $T$, tridiagonal.

Vale ressaltar que, para cada $H_{w_i}$ multiplicado à esquerda de $A$, apenas as linhas $i+1$ a $n$ têm seus elementos modificados. Simetricamente, quando se multiplica $H_{w_i}$ à direita de $A$, apenas as colunas $i+1$ a $n$ têm seus elementos modificados. (2)

Ao fim das transformações, obteremos a expressão $T=HAH^T$, com:

$$
    H^T=H_{w_1}H_{w_2}...H_{w_{n-1}}H_{w_n}\,.
$$

Para economizar tempo e memória, definimos $\bar{w}_i$ e $\bar{a}_i$, que são equivalentes a $w_i$ e $\tilde{a}_i$, respectivamente, mas sem os $i$ primeiros valores, pois são todos zero.

### Função de Implementação da Tridiagonalização de uma Matriz Real Simétrica

Feitas as considerações acima, criamos a função `tridiagonalization`, que recebe uma matriz simétrica `A` como entrada e devolve a matriz $T$ tridiagonal, representada por dois vetores, `alphas` e `betas`, que representam sua diagonal principal e sua sobrediagonal, respectivamente, e a matriz `Ht`, que representa $H^T$, descrita anteriormente.

A implementação feita se utilizada das propriedades descritas em (1) e (2) para aumentar a eficiência do código. Em cada iteração, podemos trabalhar sobre uma submatriz de $A$, sobre a qual faremos as contas, já que valores fora dela não são modificados, exceto a coluna/linha cuja maioria dos valores serão zerados. Os único valor que não são zero pertencem à diagonal principal e sua sobre/subdiagonal. Para os valores da diagonal, basta tomar $A_{i,i}$ na iteração $i$, pois seu valor não é afetado por nenhuma das transformações de Householder subsequentes. Para os valores da sobrediagonal, basta tomar $-\delta||\tilde{a}_i||$, como demonstrado anteriormente.

A implementação está no Código \ref{code:trid}, abaixo, do qual se retiraram os comentários, mantidos no arquivo original do _script_.

~~~~ {#tridiag .python .numberLines}
def tridiagonalization(A: np.array) -> Tuple[np.array, np.array, np.array]:
    A = A.copy()
    alphas = []
    betas = []

    H = np.identity(np.size(A, 0))

    for m in reversed(range(2, np.size(A, 0))):
        w_i = A[1:, 0]

        alphas.append(A[0, 0])
        betas.append(-sgn(w_i[0]) * np.sqrt(np.dot(w_i, w_i)))

        w_i[0] -= betas[-1]
        w_i2 = np.dot(w_i, w_i)

        A = A[1:, 1:]

        for Acol, Arow, Hrow in zip(np.transpose(A), A, H[:, -m:]):
            Acol -= 2 * np.dot(w_i, Acol) / w_i2 * w_i
            Arow -= 2 * np.dot(w_i, Arow) / w_i2 * w_i
            Hrow -= 2 * np.dot(w_i, Hrow) / w_i2 * w_i

    alphas.extend(np.diag(A))
    betas.append(A[1, 0])

    return (np.array(alphas), np.array(betas), H)
~~~~
**\label{code:trid}Código \ref{code:trid}:** Função que implementa a tridiagonalização de uma matriz dada, `A`.

O código segue a descrição formal apresentada anteriormente. A linha 9 define o vetor $\bar{w}_i$ da Tranformação de Householder, $H_{\bar{w}_i}$, de uma dada iteração `i` e o inicializa com $\bar{a}_i$. Sa linhas 11 e 12 adicionam os elementos calculados da diagonal principal e da sua sobrediagonal aos vetores `alphas` e `betas`, respectivamente. A linha 14 modifica o $w_i$ de acordo com a expressão $\bar{w}_i=\bar{a}_i+||\bar{a}_i||\delta e_1$. A linha 15 define uma variável auxiliar `w_i2`, equivalente a $\bar{w}_i\cdot\bar{w}_i$. A linha 17 atualiza a variável `A` para armazenar a submatriz de uma dada iteração. As linhas 19 a 22 executam as multiplicações $H_{\bar{w}_i}\bar{A}H_{\bar{w}_i}$ e $H^TH_{\bar{w}_i}$. As linhas 24 e 25 adicionam os últimos elementos da diagonal principal e da sobrediagonal da matriz resultante em `alphas` e `betas`, respectivamente.

\pagebreak

# Construção dos Testes {#sec:test}

Os aspectos matemáticos da construção de cada teste serão apresentados na Seção \ref{sec:results} junto aos resultados. Nesta seção, deseja-se detalhar a implementação _em código_ de tais testes.

<!-- ## Teste 1: Matriz Tridiagonal Simétrica a valores Constantes

### Implementação do Teste

Na seção 2.3a de [@MAT3121], é apresentada uma família de matrizes cujos autovalores e autovetores são conhecidos. É pedido que se execute o Algoritmo QR sobre algumas dessas matrizes para testar o funcionamento da implementação feita. Além disso, é pedido para que se compare o número de iterações necessárias para a convergência

\scriptsize

~~~~ {#teste1 .python .numberLines}
def teste_1():
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
~~~~
\normalsize
**\label{code:teste_1}Código \ref{code:teste_1}:** Função que implementa o teste de verificação do Algoritmo QR, tanto com deslocamento espectral quanto sem, e imprime os dados relevantes no terminal.

O código apresentado itera sobre os casos desejados (matriz de dimensões $4\times4$, $8\times8$, $16\times16$ e $32\times32$). Nas linhas 15 e 16 são criadas a diagonal principal e a sobrediagonal de acordo com o especificado em [@MAT3121], com os comprimentos adequados ao tamanho da matriz desejado. Na linha 21 é feita a execução do Algoritmo QR sem deslocamento espectral e, na linha 33, é feita a execução sem deslocamento espectral. Nas linhas 11 e 18 são adicionados os valores do número de iterações para os vetores apropriados. Na linha 45 são calculados os valores teóricos esperados dos autovetores, cujas fórmulas são definidas em [@MAT3121].

<!-- A função auxiliar `qr_1` utilizada nesse código executa o algoritmo QR de mesma maneira que a função em \ref{code:qr_algo}, calculando também os erros máximo e médio absolutos a cada iteração e os retornam em um `np.array`. Os erros máximo e médio para a `k`-ésima iteração do algoritmo foram definidos como:
$$ E_{max}^{(k)}=\max\limits_{j \in [1,n]} |\alpha^{(k)}_j-\lambda_j| $$
$$ E_{avg}^{(k)}=\frac{1}{n}\sum\limits_{j=1}^n|\alpha^{(k)}_j-\lambda_j| $$

A incorporação destes erros ao algoritmo QR, na função `qr_1` está no código \ref{code:qr_1} abaixo.

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
~~~~
\normalsize
**\label{code:qr_1}Código \ref{code:qr_1}:** Função auxiliar para cálculo dos erros médio e máximo absolutos por iteração do algoritmo. -->

<!-- ## Teste 2: Sistema Massa-Mola com 5 Massas

### Implementação do Teste

A implementação está no Código \ref{code:teste_2}, abaixo, do qual se retiraram os comentários, mantidos no arquivo original do _script_.

\scriptsize

~~~~ {#teste2 .python .numberLines}
def teste_2():
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
        
        for t in [0, 5, 10]:
            print(f"      X(t = {t:2}) = {np.array([sum([V[i][j] * Y0[j] * cos(np.sqrt(alphas_k[j]) * t) for j in range(5)]) for i in range(5)])}")
~~~~
\normalsize
**\label{code:teste_2}Código \ref{code:teste_2}:** Função que resolve a EDO equivalente ao sistema de 5 massas e 6 molas, para os três casos pedidos.

O código apresentado utiliza o Algoritmo QR para encontrar os autovalores e autovetores da matriz A da EDO do sistema. Após isso, imprime na tela os valores encontrados e o número de iterações necessárias para a convergência. Por fim, para cada $X(0)$ dado, apresenta na tela a solução $X(t)$ encontrada para o sistema e a verifica para $t=0,5,10$.

Na linha 13 são calculados os $k_n$ das molas e, nas linhas 19 e 20, esses valores são utilizados para construir a diagonal principal e a sobrediagonal da matriz $A$. Na linha 28 é executado o Algoritmo QR com deslocamento espectral.

## Teste 3: Sistema Massa-Mola com 10 Massas

### Implementação do Teste

A implementação está no Código \ref{code:teste_3}, abaixo, do qual se retiraram os comentários e a seção de geração dos gráficos, mantidos no arquivo original do _script_.

\tiny

~~~~ {#teste3 .python .numberLines}
def teste_3():
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
        
        for t in [0, 5, 10]:
            print(f"      X(t = {t:2}) = {np.array([sum([V[i][j] * Y0[j] * cos(np.sqrt(alphas_k[j]) * t) for j in range(10)]) for i in range(10)])}")
~~~~
\normalsize
**\label{code:teste_3}Código \ref{code:teste_3}:** Função que resolve a EDO equivalente ao sistema de 10 massas e 11 molas, para os três casos pedidos.

O código apresentado utiliza o Algoritmo QR para encontrar os autovalores e autovetores da matriz A da EDO do sistema. Após isso, imprime na tela os valores encontrados e o número de iterações necessárias para a convergência. Por fim, para cada $X(0)$ dado, apresenta na tela a solução $X(t)$ encontrada para o sistema e a verifica para $t=0,5,10$.

Sua implementação é análoga à do teste 2. -->

## Função Principal

A função principal do programa segue uma interface simples, com opções numeradas que o usuário pode selecionar.

\scriptsize

<!-- ~~~~ {#principal .python .numberLines}
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
~~~~
\normalsize
**\label{code:main}Código \ref{code:main}:** Função principal de execução. Implementa a CLI.

Nas linhas 4 a 26 do código, é impressa a tela principal da CLI, com o pedido da entrada de um número de 1 a 5 pelo usuário. Dependendo do valor inserido, o programa inicia uma rotina específica. As entradas 1 a 3 chamam as funções de teste descritas anteriormente. A rotina 4 é implementada na própria função principal, e seu funcionamento foi descrito na introdução. A 5ª rotina dá a opção de geração de gráficos para cada um dos 3 testes propostos em [@MAT3121]. Seu funcionamento também está descrito na introdução.

### Funções de Exibição de Gráficos

Há 4 funções criadas para produção de gráficos: `plot_1`, `plot_2`, `plot_3` e `plot_spring`. Todas foram utilizadas para produção das imagens deste relatório e são invocadas ao selecionar a rotina de exibição de gráficos da função principal. A implementação das 3 primeiras consiste basicamente na execução das funções e manipulação de dados utilizando o `matplotlib`. Além disso, são produzidas animações das simulações associadas aos problemas de molas. A função `plot_spring` é uma função auxiliar cuja implementação foi baseada na de [@Molinha]. Ela desenha, no contexto do `matplotlib`, uma mola dadas suas coordenadas, o que é bastante útil para a simulação do sistema de EDOs.

\pagebreak
# Resultados e Discussão {#sec:results}

## Resultados para o Teste 1

Conforme [@Toep], o problema de autovalores em matrizes tridiagonais simétricas cuja diagonal principal e sobrediagonal são constantes é um caso específico de aplicação da diagonalização de matrizes  Toeplitz reais de segunda ordem, que possuem uma solução analítica para seus autovalores e autovetores.

No caso do problema proposto, se a matriz de interesse $A$ é dada por um par de vetores constantes $\alpha = (2, 2, \cdots, 2)$ e $\beta=(-1, -1, \cdots, -1)$, isto é:

$$A=
    \begin{bmatrix}
        2 & -1 & & \\
        -1 & 2 & \ddots & \\
        & \ddots & \ddots & -1 \\
        & & -1 & 2
    \end{bmatrix}
$$

então os autovalores de $A$ são gerados por $\lambda_j=2\big( 1- \cos\big( \frac{j\pi}{n+1} \big) \big)$ com $j\in[1,n]$ e os respectivos autovetores são:

$$\symbf{v}_j=\bigg( \sin\bigg( \frac{j\pi}{n+1} \bigg), \sin\bigg( \frac{2j\pi}{n+1} \bigg), \cdots, \sin\bigg( \frac{nj\pi}{n+1} \bigg) \bigg) $$

Podemos utilizar tanto $\{\lambda_j\}_{j=1}^n$ quanto $\{\symbf{v}_j\}_{j=1}^n$ para avaliar se o algoritmo implementado está correto e como o erro médio dos autovalores evolui por iteração. Antes, devemos relembrar um resultado importante, retirado de [@Algelin]: 

Se $\symbf{v}_j$ é autovetor de $A$ associado ao autovalor $\lambda_j$, então seus múltiplos $k\symbf{v}_j$ também são autovetores de $A$ associados ao mesmo autovalor.

### Autovalores e Autovetores para $n = 4$

O _Algoritmo QR_ foi executado, com e sem deslocamento espectral, para uma matriz $A \in \mathbb{R}^{4\times4}$ seguindo a descrição acima com precisão $\epsilon = 10^{-6}$. Os autovalores esperados, isto é, calculados analiticamente, e obtidos pela execução do algoritmo estão dispostos na Tabela \ref{table:1} abaixo.

\setlength{\tabcolsep}{18pt}

\begin{table}[h!]
    \centering
    \begin{tabular}{|c|c|c|}
        \hline
        Autovalor & Esperado & Obtido \\
        \hline
        $\lambda_1$ & $3,618033988749895$ & $3,618033988749896$ \\
        $\lambda_2$ & $2,618033988749895$ & $2,618033988749895$ \\
        $\lambda_3$ & $1,381966011250105$ & $1,381966011250106$ \\
        $\lambda_4$ & $0,381966011250105$ & $0,381966011250105$ \\
        \hline
    \end{tabular}
    \caption{Autovalores esperados e obtidos após a execução do Algoritmo, para $n=4$.}
    \label{table:1}
\end{table}

Os autovetores da matriz são apresentados nas matrizes abaixo. A primeira matriz contém os autovetores obtidos pela execução do algoritmo. Já a segunda, os autovetores calculados analiticamente. É conveniente lembrar que cada autovetor $\symbf{v}_j$ é armazenado nas colunas de $V$.

$$
V =
\begin{bmatrix}
     0.371748 &  0.601501 &  0.601501 & 0.371748 \\
    -0.601501 & -0.371748 &  0.371748 & 0.601501 \\
     0.601501 & -0.371748 & -0.371748 & 0.601501 \\
    -0.371748 &  0.601501 & -0.601501 & 0.371748 \\
\end{bmatrix}
$$

\begin{center} Matriz dos Autovetores obtidos pelo Algoritmo, para $n=4$. \end{center}

$$
V =
\begin{bmatrix}
     0.587785 &  0.951057 &  0.951057 & 0.587785 \\
    -0.951057 & -0.587785 &  0.587785 & 0.951057 \\
     0.951057 & -0.587785 & -0.587785 & 0.951057 \\
    -0.587785 &  0.951057 & -0.951057 & 0.587785 \\
\end{bmatrix}
$$

\begin{center} Matriz dos Autovetores calculados analiticamente, para $n=4$. \end{center}

### Autovalores e Autovetores para $n = 8$

Assim como feito para $n = 4$, executou-se o _Algoritmo QR_ para $A\in\mathbb{R}^{8\times8}$, tridiagonal simétrica de diagonal e subdiagonal constantes com $\alpha_k = 2$ e $\beta_k = -1$, com e sem deslocamento espectral e com precisão $\epsilon=10^{-6}$. Apresentamos os autovalores calculados analiticamente e aqueles obtidos após a execução do método na Tabela \ref{table:2} seguinte.

\begin{table}[h!]
    \centering
    \begin{tabular}{|c|c|c|}
        \hline
        Autovalor & Esperado & Obtido \\
        \hline
        $\lambda_1$ & $3,879385241571817$ & $3,879385241571820$ \\
        $\lambda_2$ & $3,532088886237956$ & $3,532088886237959$ \\
        $\lambda_3$ & $2,9999999999999996$ & $2,999999999999889$ \\
        $\lambda_4$ & $2,3472963553338606$ & $2,347296355333974$ \\
        $\lambda_5$ & $1,6527036446661392$ & $1,652703644666068$ \\
        $\lambda_6$ & $0,9999999999999998$ & $1,000000000000071$ \\
        $\lambda_7$ & $0,4679111137620440$ & $0,467911113762044$ \\
        $\lambda_8$ & $0,12061475842818314$ & $0,120614758428183$ \\
        \hline
    \end{tabular}
    \caption{Autovalores esperados e obtidos após a execução do Algoritmo, para $n=8$.}
    \label{table:2}
\end{table}

As matrizes a seguir exibem os autovetores resultantes da aplicação do algoritmo, assim como os autovetores dados pela forma analítica, respectivamente.$$
V =
\begin{bmatrix}
     0,16123   & 0,303013 &  0,408248 & 0,464243& 0,464243 & 0,408248 & 0,303013 &0,16123 \\
     0,303013  &-0,464243 & -0,408248 &-0,16123 & 0,16123  & 0,408248 & 0,464243 &0,303013\\
     0,408248  & 0,408248 &  0,000000 &-0,408248&-0,408248 & 0,000000 & 0,408248 &0,408248\\
     0,464243  &-0,16123  &  0,408248 & 0,303013&-0,303013 &-0,408248 & 0,16123  &0,464243\\
     0,464243  &-0,16123  & -0,408248 & 0,303013& 0,303013 &-0,408248 &-0,16123  &0,464243\\
     0,408248  & 0,408248 &  0,000000 &-0,408248& 0,408248 & 0,000000 &-0,408248 &0,408248\\
     0,303013  &-0,464243 &  0,408248 &-0,16123 &-0,16123  & 0,408248 &-0,464243 &0,303013\\
     0,16123   & 0,303013 & -0,408248 & 0,464243&-0,464243 & 0,408248 &-0,303013 &0,16123 \\
\end{bmatrix}
$$

\begin{center} Matriz dos Autovetores obtidos pelo Algoritmo, para $n=8$ \end{center}

$$
V =
\begin{bmatrix}
     0,34202 &  0,642788&  0,866025&  0,984808&  0,984808&  0,866025&  0,642788&  0,34202  \\
    -0,642788& -0,984808& -0,866025& -0,34202 &  0,34202 &  0,866025&  0,984808&  0,642788 \\
     0,866025&  0,866025&  0,000000& -0,866025& -0,866025&  0,000000&  0,866025&  0,866025 \\
    -0,984808& -0,34202 &  0,866025&  0,642788& -0,642788& -0,866025&  0,34202 &  0,984808 \\
     0,984808& -0,34202 & -0,866025&  0,642788&  0,642788& -0,866025& -0,34202 &  0,984808 \\
    -0,866025&  0,866025&  0,000000& -0,866025&  0,866025&  0,000000& -0,866025&  0,866025 \\
     0,642788& -0,984808&  0,866025& -0,34202 & -0,34202 &  0,866025& -0,984808&  0,642788 \\
    -0,34202 &  0,642788& -0,866025&  0,984808& -0,984808&  0,866025& -0,642788&  0,34202 \\
\end{bmatrix}
$$

\begin{center} Matriz dos Autovetores calculados analiticamente, para $n=8$ \end{center}

### Discussão sobre os Autovalores e Autovetores Encontrados

Inicialmente, podemos comparar os autovalores obtidos. Observamos que tanto para o caso em que $n = 4$ como para $n = 8$, todos os autovalores convergiram e são muito próximos dos valores corretos. O erro médio absoluto entre o autovalor obtido e o autovalor calculado analiticamente é da ordem de $e=10^{-12}$. Julgamos que nossa implementação está concordante com o proposto e a teoria.

O Algoritmo QR é **muito** eficiente em calcular os autovalores, principalmente ao considerarmos que conseguimos atingir erros extremamente pequenos com pouquíssimas iterações: 9 para o caso $4\times4$ e 19 para o $8\times8$. Além disso, chegamos a este resultando com precisão $\epsilon=10^{-6}$ para a convergência dos $\beta_j$. O Algoritmo converge muito rapidamente e com alta precisão.

Ao analisarmos os autovetores, exibidos nas matrizes, pode-se notar que a matriz de autovetores encontrada **não é igual** àquela obtida pela forma de geração. Contudo, isso não deve ser tratado como um erro, devido à propriedade de proporcionalidade, mencionada anteriormente. 

A matriz obtida é, aproximadamente, um múltiplo da matriz calculada analiticamente. Para comparar essas grandezas, definimos uma _matriz de proporcionalidade_ $P_n\in\mathbb{R}^{n\times n}$, de modo que dada a matriz aproximada de autovetores $\tilde{V}$ e a matriz correta de autovetores $V$, definimos:
$$
\begin{matrix}
    p_{i,j}=\frac{v_{i,j}}{\tilde{v}_{i,j}} &i,j=1\dots n
\end{matrix}
$$

Tanto em $P_4$ como em $P_8$, todas as entradas da matriz são iguais. Apresentamos, em seguida, os valores dessas matrizes.

\begin{align*}
P_4 &= \begin{bmatrix}
    1,58113883 & \cdots & 1,58113883 \\
    \vdots & \ddots & \vdots \\
    1,58113883 & \cdots & 1,58113883
\end{bmatrix}
&
P_8 &= \begin{bmatrix}
    2,12132034 & \cdots & 2,12132034 \\
    \vdots & \ddots & \vdots \\
    2,12132034 & \cdots & 2,12132034
\end{bmatrix}
\end{align*}

Portanto, concluímos que os autovetores obtidos pela execução do Algoritmo QR, $\{\symbf{\tilde{v}}_j\}_{j=1}^n$, são múltiplos daqueles encontrados pela forma fechada. Geram, portanto, o mesmo espaço, são autovetores da matriz de entrada, bem como estão associados aos mesmos autovalores. O que resta discutir é a constante de proporcionalidade. Para entendê-la, devemos nos lembrar que o Algoritmo QR produz uma **base ortonormal** de autovetores, mas a fórmula dos autovetores **não** é _normalizada_. Na sequência, iremos encontrar a norma desses vetores.[^2]. Sendo: $$\symbf{v}_j=\bigg( \sin\bigg( \frac{j\pi}{n+1} \bigg), \sin\bigg( \frac{2j\pi}{n+1} \bigg), \cdots, \sin\bigg( \frac{nj\pi}{n+1} \bigg) \bigg)$$
temos: $$\lVert \symbf{v}_j \rVert=\sqrt{\sin^2\bigg( \frac{j\pi}{n+1} \bigg)+\sin^2\bigg( \frac{2j\pi}{n+1} \bigg)+\cdots+\sin^2\bigg( \frac{nj\pi}{n+1} \bigg)}=\sqrt{\sum_{k=1}^n\sin^2\bigg( \frac{kj\pi}{n+1} \bigg)}$$

Iremos utilizar que $\sin^2(x)=\frac{1}{2}(1-\cos(2x))$ para escrever o somatório como: $$\sum_{k=1}^n\sin^2\bigg( \frac{kj\pi}{n+1} \bigg)=\sum_{k=1}^n\frac{1}{2}\bigg(1-\cos\bigg(2\frac{kj\pi}{n+1}\bigg)\bigg)$$

Podemos utilizar a linearidade para passar o fator multiplicativo para fora e separar em dois somatórios, sendo $\sum\limits_{k=1}^n1=n$:$$\sum_{k=1}^n\frac{1}{2}\bigg(1-\cos\bigg(2\frac{kj\pi}{n+1}\bigg)\bigg)=\frac{n}{2}-\frac{1}{2}\sum_{k=1}^n\cos\bigg(kj\frac{2\pi}{n+1}\bigg)$$

[^2]: Nos passos que se seguem, $i$ é a unidade imaginária, isto é, $i=\sqrt{-1}$.

Da fórmula de Euler, temos que $e^{i\theta}=\cos(\theta)+i\sin(\theta)$. Conforme [@Unity], este somatório de cossenos é a parte real das raízes $n+1$-ésimas da unidade, com exceção do próprio $1$. Podemos explorar melhor este somatório, fazendo: $$ \sum_{k=1}^n\cos\bigg(kj\frac{2\pi}{n+1}\bigg) = \sum_{k=0}^n\cos\bigg(kj\frac{2\pi}{n+1}\bigg) - 1$$

e, aplicando a fórmula de Euler: $$\sum_{k=0}^n\cos\bigg(kj\frac{2\pi}{n+1}\bigg) = \mathfrak{Re} \bigg\{ \sum_{k=0}^n\exp\bigg(kji\frac{2\pi}{n+1}\bigg) \bigg\} =\mathfrak{Re} \bigg\{ \sum_{k=0}^n\exp\bigg(ji\frac{2\pi}{n+1}\bigg)^k \bigg\} $$

Que é a soma de uma progressão geométrica, cujo resultado é conhecido: $$\mathfrak{Re} \bigg\{ \sum_{k=0}^n\exp\bigg(j\frac{2\pi}{n+1}\bigg)^k \bigg\} = \mathfrak{Re} \Bigg\{ \frac{1-\exp\big(\frac{j2\pi i}{n+1}\big)^{n+1}}{1-\exp\big(\frac{j2\pi i}{n+1}\big)} \Bigg\} = \mathfrak{Re} \Bigg\{ \frac{1-\exp\big(j2\pi i\big)}{1-\exp\big(\frac{j2\pi i}{n+1}\big)} \Bigg\}$$

Como $\exp(j2\pi i) = \cos(j2\pi) + i\sin(j2\pi) = 1$: $$\mathfrak{Re} \Bigg\{ \frac{1-\exp\big(j2\pi i\big)}{1-\exp\big(\frac{j2\pi i}{n+1}\big)} \Bigg\} = \mathfrak{Re}\{0\}=0$$

Portanto, teremos que o somatório original vale: $$ \sum_{k=1}^n\cos\bigg(kj\frac{2\pi}{n+1}\bigg) = \sum_{k=0}^n\cos\bigg(kj\frac{2\pi}{n+1}\bigg) - 1 = 0 -1 = -1$$

e, desta forma, como queríamos mostrar: $$\sum_{k=1}^n\frac{1}{2}\bigg(1-\cos\bigg(2\frac{kj\pi}{n+1}\bigg)\bigg) = \frac{n}{2}-\frac{1}{2}\sum_{k=1}^n\cos\bigg(kj\frac{2\pi}{n+1}\bigg) = \frac{n}{2}-\frac{1}{2}\times(-1)$$

Logo: $$\sqrt{\sum_{k=1}^n\sin^2\bigg( \frac{kj\pi}{n+1} \bigg)} = \sqrt{\frac{n+1}{2}}$$ Concluímos, então, que: $$\lVert \symbf{v}_j \rVert = \sqrt{\frac{n+1}{2}}$$
\begin{flushright}
$\Box$
\end{flushright}

Voltando ao questionamento inicial, encontramos a origem dos valores encontrados nas matrizes de proporcionalidade. Considerando os casos em que $n = 4$ e $n=8$, devemos ter:

* Para $P_4$, $p_{i,j} = \frac{v_{i,j}}{\tilde{v}_{i,j}} = \lVert \symbf{v}_j \rVert=\sqrt{\frac{5}{2}} \approx1,58113883$. 
* Para $P_8$, $p_{i,j} = \frac{v_{i,j}}{\tilde{v}_{i,j}} = \lVert \symbf{v}_j \rVert=\sqrt{\frac{9}{2}} \approx2,12132034$. 

Que são exatamente os valores encontrados. Não fizemos a mesma análise _no relatório_ para as demais matrizes, com $n=16$ e $n=32$ pois são muito grandes. Contudo, tais análises foram feitas em código, são obtíveis pela execução da CLI do programa e seguem exatamente a mesma tendência.

### Comparação do Número de Iterações

Executamos o Algoritmo QR, com e sem deslocamento espectral, para os casos em que $n=4,8,16$ e $32$, coletando o número de iterações necessárias para convergência. Sumarizamos os dados na Tabela \ref{table:3} abaixo.

\begin{table}[h!]
    \centering
    \begin{tabular}{|c|c|c|}
        \hline
        & \multicolumn{2}{|c|}{Número de Iterações} \\
        $n$ & \multicolumn{1}{c}{Com Deslocamento} & \multicolumn{1}{c|}{Sem Deslocamento}\\
        \hline
        4  & 9 & 45 \\
        8  & 19 & 143 \\
        16 & 35 & 473 \\
        32 & 66 & 1600 \\
        \hline
    \end{tabular}
    \caption{Número de iterações para diferentes tamanhos de matriz, com e sem deslocamento espectral.}
    \label{table:3}
\end{table}

Para os valores tabelados, é possível observar uma correlação numérica entre o valor de $n$ e o número de iterações necessárias para convergência. Para o uso de deslocamento espectral, todos os valores apresentados se aproximam de $2n$, com erro absoluto de 3, no máximo. Para a execução do algoritmo sem deslocamento espectral, todos os valores apresentados estão entre $1,5n^2$ e $3n^2$.

A fim de desenvolver uma análise mais profunda a respeito dessa correlação, executamos o Algoritmo QR para $n$ variando de $3$ a $128$, com e sem deslocamento espectral, e coletamos o número de iterações para cada teste. Com tais valores, foi possível produzir os gráficos da Figura \ref{fig:t1} abaixo.

\begin{figure}[h]
    \centering
    \includegraphics[width = \linewidth]{fig1}
    \caption{Gráficos do Número de Iterações necessárias para convergência com relação ao tamanho da matriz de entrada para o Algoritmo QR com e sem deslocamento espectral.}
    \label{fig:t1}
\end{figure}

Observamos, primeiramente, que, para todos os valores, o número de iterações quando utilizado deslocamento espectral é muito menor quando comparado ao caso sem utilização. A aceleração provida pelo deslocamento espectral, em particular quando feita pela heurística de Wilkinson, torna a convergência muito mais rápida. Se compararmos os valores encontrados, temos que, para $n = 128$, com deslocamento, o número de iterações é muito menor que 1000, ao passo que, sem deslocamento, são necessárias quase 20000 iterações para convergência. Nota-se o salto vertiginoso.

Constata-se que o número de iterações para o caso com deslocamento espectral cresce aproximadamente linear com o tamanho da matriz. De acordo com [@On], este comportamento é esperado pois o número de iterações é proporcional ao número de autovalores da matriz. Como o número de autovalores é igual ao tamanho da matriz, para o problema com que lidamos, observamos a tendência de $O(n)$ iterações com deslocamento espectral.

A convergência para o caso sem deslocamento espectral é muito mais lenta, sendo comparável a $O(n^2)$ ou $O(n^3)$ iterações. Sobre os aspectos gerais do Algoritmo QR, [@Tiozao] menciona que o algoritmo é estável, converge cubicamente para a maioria dos casos (quadraticamente, caso contrário), e tem um custo computacional de $O(n^3)$.

### Evolução do Erro por Iteração

Utilizando a definição de erro médio por iteração $E_{avg}^{(k)}$, definida anteriormente como: $$E_{avg}^{(k)}=\frac{1}{n}\sum\limits_{j=1}^n|\alpha^{(k)}_j-\lambda_j|$$

Executamos o Algoritmo QR com deslocamento espectral para $n = 512$, calculando, a cada iteração, o erro por autovalor, dado que temos o valor correto dos autovalores. A partir dos pontos obtidos, foi criada a Figura \ref{fig:t2}.

\begin{figure}[h]
    \centering
    \includegraphics[width = \linewidth]{fig2}
    \caption{Evolução do Erro Médio dos Autovalores por Iteração.}
    \label{fig:t2}
\end{figure}

Na Figura, podemos observar duas importantes características do algoritmo. Primeiramente, devemos dar especial atenção às últimas iterações, pois delas podemos melhor analisar a convergência de todos os autovalores, uma vez que converge-se primeiro $\alpha_j$ para, após disso, convergir $\alpha_{j-1}$. Assim, na maior parte das iterações o erro médio está na ordem de $10^{-4}$, indo rapidamente para $10^{-12}$ conforme as submatrizes ficam menores. Além disso, observamos que a convergência para cada autovalor é muito ágil, bastando poucas iterações por autovalor para convergir, o que sustenta o que se desenvolveu anteriormente sobre a complexidade do número de iterações.

## Resultados para os Testes 2 e 3

Desejamos solucionar o sistema mecânico composto por $m$ massas e $m+1$ molas. Todas as massas têm $2$ kg e as constantes elásticas das molas são dadas por $k_i$ N/m, $i=1,\dots,n+1$\. O deslocamento de cada massa em relação ao equilíbrio é dado por $x_i(t)$ para um determinado instante $t$. As massas se movimentam sem quaisquer perdas sobre uma superfície plana sem amortecimento, estando o sistema anexo a dois anteparos. Desta forma, considerando a dinâmica do sistema, podemos descrevê-lo pelo sistema de equações diferenciais $$X''(t)+AX(t)=0$$ Sendo $X(t)=(x_1(t),\dots,x_n(t))^T$ e $A$, que chamaremos de _matriz dos coeficientes_, a matriz: $$A = \frac{1}{m}
    \begin{bmatrix}
        k_1+k_2 & -k_2 & & \\
        -k_2 & k_2+k_3 & \ddots & \\
        & \ddots & \ddots & -k_n \\
        & & -k_n & k_n+k_{n+1}
    \end{bmatrix}
$$

que é real e tridiagonal simétrica. Vamos considerar a decomposição de $A$ em sua matriz de autovetores e matriz diagonal de autovalores, isto é, iremos fatorar $A=Q\Lambda Q^T$, sendo $Q$ a matriz ortogonal cujas colunas são autovetores de $A$ e $\Lambda$ a matriz diagonal cujas entradas são autovalores de $A$. Podemos, portanto, reescrever o sistema de equações diferenciais como: $$X''(t)+Q\Lambda Q^T X(t)=0$$ Multiplicando o sistema por $Q^T$ e fazendo a substituição $Y(t)=Q^T X(t)$, temos: $$Y''(t)+\Lambda Y(t)=0$$ O sistema acima é mais interessante que o primeiro em relação a sua solução, pois é composto por `n` equações diferenciais desacopladas cuja solução é conhecida e intrinsicamente relacionada aos autovalores de $A$. O sistema é do tipo: $$
\begin{bmatrix}
    y''_1(t) \\ y''_2(t) \\ \vdots \\ y''_n(t)
\end{bmatrix} + 
\begin{bmatrix}
    \lambda_1 & & &\\
    & \lambda_2 & &\\
    & & \ddots &\\
    & & & \lambda_n\\
\end{bmatrix}\begin{bmatrix}
    y_1(t) \\ y_2(t) \\ \vdots \\ y_n(t)
\end{bmatrix} =0
$$

Portanto, devemos resolver EDOs da forma $y''_j(t)=-\lambda_jy_j(t)$. Sua solução é conhecida e dada por:$$
\begin{matrix}
    y_j(t)=a_j\cos(\omega_j t)+b_j\sin(\omega_j t) & \omega_j=\sqrt{\lambda_j}
\end{matrix}
$$

Vale ressaltar que, ao associar à solução-geral do sistema um conjunto de condições inicias (c.i.s), a solução fica unicamente definida para o problema em questão. Por exemplo, pode-se tomar $X(0)$ e $X'(0)$. Nota-se que condições iniciais em $X$ não são diretamente operáveis em $Y$ e devem passar pela transformação original.

Em ambos os sistemas massa-mola a se resolver, temos que $X'(0)=0$.  $Q$ é inversível, logo $Q$ é bijetora. Portanto, $0_X\mapsto 0_Y$ e, consequentemente, $Y'(0)=0$. $Y'(t)$ é dada por suas componentes $y'_j(t)=-a_j\omega_j\sin(\omega_j t)+b_j\omega_j\cos(\omega_j t)$. Como $Y'(0)=0$, $y'_j(0)=-a_j\omega_j\sin(0)+b_j\omega_j\cos(0)=b_j\omega_j=0\iff b_j=0$. Com isso, resta apenas o termo cossenoidal nas funções de $Y$, isto é, $y_j(t)=a_j\cos(\omega_j t)$. 

Para encontrar $Y(0)$, devemos fazer $Y(0)=Q^TX(0)$, sendo $Y(0)=(a_1, a_2, \dots, a_n)^T$ e, portanto: $(a_1, a_2, \dots, a_n)^T = Q^TX(0)$. Definido $Y(t)$ pelo problema de c.i.s, podemos reconstruir $X(t)$ fazendo $X(t)=QY(t)$.

### Solução para o Teste 2

Para o teste 2, desejamos solucionar o sistema mecânico composto por $5$ massas e $6$ molas. Todas as massas têm $2$ kg e as constantes elásticas das molas são dadas por $k_i=(40+2i)$ N/m, $i=1,\dots,6$\. Sua matriz $A$ é, portanto: $$A =
    \begin{bmatrix}
        43 & -22 &   0 &   0 &   0 \\
        -22 &  45 & -23 &   0 &   0 \\
        0 & -23 &  47 & -24 &   0 \\
        0 &   0 & -24 &  49 & -25 \\
        0 &   0 &   0 & -25 &  51 \\
    \end{bmatrix}
$$

Executando o Algoritmo QR com deslocamento espectral sobre $A$, encontramos seus autovalores, que podem ser transformados nas frequências do sistema. Os autovalores de $A$, bem como sua conversão em frequência, estão na Tabela \ref{table:4}.

\begin{table}[h!]
    \centering
    \begin{tabular}{|c|c|c|}
        \hline
        $j$ & $\lambda_j$ & $\omega_j$ \\
        \hline
            $1$ & $88,445006$ & $9,404520$ \\
            $2$ & $70,113831$ & $8,373400$ \\
            $3$ & $46,773186$ & $6,839093$ \\
            $4$ & $23,398633$ & $4,837213$ \\
            $5$ & $6,269344$  & $2,503866$ \\
        \hline
    \end{tabular}
    \caption{Autovalores para a matriz do teste 2.}
    \label{table:4}
\end{table}

Os autovetores de $A$ compõem os _modos naturais de vibração_ do sistema. Isto é, quando o deslocamento inicial é paralelo a um dos autovetores de $A$, todas as massas entram em ressonância, vibrando com a **mesma frequência**! A matriz abaixo exibe os autovetores para este teste.

$$
Q =
\begin{bmatrix}
     0,189335 &  0,474976 &  0,598822 &  0,532504 & 0,310486 \\
     0,391105 & -0,585383 & -0,102703 &  0,474446 & 0,518379 \\
     0,557661 &  0,184858 & -0,564869 & -0,063757 & 0,575934 \\
     0,588202 &  0,38296  &  0,093085 & -0,517375 & 0,480644 \\
     0,392711 & -0,500893 &  0,550565 & -0,468614 & 0,268632 \\
\end{bmatrix}
$$ \begin{center} Matriz dos Autovetores para o problema do Teste 2 \end{center}

**Primeiro Conjunto de Condições Iniciais**  Ao utilizar as condições iniciais $X(0)=(−2,−3,−1,−3,−1)^T$, teremos $Y(0)=Q^TX(0)=(1.608881, -0.026647, -1.154489, -0.403849, -4.462606)^T$. Como $Y(0)=(a_1, a_2, a_3, a_4, a_5)^T$, as funções ficam definidas como: $$Y(t)=\begin{bmatrix}
     1.608881  \cos (9,404520 t) \\
     -0.026647 \cos(8,373400 t) \\
     -1.154489 \cos(6,839093 t) \\
     -0.403849 \cos(4,837213 t) \\
     -4.462606 \cos(2,503866 t) \\
\end{bmatrix}$$

Desta forma, recuperamos $X(t)$ por $X(t)=QY(t)$, do qual obtemos: $$\begin{bmatrix}
     x_1(t) \\ x_2(t) \\ x_3(t) \\ x_4(t) \\ x_5(t)
\end{bmatrix}=\begin{bmatrix}
 0.304617 & - 0.012657 & - 0.691334 & - 0.215051 & - 1.385575 \\
-0.629242 & 0.015599 & 0.118570 & - 0.191604 & - 2.313322 \\
 0.897210 & - 0.004926 & 0.652134 & 0.025748 & - 2.570167 \\
-0.946347 & - 0.010205 & - 0.107466 & 0.208941 & - 2.144924 \\
 0.631824 &  0.013347 & - 0.635621 & 0.189249 & - 1.198800 \\
\end{bmatrix}\begin{bmatrix}
     \cos(9.404520 t) \\ \cos(8.373400 t) \\ \cos(6.839093 t) \\ \cos(4.837213 t) \\ \cos(2.503866 t) 
\end{bmatrix}
$$

A Figura \ref{fig:t3} exibe a evolução da solução, com o deslocamento de cada massa por 10 segundos.

\begin{figure}[h]
    \centering
    \includegraphics[width = \linewidth]{fig3}
    \caption{Evolução da Solução para o primeiro conjunto de $X(0)$.}
    \label{fig:t3}
\end{figure}

O deslocamento resultante de cada mola é a superposição de componentes cossenodais de diferentes frequências. Observamos um evidente comportamento vinculado entre as oscilações de cada massa, além de percebermos nas oscilações que o harmônico de menor frequência é o mais relevante, de maior amplitude, sendo notável sua contribuição para o formato do sinal. Para as massas conectadas a molas de maior coeficiente elástico, percebe-se maior tendência na mudança de sentido do movimento. A massa do meio é aquela com maior amplitude de deslocamento.

\pagebreak

**Segundo Conjunto de Condições Iniciais**  Neste caso, as condições iniciais das EDOs são $X(0)=(1, 10, -4, 3, -2)^T$, assim $Y(0)=Q^TX(0)=(-8.502389, -3.967618,  1.009392,  4.917091,  4.095209)^T$. Teremos que o vetor de funções fica $$Y(t)=\begin{bmatrix}
     -8.50238\cos (9,404520 t) \\
     -3.967618\cos(8,373400 t) \\
     1.009392\cos(6,839093 t) \\
     4.917091\cos(4,837213 t) \\
     4.095209\cos(2,503866 t) \\
\end{bmatrix}$$

Assim, encontramos a solução $X(t)$ por $X(t)=QY(t)$, resultando em: $$\begin{bmatrix}
     x_1(t) \\ x_2(t) \\ x_3(t) \\ x_4(t) \\ x_5(t)
\end{bmatrix}=\begin{bmatrix}
-1.609797 & - 1.884524 &  0.604447 &  2.618370 &  1.271504 \\
3.325329 &  2.322575 & - 0.103668 &  2.332893 &  2.122871 \\
-4.741452 & - 0.733445 & - 0.570174 & - 0.313499 &  2.358570 \\
5.001122 & - 1.519437 &  0.093960 & - 2.543981 &  1.968336 \\
-3.338978 &  1.987353 &  0.555736 & - 2.304216 &  1.100105 \\
\end{bmatrix}\begin{bmatrix}
     \cos(9.404520 t) \\ \cos(8.373400 t) \\ \cos(6.839093 t) \\ \cos(4.837213 t) \\ \cos(2.503866 t) 
\end{bmatrix}
$$

Apresentamos a evolução do sistema abaixo na Figura \ref{fig:t4}

\begin{figure}[h]
    \centering
    \includegraphics[width = \linewidth]{fig4}
    \caption{Evolução da Solução para o segundo conjunto de $X(0)$.}
    \label{fig:t4}
\end{figure}

Podemos comparar estas condições iniciais com as do primeiro conjunto. Notamos que diferentes c.i.s. implicam uma evolução do sistema consideravelmente distinta. A alta amplitude do deslocamento inicial produz oscilações proporcionalmente vigorosas e relativamente rápidas. Essa consideração de rapidez é refletida pelo fato de a componente espectral de maior amplitude ser aquela de maior frequência. A massa central tem o padrão mais comportado de oscilação. Por fim, notamos que as massas conectadas a molas de maiores constantes elásticas tem comportamento oscilatório mais acentuado, bem como maior aceleração.

**Terceiro Conjunto de Condições Iniciais**  Deve-se utilizar como c.i. o autovetor associado ao autovalor de maior frequência. Portanto, $X(0)=( 0.189335, -0.391105,  0.557661, -0.588202, 0.392711)^T$. Temos portanto $Y(0)=Q^TX(0)=(1, 0, 0, 0, 0)^T$. Logo, $X(t)$ é o vetor $$X(t)=\begin{bmatrix}
 0.189335 \cos(9.404520 t) \\
-0.391105 \cos(9.404520 t) \\
0.557661 \cos(9.404520 t) \\
-0.588202 \cos(9.404520 t) \\
0.392711 \cos(9.404520 t) \\
\end{bmatrix}$$

O gráfico do deslocamento das massas no tempo está na Figura \ref{fig:t5} abaixo.

\begin{figure}[h]
    \centering
    \includegraphics[width = \linewidth]{fig5}
    \caption{Evolução da Solução para o terceiro conjunto de $X(0)$.}
    \label{fig:t5}
\end{figure}

Como o vetor de deslocamento inicial é um autovetor da matriz de coeficientes do sistema, temos que as massas entram em ressonância e vibram à mesma frequência. Isso ocorre pois, nesta situação, $X(0)$ é um modo natural de vibração do sistema e a amplitude das demais componentes cossenodais é nula. O valor encontrado para $Y(0)$ é o vetor $X(0)$ na base de autovetores de $A$, sendo, nesta situação, paralelo a um dos vetores que formam a base.

\pagebreak

### Solução para o Teste 3

Neste teste, soluciona-se um sistema de $10$ massas e $11$ molas. Igualmente, as massas têm $2$ kg, mas as constantes são dadas por $k_i=(40+2(-1)^i)$ N/m, $i=1,\dots,11$\. A matriz de coeficientes é: $$A =
    \begin{bmatrix}
  40& -21&   0&   0&   0&   0&   0&   0&   0&   0\\
 -21&  40& -19&   0&   0&   0&   0&   0&   0&   0\\
   0& -19&  40& -21&   0&   0&   0&   0&   0&   0\\
   0&   0& -21&  40& -19&   0&   0&   0&   0&   0\\
   0&   0&   0& -19&  40& -21&   0&   0&   0&   0\\
   0&   0&   0&   0& -21&  40& -19&   0&   0&   0\\
   0&   0&   0&   0&   0& -19&  40& -21&   0&   0\\
   0&   0&   0&   0&   0&   0& -21&  40& -19&   0\\
   0&   0&   0&   0&   0&   0&   0& -19&  40& -21\\
   0&   0&   0&   0&   0&   0&   0&   0& -21&  40\\
    \end{bmatrix}
$$

A aplicação dos Algoritmo QR em $A$ retorna seus autovalores, cujas raízes são as frequências de vibração do sistema massa-mola. Os autovalores e as respectivas frequências são exibidos na Tabela \ref{table:5} abaixo.

\begin{table}[h!]
    \centering
    \begin{tabular}{|c|c|c|}
        \hline
        $j$ & $\lambda_j$ & $\omega_j$ \\
        \hline
            $1$ & $78.398748$ &  $8.854307$ \\
            $2$ & $73.729928$ &  $8.586613$ \\
            $3$ & $66.393758$ &  $8.148237$ \\
            $4$ & $57.063951$ &  $7.554068$ \\
            $5$ & $47.001372$  & $6.855755$ \\
            $6$ & $32.998628$  & $5.744443$ \\
            $7$ & $22.936049$  & $4.789160$ \\
            $8$ & $13.606242$  & $3.688664$ \\
            $9$ & $ 6.270072$  & $2.504011$ \\
            $10$& $ 1.601252$  & $1.265406$ \\
        \hline
    \end{tabular}
    \caption{Autovalores para a matriz do teste 2.}
    \label{table:5}
\end{table}

\pagebreak
Os modos naturais de vibração do sistema são constituídos pelos autovetores de $A$. A matriz de autovetores para o teste atual está abaixo, estando armazenados em suas colunas.

$$
Q =
\begin{bmatrix}
 0.125245 &  0.240074 &  0.334551 & \cdots &  0.398212 & 0.334551  & 0.240074  & 0.125245 \\
-0.229012 & -0.385603 & -0.420479 & \cdots &  0.323575 & 0.420479  & 0.385603  & 0.229012 \\
 0.324402 &  0.419201 &  0.214339 & \cdots & -0.149526 & 0.214339  & 0.419201  & 0.324402 \\
-0.385971 & -0.324437 &  0.111042 & \cdots & -0.414258 & -0.111042 &  0.324436 &  0.385971 \\
 0.421493 &  0.112631 & -0.391154 & \cdots & -0.206781 & -0.391154 &  0.112631 &  0.421493 \\
-0.421493 &  0.112631 &  0.391154 & \cdots &  0.206781 & -0.391154 & -0.112631 &  0.421493 \\
 0.385971 & -0.324437 & -0.111042 & \cdots &  0.414258 & -0.111042 & -0.324437 &  0.385971 \\
-0.324402 &  0.419201 & -0.214339 & \cdots &  0.149526 & 0.214339  &-0.419201  & 0.324402 \\
 0.229012 & -0.385603 &  0.420479 & \cdots & -0.323575 & 0.420479  &-0.385603  & 0.229012 \\
-0.125245 &  0.240074 & -0.334551 & \cdots & -0.398212 & 0.334551  &-0.240074  & 0.125245 \\
\end{bmatrix}
$$ \begin{center} Matriz dos Autovetores para o problema do Teste 3 \end{center}

**Primeiro Conjunto de Condições Iniciais**  Ao utilizar as condições iniciais $X(0)=(−2,−3,−1,−3,−1,$ $−2,−3,−1,−3,−1)^T$, teremos $Y(0)=Q^TX(0)=(0.296248  ,2.363722 ,-0.725705 ,-0.819342 ,-0.05771,$ $-0.720211,-0.604993, -2.11549  , -0.127443 , -5.978921)^T$. Utilizando as condições iniciais dadas, teremos: $$Y(t)=\begin{bmatrix}
0.2962480\cos(8.854307 t) \\  
2.3637220\cos(8.586613 t) \\ 
-0.725705\cos(8.148237 t) \\ 
-0.819342\cos(7.554068 t) \\ 
-0.057710\cos(6.855755 t) \\  
-0.720211\cos(5.744443 t) \\ 
-0.604993\cos(4.789160 t) \\
-2.115490\cos(3.688664 t) \\  
-0.127443\cos(2.504011 t) \\ 
-5.978921\cos(1.265406 t) \\
\end{bmatrix}$$

Pode-se reconstruir $X(t)$ utilizando $X(t)=QY(t)$. Por brevidade, omitiremos a matriz, que pode ser vista na execução da CLI. A Figura \ref{fig:t6} exibe a evolução da solução para 10 massas, com o deslocamento de cada uma, por 10 segundos.

\begin{figure}[H]
    \centering
    \includegraphics[width = \linewidth]{fig6}
    \includegraphics[width = \linewidth]{fig9}
    \caption{Evolução da Solução para o primeiro conjunto de $X(0)$.}
    \label{fig:t6}
\end{figure}

Como anteriormente, as soluções do sistema são superposições de cossenos. Observamos que existem pares de massas com oscilações similares. Por exemplo, os deslocamentos das massas 5 e 6 têm formatos muito parecidos, bem como os de 3 e 8. Além disso, massas em posições centrais tendem a oscilar mais regularmente que aquelas nas pontas, apesar de possuírem amplitude maior. A componente espectral de menor frequência é a mais relevante para a oscilação, pois tem a maior das amplitudes, o que é muito perceptível. Dada a pequena diferença entre os coeficientes das molas utilizadas, é possível notar uma diminuição na variação das amplitudes dos deslocamentos das massas, se comparadas às do teste 2. Todas as massas possuem deslocamento absoluto em torno de 2 a 4 centímetros.

**Segundo Conjunto de Condições Iniciais**  Para este conjunto de condições iniciais, temos $X(0)=(1, 10, -4, 3, -2, 1, 10, -4, 3, -2)^T$, gerando $Y(0)=Q^TX(0)=(  0.209968,-12.936835, -1.543527,  2.183659,$ $-2.489729,  0.853680,  6.979812,  2.364570,  0.810496,  4.852832)^T$. A função em $Y$ se torna: $$Y(t)=\begin{bmatrix}
  0.209968\cos(8.854307 t) \\
-12.936835\cos(8.586613 t) \\
 -1.543527\cos(8.148237 t) \\
  2.183659\cos(7.554068 t) \\
 -2.489729\cos(6.855755 t) \\
  0.853680\cos(5.744443 t) \\
  6.979812\cos(4.789160 t) \\
  2.364570\cos(3.688664 t) \\
  0.810496\cos(2.504011 t) \\
  4.852832\cos(1.265406 t) \\
\end{bmatrix}$$

A partir das c.i.s e dos autovalores do sistema, encontramos a solução $X(t)$ por $X(t)=QY(t)$. Dela, pode-se analisar a oscilação das massas e a evolução do sistema. Traçaram-se gráficos para essa solução, que se encontram na Figura \ref{fig:t7} abaixo.

\begin{figure}[H]
    \centering
    \includegraphics[width = \linewidth]{fig7}
    \includegraphics[width = \linewidth]{fig10}
    \caption{Evolução da Solução para o segundo conjunto de $X(0)$.}
    \label{fig:t7}
\end{figure}

Comparando com a solução anterior, observamos que a amplitude das oscilações é muito maior, o que reflete as condições iniciais. Também observamos que existem pares de massas com funções-deslocamento similares. As massas centrais e da bordam oscilam muito, ao passo que naquelas de posições intermediárias predominam componentes de baixa frequência, embora a amplitude de oscilação seja maior. Os formatos dos deslocamentos das massas são muito mais similares a cossenos puros, ao passo que, na solução anterior, há um perfil similar ao de um batimento.

\pagebreak

**Terceiro Conjunto de Condições Iniciais** Utilizaremos, para o tempo $t=0$, o autovetor associado à oscilação de maior frequência, i.e., cujo autovalor tem maior módulo.

Encontramos $X(0)=(0.125245, -0.229012,  0.324402, -0.385971,  0.421493, -0.421493,  0.385971, -0.324402,$ $0.229012, -0.125245)^T$, a partir do qual, por ser autovetor de $A$, na base ortonormal dos autovetores de $A$, obtemos $Y(0)=Q^TX(0)=(1, 0, \cdots, 0, 0)^T$. Cada $y_j(t)$ é composto por um único cosseno, cuja frequência é o autovalor de maior módulo: $$X(t)=\begin{bmatrix}
         0.125245 \cos(8.854307 t) \\
        -0.229012 \cos(8.854307 t) \\
         0.324402 \cos(8.854307 t) \\
        -0.385971 \cos(8.854307 t) \\
         0.421493 \cos(8.854307 t) \\
        -0.421493 \cos(8.854307 t) \\
         0.385971 \cos(8.854307 t) \\
        -0.324402 \cos(8.854307 t) \\
         0.229012 \cos(8.854307 t) \\
        -0.125245 \cos(8.854307 t)
\end{bmatrix}$$

Como $X(0)$ é autovetor, apenas um termo de $Y(t)$ não é nulo. Foram construídos os gráficos do deslocamento para cada massa por 10 segundos, exibidos na Figura \ref{fig:t8} abaixo.

Tal qual esperado, as massas oscilam todas sob mesma frequência e com amplitude igual ao deslocamento inicial. Pode-se notar um padrão alternante entre os valores iniciais dos cossenos. Para cada massa de índice ímpar, seu deslocamento inicial é positivo e, para massas com índice par, o deslocamento inicial é negativo. Além disso, há uma tendência de aumento da amplitude conforme a massa se afasta dos anteparos.

Na implementação destes problemas, como não é possível trabalhar simbolicamente com os cossenos, armazenamos as frequências de oscilação na forma de autovalores. Quando necessária o cálculo de $X(t)$ para um determinado $t$, tomamos as frequências, calculando a raiz, termo a termo, do vetor de autovalores, compondo seu produto com o tempo aos cossenos. Assim, construímos um vetor $C(t)$ dado por $c_j=cos(\sqrt{\lambda_j} t)$, que é aplicado na forma $X(t)=QY(0)\circ C(t)$, em que $\circ:\mathbb{R}^{n\times n}\times\mathbb{R}^{n\times n}\rightarrow\mathbb{R}^{n\times n}$ é o produto de Hadamard.

\begin{figure}[H]
    \centering
    \includegraphics[width = \linewidth]{fig8}
    \includegraphics[width = \linewidth]{fig11}
    \caption{Evolução da Solução para o terceiro conjunto de $X(0)$.}
    \label{fig:t8}
\end{figure}

Todos os testes envolvendo molas, para diferentes condições iniciais, podem ser vistos **com animação** pela execução da rotina 5 na CLI. A animação exibe os gráficos, atualizados tempo a tempo, e as massas oscilando com o tempo. A posição inicial de cada massa (à qual se aplica o deslocamento) foi escolhida para corresponder à melhor visualização, pois não interfere na solução.  -->

\pagebreak

# Referências {-}

\setstretch{1}
