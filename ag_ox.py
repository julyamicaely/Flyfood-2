import random as rd

def ler_arquivo_matriz(nome_arquivo):
    with open(nome_arquivo, "r") as arquivo:
        return arquivo.read()

def formular_matriz(conteudo_matriz):
    linhas = conteudo_matriz.strip().split("\n")
    n_linhas, n_colunas = map(int, linhas[0].split())
    matriz = []
    for i in range(1, n_linhas + 1):
        valores = linhas[i].split()
        matriz.append(valores)
    return matriz

def encontrar_pontos(matriz):  # encontrar pontos de interesse
    loc_pontos = {}
    for i in range(len(matriz)):
        for j in range(len(matriz[0])):
            valor = matriz[i][j]
            if valor != '0':
                loc_pontos[valor] = (i, j)
    return loc_pontos

def criar_populacao(tam_populacao, ind_visitaveis):
    populacao = []
    pontos_meio= [p for p in ind_visitaveis if p != 'R'] # Remove 'R' para embaralhar apenas os intermediários
    for _ in range(tam_populacao):
        individuo = ['R'] + rd.sample(pontos_meio, len(pontos_meio)) + ['R']
        populacao.append(individuo)
    return populacao

def dist_manhattan(tupleA, tupleB):
    dist = abs(tupleA[0] - tupleB[0]) + abs(tupleA[1] - tupleB[1])
    return dist

def calcular_aptidoes(loc_pontos, rotas): # calcular as distâncias do ponto de partida até os outros pontos
    distancias = []
    for rota in rotas:
        soma = 0
        for a, b in zip(rota, rota[1:]):
            soma += dist_manhattan(loc_pontos[a], loc_pontos[b])
        distancias.append((rota, soma))
    return distancias

def selecao_por_ranqueamento(populacao, aptidoes, k=1):
    # Ordena por aptidão (menor distância = melhor)
    ordenado = sorted(zip(populacao, aptidoes), key=lambda x: x[1])
    n = len(ordenado)

    # Atribui probabilidade de forma linear crescente (melhor tem mais chance)
    # Ex: se n = 4 → pesos: [1, 2, 3, 4]
    pesos = list(range(1, n + 1))
    total = sum(pesos)
    probs = [p / total for p in pesos]

    # Acumulada
    acumulada = []
    atual = 0
    for p in probs:
        atual += p
        acumulada.append(atual)

    # Seleciona k indivíduos com base em roleta por rank
    selecionados = []
    for _ in range(k):
        r = rd.random()
        for i, limite in enumerate(acumulada):
            if r <= limite:
                selecionados.append(ordenado[i][0])
                break
    return selecionados[0] if k == 1 else selecionados
        
def selecionar_nova_geracao(populacao, aptidoes, tamanho_nova_geracao):
    nova_geracao = []
    for _ in range(tamanho_nova_geracao):
        individuo = selecao_por_ranqueamento(populacao, aptidoes)
        nova_geracao.append(individuo)
    return nova_geracao

def ox_crossover(pai1, pai2):
    tam = len(pai1)
    filho1 = ['R'] + [None] * (tam - 2) + ['R']
    filho2 = ['R'] + [None] * (tam - 2) + ['R']

    cx1, cx2 = sorted(rd.sample(range(1, tam - 1), 2))

    filho1[cx1:cx2+1] = pai1[cx1:cx2+1]
    filho2[cx1:cx2+1] = pai2[cx1:cx2+1]

    def preencher(filho, outro_pai):
        genes_faltantes = [gene for gene in outro_pai[1:-1] if gene not in filho]
        pos = (cx2 + 1) if (cx2 + 1) < (tam - 1) else 1 

        for gene in genes_faltantes:
            while filho[pos] is not None:
                pos += 1
                if pos == tam - 1:
                    pos = 1
            filho[pos] = gene
        return filho

    filho1 = preencher(filho1, pai2)
    filho2 = preencher(filho2, pai1)

    return filho1, filho2

def mutacao_swap(rota, taxa=0.1):
    nova = rota[:]
    if rd.random() < taxa:
        i, j = rd.sample(range(1, len(nova)-1), 2)
        nova[i], nova[j] = nova[j], nova[i]
    return nova

if __name__ == "__main__":
    conteudo = ler_arquivo_matriz("arquivo_4.txt")
    matriz = formular_matriz(conteudo)
    pontos = encontrar_pontos(matriz)

    tam_popula = 40
    nome_pontos = list(pontos.keys())
    populacao = criar_populacao(tam_popula, nome_pontos)
    aptidoes = calcular_aptidoes(pontos, populacao)

    num_geracoes = 70
    for geracao in range(num_geracoes):
        apt_geracao = calcular_aptidoes(pontos, populacao)
        valores_aptidoes = [aptidao for _, aptidao in apt_geracao]
        aptidoes_invertidas = [1/(a+1e-6) for a in valores_aptidoes]

        melhor_rota, melhor_aptidao = min(apt_geracao, key=lambda x: x[1])  
        print(f"Geração {geracao+1} - Melhor rota: {melhor_rota} | Distância: {melhor_aptidao}")

        nova_popula = []
        nova_popula.append(melhor_rota[:])
        while len(nova_popula) < tam_popula:
            pai1 = selecao_por_ranqueamento(populacao, aptidoes_invertidas)
            pai2 = selecao_por_ranqueamento(populacao, aptidoes_invertidas)
            while pai2 == pai1:
                pai2 = selecao_por_ranqueamento(populacao, aptidoes_invertidas)
            filho1, filho2 = ox_crossover(pai1, pai2)
            filho1 = mutacao_swap(filho1, 0.01)
            filho2 = mutacao_swap(filho2, 0.02)
            nova_popula.append(filho1)
            if len(nova_popula) < tam_popula:
                nova_popula.append(filho2)
        
        populacao = nova_popula
