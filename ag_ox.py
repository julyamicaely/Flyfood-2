import random as rd
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def ler_arquivo_tsp(nome_arquivo):
    objArq = open("edgesbrasil58.tsp")
    distancias = {}

    for i in range(1, 58): #linhas 1 a 57 pois a 58 nao terá aresta
        linha = objArq.readline() #le só uma linha do arquivo
        #transformando a linha em lista de strings:
        lista = linha.split() #obs: lista de strings (não int)

        for j in range(i+1, 59): #colunas i+1 a 58
            if len(lista) > 0:
                peso = int(lista.pop(0)) #obs: peso int, poderia ser float em outro problema
            else:
                print(f"Erro! linha {i} do arquivo não possui elementos suficientes")
                exit()
            #gravando a aresta em (i, j) e (j, i):
            distancias[(i,j)] = peso
            distancias[(j,i)] = peso
    objArq.close()
    return distancias

def calcular_aptidoes_tsp(distancias, rotas):
    distancias_totais = []
    for rota in rotas:
        soma = 0
        for a, b in zip(rota, rota[1:]):
            soma += distancias[(a, b)]
        distancias_totais.append((rota, soma))
    return distancias_totais

def custo_caminho_tsp(permutacao, dicDistancias):
	#ex: permutacao = [5, 14, 2, 3, 7, ...]
	soma = 0
	for i in range(len(permutacao)-1):
		a = permutacao[i]
		b = permutacao[i+1]
		if (a,b) in dicDistancias:
			soma += dicDistancias[(a,b)]
		else:
			print(f"Erro! ({a},{b}) não existe no dicionario!")
			exit()
	soma += dicDistancias[(permutacao[-1],permutacao[0])]
	return soma

def inicializa_populacao_tsp(tamanho, qtdeCidades):
	#criando uma lista com "tamanho" permutacoes aleatorias de cidades:
	lista = []
	for i in range(tamanho):
		individuo = list(range(1, qtdeCidades+1))
		rd.shuffle(individuo)
		lista.append(individuo)
	return lista

def calcula_aptidao_tsp(populacao):
	listaAptidao = []
	for elem in populacao:
		listaAptidao.append(custoCaminho(elem, distancias))
	return listaAptidao


def ox_crossover_tsp(pai1, pai2):
    tam = len(pai1)
    filho1 = [None] * tam
    filho2 = [None] * tam

    cx1, cx2 = sorted(rd.sample(range(tam), 2))

    filho1[cx1:cx2+1] = pai1[cx1:cx2+1]
    filho2[cx1:cx2+1] = pai2[cx1:cx2+1]

    def preencher(filho, outro_pai):
        pos = (cx2 + 1) % tam
        for gene in outro_pai:
            if gene not in filho:
                while filho[pos] is not None:
                    pos = (pos + 1) % tam
                filho[pos] = gene
        return filho

    filho1 = preencher(filho1, pai2)
    filho2 = preencher(filho2, pai1)

    return filho1, filho2

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

def encontrar_pontos(matriz):
    loc_pontos = {}
    for i in range(len(matriz)):
        for j in range(len(matriz[0])):
            valor = matriz[i][j]
            if valor != '0':
                loc_pontos[valor] = (i, j)
    return loc_pontos

def criar_populacao(tam_populacao, ind_visitaveis):
    populacao = []
    pontos_meio= [p for p in ind_visitaveis if p != 'R'] 
    for _ in range(tam_populacao):
        individuo = ['R'] + rd.sample(pontos_meio, len(pontos_meio)) + ['R']
        populacao.append(individuo)
    return populacao

def dist_manhattan(tupleA, tupleB):
    dist = abs(tupleA[0] - tupleB[0]) + abs(tupleA[1] - tupleB[1])
    return dist

def calcular_aptidoes(loc_pontos, rotas):
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

    pesos = list(range(1, n + 1))
    total = sum(pesos)
    probs = [p / total for p in pesos]

    acumulada = []
    atual = 0
    for p in probs:
        atual += p
        acumulada.append(atual)

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

def main_matriz(nome_arquivo):   
    conteudo = ler_arquivo_matriz(nome_arquivo)
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
    else:
        return ' '.join(melhor_rota), melhor_aptidao

def main_tsp(nome_arquivo):
    distancias = ler_arquivo_tsp(nome_arquivo)
    qtde_cidades = 58
    tam_popula = 40

    populacao = inicializa_populacao_tsp(tam_popula, qtde_cidades)

    num_geracoes = 70
    for geracao in range(num_geracoes):
        aptidoes = [custo_caminho_tsp(ind, distancias) for ind in populacao]
        aptidoes_invertidas = [1/(a+1e-6) for a in aptidoes]

        melhor_idx = aptidoes.index(min(aptidoes))
        melhor_rota = populacao[melhor_idx]
        melhor_aptidao = aptidoes[melhor_idx]
        print(f"Geração {geracao+1} | {melhor_rota}: {melhor_aptidao}")

        nova_popula = [melhor_rota[:]]
        while len(nova_popula) < tam_popula:
            pai1 = selecao_por_ranqueamento(populacao, aptidoes_invertidas)
            pai2 = selecao_por_ranqueamento(populacao, aptidoes_invertidas)
            while pai2 == pai1:
                pai2 = selecao_por_ranqueamento(populacao, aptidoes_invertidas)
            filho1, filho2 = ox_crossover_tsp(pai1, pai2)
            filho1 = mutacao_swap(filho1, 0.3)
            filho2 = mutacao_swap(filho2, 0.4)
            nova_popula.append(filho1)
            if len(nova_popula) < tam_popula:
                nova_popula.append(filho2)
        populacao = nova_popula

def algoritmo_genetico(distancias, qtde_cidades=58, tam_popula=30, num_geracoes=70):
    populacao = inicializa_populacao_tsp(tam_popula, qtde_cidades)
    melhores_p_geracao = []

    for geracao in range(num_geracoes):
        aptidoes = [custo_caminho_tsp(ind, distancias) for ind in populacao]
        melhor_aptidao = min(aptidoes)
        melhores_p_geracao.append(melhor_aptidao)

        aptidoes_invertidas = [1/(a+1e-6) for a in aptidoes]
        melhor_idx = aptidoes.index(melhor_aptidao)
        melhor_rota = populacao[melhor_idx]

        nova_popula = [melhor_rota[:]]
        while len(nova_popula) < tam_popula:
            pai1 = selecao_por_ranqueamento(populacao, aptidoes_invertidas)
            pai2 = selecao_por_ranqueamento(populacao, aptidoes_invertidas)
            while pai2 == pai1:
                pai2 = selecao_por_ranqueamento(populacao, aptidoes_invertidas)
            filho1, filho2 = ox_crossover_tsp(pai1, pai2)
            filho1 = mutacao_swap(filho1, 0.1)
            filho2 = mutacao_swap(filho2, 0.2)
            nova_popula.append(filho1)
            if len(nova_popula) < tam_popula:
                nova_popula.append(filho2)
        populacao = nova_popula

    return populacao, melhor_rota, melhores_p_geracao

def plotar_resultados_tsp():
    distancias = ler_arquivo_tsp("edgesbrasil58.tsp")
    qtde_cidades = 58
    geracoes = 70
    num_execucoes = 30

    dados_execucoes = {
        "Execução": [],
        "Custo": []
    }
    melhores_globais = []

    for i in range(num_execucoes):
        print(f"Execução {i+1}/{num_execucoes}")
        _, _, melhores_p_geracao = algoritmo_genetico(distancias, qtde_cidades=qtde_cidades, num_geracoes=geracoes)
        dados_execucoes["Execução"].extend([f"Exec {i+1}"] * geracoes)
        dados_execucoes["Custo"].extend(melhores_p_geracao)
        melhores_globais.append(min(melhores_p_geracao))

    df_execucoes = pd.DataFrame(dados_execucoes)

    plt.figure(figsize=(12, 6))
    sns.boxplot(x="Execução", y="Custo", data=df_execucoes, showfliers=False, color="lightblue")
    sns.stripplot(x="Execução", y="Custo", data=df_execucoes, jitter=0.25, size=2, alpha=0.5, color='black')
    ymin = df_execucoes["Custo"].min()
    offset = (df_execucoes["Custo"].max() - ymin) * 0.3

    for i, melhor in enumerate(melhores_globais):
        plt.text(i, ymin - offset, f'{melhor}', ha='center', va='top', fontsize=8, color='darkred', rotation=90)

    plt.title("Distribuição dos Melhores Resultados por Execução - TSP Brasil 58", fontsize=14)
    plt.ylabel("Melhor Custo por Geração")
    plt.xlabel("Execução")
    plt.xticks(rotation=90)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

def plotar_rota(pontos, rota_str, dronometros_da_rota):
    coords = [pontos['R']] + [pontos[letra] for letra in rota_str.split()] + [pontos['R']]
    fig, ax = plt.subplots()

    linhas = max(coord[0] for coord in pontos.values()) + 1
    colunas = max(coord[1] for coord in pontos.values()) + 1
    for i in range(linhas):
        ax.axhline(i, color='lightgray', linewidth=0.5)
    for j in range(colunas):
        ax.axvline(j, color='lightgray', linewidth=0.5)
 
    for label, (x, y) in pontos.items():
        ax.plot(y + 0.5, linhas - x - 0.5, 'o', markersize=10)
        ax.text(y + 0.5, linhas - x - 0.5, label, ha='center', va='center', color='white', bbox=dict(facecolor='black', boxstyle='circle'))

    for i in range(len(coords) - 1):
        x1, y1 = coords[i]
        x2, y2 = coords[i + 1]
        
        path = [
            (y1 + 0.5, linhas - x1 - 0.5),
            (y2 + 0.5, linhas - x1 - 0.5),
            (y2 + 0.5, linhas - x2 - 0.5)
        ]
        ax.plot(*zip(*path), color='blue')

    ax.set_aspect('equal')
    ax.set_xticks(range(colunas))
    ax.set_yticks(range(linhas))
    ax.set_xlim(0, colunas)
    ax.set_ylim(0, linhas)
    ax.invert_yaxis()
    plt.grid(True)
    # remove R da rota_str
    rota_str = rota_str.replace('R', '').strip()
    plt.title(f"Melhor rota: {rota_str}\nDronometros: {dronometros_da_rota}", fontsize=12)
    plt.show()

if __name__ == "__main__":
    nome_arquivo_matriz = "arquivo_5.txt"
    nome_arquivo_tsp = "edgesbrasil58.tsp"

    # Descomente a linha abaixo para executar o algoritmo genético com a matriz
    # main_matriz()
    
    # Descomente a linha abaixo para executar o TSP com o arquivo especificado
    #main_tsp("edgesbrasil58.tsp")

    # Executar o algoritmo genético e plotar os resultados do tsp
    #plotar_resultados_tsp()

    # Descomente a linha abaixo para plotar a rota a partir do arquivo de matriz
    rota_str, dronometros_da_rota = main_matriz(nome_arquivo_matriz)
    plotar_rota(encontrar_pontos(formular_matriz(ler_arquivo_matriz(nome_arquivo_matriz))), rota_str, dronometros_da_rota)