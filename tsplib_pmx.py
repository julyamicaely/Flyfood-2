import random as rd
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def ler_arquivo(arquivo="edgesbrasil58.tsp"):
    objArq = open(arquivo)
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

#funcao que retorna o custo total do caminho:
def custoCaminho(permutacao, dicDistancias):
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
	
def inicializaPopulacao(tamanho, qtdeCidades):
	import random
	#criando uma lista com "tamanho" permutacoes aleatorias de cidades:
	lista = []
	for i in range(tamanho):
		individuo = list(range(1, qtdeCidades+1))
		random.shuffle(individuo)
		lista.append(individuo)
	return lista

def calculaAptidao(populacao, distancias):
	listaAptidao = []
	for elem in populacao:
		listaAptidao.append(custoCaminho(elem, distancias))
	return listaAptidao

def roleta(populacao, aptidoes):
    soma = sum(aptidoes)
    # Em problemas de minimizacao, menor custo e melhor, entao iinverter
    inv_apt = [1/a for a in aptidoes]
    soma_inv = sum(inv_apt)
    proba = [a/soma_inv for a in inv_apt]

    acumulada = []
    atual = 0
    for p in proba:
        atual += p
        acumulada.append(atual)
    r = rd.random()
    for i, limite in enumerate(acumulada):
        if r < limite:
            return populacao[i]

def selecionar_nova_geracao(populacao, aptidoes, tamanho_nova_geracao):
    nova_geracao = []
    for _ in range(tamanho_nova_geracao):
        individuo = roleta(populacao, aptidoes)
        nova_geracao.append(individuo)
    return nova_geracao

def pmx_crossover(p1, p2):
    size = len(p1)
    c1, c2 = [None]*size, [None]*size
    cx1, cx2 = sorted(rd.sample(range(1, size-1), 2))
    c1[cx1:cx2+1] = p1[cx1:cx2+1]
    c2[cx1:cx2+1] = p2[cx1:cx2+1]

    def completar_filho(filho, pais):
        for i in range(cx1, cx2+1):
            if pais[i] not in filho:
                pos = i
                val = pais[i]
                while True:
                    mapped = filho[pos]
                    pos = pais.index(mapped)
                    if filho[pos] is None:
                        filho[pos] = val
                        break
        for i in range(len(filho)):
            if filho[i] is None:
                filho[i] = pais[i]
        return filho

    return completar_filho(c1, p2), completar_filho(c2, p1)

def mutacao_swap(rota, taxa=0.1):
    nova = rota[:]
    if rd.random() < taxa:
        i, j = rd.sample(range(1, len(nova)-1), 2)
        nova[i], nova[j] = nova[j], nova[i]
    return nova

def algoritmo_genetico(distancias, qtde_cidades, tamanho_pop=40, geracoes=70):
    populacao = inicializaPopulacao(tamanho_pop, qtde_cidades)
    melhores_p_geracao = []
    for g in range(geracoes):
        aptidoes = calculaAptidao(populacao, distancias)
        nova_geracao = []
        melhor = min(aptidoes)
        melhores_p_geracao.append(melhor)

        while len(nova_geracao) < tamanho_pop:
            pai1 = roleta(populacao, aptidoes)
            pai2 = roleta(populacao, aptidoes)
            filho1, filho2 = pmx_crossover(pai1, pai2)
            filho1 = mutacao_swap(filho1, taxa=0.01)
            filho2 = mutacao_swap(filho2, taxa=0.02)
            nova_geracao.extend([filho1, filho2])

        populacao = nova_geracao[:tamanho_pop]
    
    aptidoes = calculaAptidao(populacao, distancias)
    melhor = min(zip(populacao, aptidoes), key=lambda x: x[1])
    return melhor[0], melhor[1], melhores_p_geracao

if __name__ == "__main__":
    distancias = ler_arquivo("edgesbrasil58.tsp")
    
    
    dados_execucoes = {
        "Execução": [],
        "Custo": []
    }

    melhores_globais = []
    geracoes = 70
    num_execucoes = 30
    for i in range(num_execucoes):
        print(f"Execução {i+1}/{num_execucoes}")
        _, _, melhores_p_geracao = algoritmo_genetico(distancias, qtde_cidades=58)
        
        dados_execucoes["Execução"].extend([f"Exec {i+1}"] * geracoes)
        dados_execucoes["Custo"].extend(melhores_p_geracao)

        melhores_globais.append(min(melhores_p_geracao))

    # Criando o DataFrame com todos os dados
    df_execucoes = pd.DataFrame(dados_execucoes)

    # Plotagem
    plt.figure(figsize=(12, 6))
    sns.boxplot(x="Execução", y="Custo", data=df_execucoes, showfliers=False, color="lightblue")
    sns.stripplot(x="Execução", y="Custo", data=df_execucoes, jitter=0.25, size=2, alpha=0.5, color='black')

    # Adiciona os valores dos melhores globais
    for i, melhor in enumerate(melhores_globais):
        plt.text(i, melhor - 20, f'{melhor}', ha='center', va='bottom', fontsize=8, color='darkred', rotation=90)

    plt.title("Distribuição dos Melhores Resultados por Execução -  TSP Brasil 58", fontsize=14)
    plt.ylabel("Melhor Custo por Geração")
    plt.xlabel("Execução")
    plt.xticks(rotation=90)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

