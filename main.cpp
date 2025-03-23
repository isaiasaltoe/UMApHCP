#include <iostream>
#include <cmath>
#include <float.h>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include <vector>
#include <random>
#include "func.h"



int main() {
    // Ler dados e calcular distâncias
    ler_dados("inst200.txt");
    calculo_distancias();

    // Parâmetros do algoritmo genético
    int p = 50; // Número de hubs
    int pop_size = POP_SIZE; // Tamanho da população
    int max_gen = MAX_GEN; // Número máximo de gerações

    // Executar algoritmo genético
    Individuo melhor = algoritmo_genetico(p, pop_size, max_gen);

    // Exibir a melhor solução encontrada
    idSolucao solucao;
    solucao.numNos = numNos;
    solucao.p = p;
    solucao.hubs = melhor.hubs;
    solucao.rotas = nullptr; // Não estamos usando rotas aqui
    solucao.fo = melhor.fitness;

    imprimir_solucao(solucao);

    // Liberar memória
    liberar_individuo(&melhor);
    for (int i = 0; i < numNos; i++) delete[] distancias[i];
    delete[] distancias;
    delete[] nos;

    return 0;
}

void ler_dados( const char* arq) {
    FILE* f = fopen(arq, "r");
    if (f == NULL) {
        printf("Erro ao abrir o arquivo.\n");
        return;
    }
    fscanf(f, "%d", &numNos);
    printf("Número de nós: %d\n", numNos);

    nos = new idNo[numNos]; 
    for (int i = 0; i < numNos; i++) {
        fscanf(f, "%lf %lf", &nos[i].x, &nos[i].y);
    }
    fclose(f);
}

void calculo_distancias() {
    distancias = new double*[numNos]; 
    for (int i = 0; i < numNos; i++) {
        distancias[i] = new double[numNos]; 
    }

    double dx, dy;
    for (int i = 0; i < numNos; i++) {
        for (int j = 0; j < numNos; j++) {
            dx = nos[i].x - nos[j].x;
            dy = nos[i].y - nos[j].y;
            distancias[i][j] = sqrt(dx * dx + dy * dy);
        }
    }
}

// Função para gerar um indivíduo aleatório
Individuo gerar_individuo(int p) {
    Individuo ind;
    ind.hubs = new int[p];

    std::vector<int> candidatos(numNos);
    for (int i = 0; i < numNos; i++) candidatos[i] = i;

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(candidatos.begin(), candidatos.end(), g);

    for (int i = 0; i < p; i++) {
        ind.hubs[i] = candidatos[i];
    }

    ind.fitness = calcular_custo_maximo(p, ind.hubs, 1.0, 0.75, 1.0, nullptr);

    return ind;
}

// Função para gerar a população inicial
void gerar_populacao(Individuo* populacao, int pop_size, int p) {
    for (int i = 0; i < pop_size; i++) {
        populacao[i] = gerar_individuo(p);
    }
}

// Função de seleção por torneio
Individuo selecao_torneio(Individuo* pop, int pop_size) {
    int idx1 = rand() % pop_size;
    int idx2 = rand() % pop_size;
    return (pop[idx1].fitness < pop[idx2].fitness) ? pop[idx1] : pop[idx2];
}

// Função de crossover (combinação de dois pais)
Individuo crossover(Individuo pai1, Individuo pai2, int p) {
    Individuo filho;
    filho.hubs = new int[p];

    int ponto = rand() % p;
    for (int i = 0; i < ponto; i++) filho.hubs[i] = pai1.hubs[i];
    for (int i = ponto; i < p; i++) filho.hubs[i] = pai2.hubs[i];

    filho.fitness = calcular_custo_maximo(p, filho.hubs, 1.0, 0.75, 1.0, nullptr);
    return filho;
}

// Função de mutação (alteração aleatória de um hub)
void mutacao(Individuo* ind, int p) {
    if ((rand() / (double)RAND_MAX) < TAXA_MUTACAO) {
        int idx = rand() % p;
        int novo_hub;
        do {
            novo_hub = rand() % numNos;
        } while (std::find(ind->hubs, ind->hubs + p, novo_hub) != ind->hubs + p);
        ind->hubs[idx] = novo_hub;
    }
}

// Algoritmo genético
Individuo algoritmo_genetico(int p, int pop_size, int max_gen) {
    Individuo* populacao = new Individuo[pop_size];
    gerar_populacao(populacao, pop_size, p);

    Individuo melhor = populacao[0];

    for (int gen = 0; gen < max_gen; gen++) {
        Individuo* nova_populacao = new Individuo[pop_size];

        for (int i = 0; i < pop_size; i++) {
            Individuo pai1 = selecao_torneio(populacao, pop_size);
            Individuo pai2 = selecao_torneio(populacao, pop_size);
            Individuo filho = crossover(pai1, pai2, p);
            mutacao(&filho, p);
            nova_populacao[i] = filho;

            if (filho.fitness < melhor.fitness) {
                melhor = filho;
            }
        }

        for (int i = 0; i < pop_size; i++) liberar_individuo(&populacao[i]);
        delete[] populacao;

        populacao = nova_populacao;
    }

    for (int i = 0; i < pop_size; i++) liberar_individuo(&populacao[i]);
    delete[] populacao;

    return melhor;
}

// Função para liberar memória de um indivíduo
void liberar_individuo(Individuo* ind) {
    if (ind == nullptr || ind->hubs == nullptr) return;
    delete[] ind->hubs;
    ind->hubs = nullptr;
}

// Função para calcular o custo máximo (fitness)
double calcular_custo_maximo(int p, int* hubs, float beta, float alpha, float lambda, idRota** rotas) {
    double maxCusto = 0.0;
    double custo_hub, custo_1hub, custo_2hubs;
    int origemHub, destinoHub;

    for (int i = 0; i < numNos; i++) {
        for (int j = 0; j < numNos; j++) {
            double menorCusto = DBL_MAX;
            origemHub = -1;
            destinoHub = -1;

            for (int k = 0; k < p; k++) {
                int hub_k = hubs[k];
                custo_hub = beta * distancias[i][hub_k] + lambda * distancias[hub_k][j];
                if (custo_hub < menorCusto) {
                    menorCusto = custo_hub;
                    origemHub = hub_k;
                    destinoHub = hub_k;
                }

                for (int l = 0; l < p; l++) {
                    int hub_l = hubs[l];
                    custo_2hubs = beta * distancias[i][hub_k] + alpha * distancias[hub_k][hub_l] + lambda * distancias[hub_l][j];
                    if (custo_2hubs < menorCusto) {
                        menorCusto = custo_2hubs;
                        origemHub = hub_k;
                        destinoHub = hub_l;
                    }
                }
            }

            if (rotas != nullptr) {
                rotas[i][j] = {i, origemHub, destinoHub, j, menorCusto};
            }

            if (menorCusto > maxCusto) {
                maxCusto = menorCusto;
            }
        }
    }

    return maxCusto;
}

// Função para imprimir a solução
void imprimir_solucao(idSolucao solucao) {
    printf("n: %d p: %d\n", solucao.numNos, solucao.p);
    printf("FO: %.2f\n", solucao.fo);

    printf("HUBS: [");
    for (int i = 0; i < solucao.p; i++) {
        printf("%d", solucao.hubs[i]);
        if (i < solucao.p - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}