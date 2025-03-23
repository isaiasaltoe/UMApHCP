#include <iostream>
#include <cmath>
#include <float.h>
#include <cstdlib>
#include <cstdio>
#include "func.h"
#include <time.h>
#include <vector>
#include <random>
#include <algorithm>


int main() {
    printf("Iniciando programa...\n");
    srand(time(NULL));

    printf("Lendo dados da instância...\n");
    ler_dados("inst200.txt");
    calculo_distancias();

    printf("Iniciando algoritmo genético...\n");
    int p = 50;
    int pop_size = POP_SIZE;
    int max_gen = MAX_GEN;
    double tempo_limite = 300;

    Individuo melhor = algoritmo_genetico(p, pop_size, max_gen, tempo_limite);

    printf("Execução do algoritmo concluída!\n");
    printf("Melhor FO: %.2f\n", melhor.fitness);
    printf("Hubs selecionados: ");
    for (int i = 0; i < p; i++) {
        printf("%d ", melhor.hubs[i]);
    }
    printf("\n");

    liberar_individuo(&melhor);
    for (int i = 0; i < numNos; i++) delete[] distancias[i];
    delete[] distancias;
    delete[] nos;

    return 0;
}

void ler_dados(const char* arq) {
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

/* 
int* selecionar_hubs(int p) {
    int* hubs = (int*)malloc(p * sizeof(int));
    hubs[0] = 3;
    hubs[1] = 5;
    hubs[2] = 13;
    hubs[3] = 16;
    return hubs;
}
*/

void liberar_individuo(Individuo* ind) {
    delete[] ind->hubs;
}



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

void gerar_populacao(Individuo* populacao, int pop_size, int p) {
    for (int i = 0; i < pop_size; i++) {
        populacao[i] = gerar_individuo(p);
    }
}


Individuo selecao_torneio(Individuo* pop, int pop_size) {
    int idx1 = rand() % pop_size;
    int idx2 = rand() % pop_size;
    return (pop[idx1].fitness < pop[idx2].fitness) ? pop[idx1] : pop[idx2];
}

Individuo crossover(Individuo pai1, Individuo pai2, int p) {
    Individuo filho;
    filho.hubs = new int[p];

    int ponto = rand() % p;
    for (int i = 0; i < ponto; i++) filho.hubs[i] = pai1.hubs[i];
    for (int i = ponto; i < p; i++) filho.hubs[i] = pai2.hubs[i];

    filho.fitness = calcular_custo_maximo(p, filho.hubs, 1.0, 0.75, 1.0, nullptr);
    return filho;
}

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

Individuo algoritmo_genetico(int p, int pop_size, int max_gen, double tempo_limite) {
    clock_t inicio = clock();
    Individuo* populacao = new Individuo[pop_size];
    gerar_populacao(populacao, pop_size, p);

    Individuo melhor = populacao[0];

    for (int gen = 0; gen < max_gen; gen++) {
        // Verificar tempo limite
        double tempo_gasto = (double)(clock() - inicio) / CLOCKS_PER_SEC;
        if (tempo_gasto >= tempo_limite) break;

        Individuo* nova_populacao = new Individuo[pop_size];

        for (int i = 0; i < pop_size / 2; i++) {
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

    // Liberar memória
    for (int i = 0; i < pop_size; i++) liberar_individuo(&populacao[i]);
    delete[] populacao;

    return melhor;
}

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

            rotas[i][j] = {i, origemHub, destinoHub, j, menorCusto};
            if (menorCusto > maxCusto) {
                maxCusto = menorCusto;
            }
        }
    }

    return maxCusto;
}

/* 
idSolucao Construir_Solucao_inicial(int p) {
    clock_t h;
    double tempo;
    const int rep_hubs = 1000;
    const int rep_fo = 1000;

    idSolucao solucao;
    solucao.numNos = numNos;
    solucao.p = p;

    // Medir tempo da seleção de hubs (1 vez)
    h = clock();  
    solucao.hubs = selecionar_hubs(p); // Já usa new internamente
    h = clock() - h;
    tempo = (double)h / CLOCKS_PER_SEC;
    printf("Tempo da seleção de hubs: %.5f segundos\n", tempo);

    // Alocar memória para as rotas usando new
    solucao.rotas = new idRota*[numNos]; // Alocar array de ponteiros
    for (int i = 0; i < numNos; i++) {
        solucao.rotas[i] = new idRota[numNos]; // Alocar cada linha da matriz
    }

    // Medir tempo da função objetivo (1 vez)
    h = clock();
    solucao.fo = calcular_custo_maximo(p, solucao.hubs, 1.0, 0.75, 1.0, solucao.rotas);
    h = clock() - h;
    tempo = (double)h / CLOCKS_PER_SEC;
    printf("Tempo do cálculo da função objetivo: %.5f segundos\n", tempo);

    return solucao;
}
*/


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

   
    printf("OR H1 H2 DS CUSTO\n");
    for (int i = 0; i < solucao.numNos; i++) {
        for (int j = 0; j < solucao.numNos; j++) {
            printf("%d %d %d %d %.2f\n",
                solucao.rotas[i][j].OR,
                solucao.rotas[i][j].H1,
                solucao.rotas[i][j].H2,
                solucao.rotas[i][j].DS,
                solucao.rotas[i][j].custo);
        }
    }
}

void gravar_solucao(idSolucao solucao, const char* arquivo_saida) {
    FILE* f = fopen(arquivo_saida, "w");
    if (f == NULL) {
        printf("Erro ao abrir o arquivo.\n");
        return;
    }
    fprintf(f, "n: %d p: %d\n", solucao.numNos, solucao.p);
    fprintf(f, "FO: %.2f\n", solucao.fo);
    fprintf(f, "HUBS: [");
    for (int i = 0; i < solucao.p; i++) {
        fprintf(f, "%d", solucao.hubs[i]);
        if (i < solucao.p - 1) {
            fprintf(f, ", ");
        }
    }
    fprintf(f, "]\n");
    fprintf(f, "OR H1 H2 DS CUSTO\n");
    for (int i = 0; i < solucao.numNos; i++) {
        for (int j = 0; j < solucao.numNos; j++) {
            fprintf(f, "%d %d %d %d %.2f\n",
                solucao.rotas[i][j].OR,
                solucao.rotas[i][j].H1,
                solucao.rotas[i][j].H2,
                solucao.rotas[i][j].DS,
                solucao.rotas[i][j].custo);
        }
    }
    fclose(f);
}

idSolucao clonar_solucao(idSolucao solucao) {
    idSolucao clone_solucao;

    clone_solucao.fo = solucao.fo;
    clone_solucao.numNos = solucao.numNos;
    clone_solucao.p = solucao.p;

    if (solucao.numNos <= 0) {
        printf("Erro: Número de nós inválido (%d).\n", solucao.numNos);
        exit(1);
    }

    clone_solucao.hubs = new int[solucao.p]; // Substituir malloc por new
    for (int i = 0; i < solucao.p; i++) {
        clone_solucao.hubs[i] = solucao.hubs[i];
    }

    clone_solucao.rotas = new idRota*[solucao.numNos]; // Substituir malloc por new
    for (int i = 0; i < solucao.numNos; i++) {
        clone_solucao.rotas[i] = new idRota[solucao.numNos]; // Substituir malloc por new
        for (int j = 0; j < solucao.numNos; j++) {
            clone_solucao.rotas[i][j] = solucao.rotas[i][j];
        }
    }

    return clone_solucao;
}

void liberar_solucao(idSolucao *solucao) {
    if (solucao == nullptr) return;

    delete[] solucao->hubs; // Substituir free por delete[]

    if (solucao->rotas != nullptr) {
        for (int i = 0; i < solucao->numNos; i++) {
            delete[] solucao->rotas[i]; // Substituir free por delete[]
        }
        delete[] solucao->rotas; // Substituir free por delete[]
    }

    solucao->hubs = nullptr;
    solucao->rotas = nullptr;
}

void ler_solucao(idSolucao *solucao, const char* arquivo_entrada) {
    FILE* f = fopen(arquivo_entrada, "r");
    if (f == nullptr) {
        printf("Erro ao abrir o arquivo: %s\n", arquivo_entrada);
        return;
    }

    fscanf(f, "n: %d p: %d\n", &solucao->numNos, &solucao->p);

    solucao->hubs = new int[solucao->p]; // Substituir malloc por new
    for (int i = 0; i < solucao->p; i++) {
        fscanf(f, "%d", &solucao->hubs[i]);
        if (i < solucao->p - 1) {
            fscanf(f, ", ");
        }
    }
    fscanf(f, "]\n");

    solucao->rotas = new idRota*[solucao->numNos]; // Substituir malloc por new
    for (int i = 0; i < solucao->numNos; i++) {
        solucao->rotas[i] = new idRota[solucao->numNos]; // Substituir malloc por new
    }

    fscanf(f, "OR H1 H2 DS CUSTO\n");
    for (int i = 0; i < solucao->numNos; i++) {
        for (int j = 0; j < solucao->numNos; j++) {
            fscanf(f, "%d %d %d %d %lf\n",
                &solucao->rotas[i][j].OR,
                &solucao->rotas[i][j].H1,
                &solucao->rotas[i][j].H2,
                &solucao->rotas[i][j].DS,
                &solucao->rotas[i][j].custo);
        }
    }

    fclose(f);
}