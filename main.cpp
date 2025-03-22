#include <iostream>
#include <cmath>
#include <float.h>
#include <cstdlib>
#include <cstdio>
#include "func.h"
#include <time.h>
#include <algorithm>


int main() {
    //srand(time(NULL));
    ler_dados("inst200.txt");
    calculo_distancias();
    
    // parametro de numero de hubs
    int p = 50;
  

    idSolucao solucao = Construir_Solucao_inicial(p);

    imprimir_solucao(solucao);
    liberar_solucao(&solucao);
     

    for (int i = 0; i < numNos; i++) {
        free(distancias[i]);
    }
    free(distancias);
    free(nos);
     
    return 0;
}

void ler_dados(char* arq) {
    FILE* f = fopen(arq, "r");
    if (f == NULL) {
        printf("Erro ao abrir o arquivo.\n");
        return;
    }
    fscanf(f, "%d", &numNos);
    printf("Número de nós: %d\n", numNos);

    nos = (idNo*)malloc(numNos * sizeof(idNo));
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

int* selecionar_hubs(int p) {
    double* somaDistancias = (double*)malloc(numNos * sizeof(double));
    int* hubs = (int*)malloc(p * sizeof(int));

    for (int i = 0; i < numNos; i++) {
        somaDistancias[i] = 0.0;
        for (int j = 0; j < numNos; j++) {
            somaDistancias[i] += distancias[i][j];
        }
    }

    int* indices = (int*)malloc(numNos * sizeof(int));
    for (int i = 0; i < numNos; i++) {
        indices[i] = i;
    }

    std::sort(indices, indices + numNos, [&somaDistancias](int i, int j) {
        return somaDistancias[i] > somaDistancias[j];
    });

    
    int mid = p / 2;
    for (int k = 0; k < mid; k++) {
        hubs[k] = indices[k];  
    }

    for (int k = mid; k < p; k++) {
        hubs[k] = indices[numNos - 1 - (k - mid)];  
    }

    free(somaDistancias);
    free(indices);

    return hubs;
}


double calcular_custo_maximo(int p, int* hubs, float beta, float alpha, float lambda, idRota** rotas) {
    double maxCusto = 0.0;

    
    for (int i = 0; i < numNos; i++) {
        for (int j = 0; j < numNos; j++) {
            double menorCusto = DBL_MAX;
            int origemHub = -1;
            int destinoHub = -1;

            
            if (i == j) {
                for (int k = 0; k < p; k++) {
                    int hub_k = hubs[k];
                    double custo_hub = beta * distancias[i][hub_k] + lambda * distancias[hub_k][i];
                    if (custo_hub < menorCusto) {
                        menorCusto = custo_hub;
                        origemHub = hub_k;
                        destinoHub = hub_k; 
                    }
                }
            } else {
                
                for (int k = 0; k < p; k++) {
                    int hub_k = hubs[k];
                    double custo_1hub = beta * distancias[i][hub_k] + lambda * distancias[hub_k][j];
                    if (custo_1hub < menorCusto) {
                        menorCusto = custo_1hub;
                        origemHub = hub_k;
                        destinoHub = hub_k;
                    }

                    for (int l = 0; l < p; l++) {
                        int hub_l = hubs[l];
                        double custo_2hubs = beta * distancias[i][hub_k] + 
                                             alpha * distancias[hub_k][hub_l] + 
                                             lambda * distancias[hub_l][j];

                        if (custo_2hubs < menorCusto) {
                            menorCusto = custo_2hubs;
                            origemHub = hub_k;
                            destinoHub = hub_l;
                        }
                    }
                }
            }
             
        
            
            rotas[i][j].OR = i;
            rotas[i][j].DS = j;
            rotas[i][j].H1 = origemHub;
            rotas[i][j].H2 = destinoHub;
            rotas[i][j].custo = menorCusto;

            if (menorCusto > maxCusto) {
                maxCusto = menorCusto;
            }
        }

        
    }
    return maxCusto;
}

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
    solucao.hubs = selecionar_hubs(p);
    h = clock() - h;
    tempo = (double)h / CLOCKS_PER_SEC;
    printf("Tempo da seleção de hubs: %.5f segundos\n", tempo);
     
    
     /*
    //Medir tempo da seleção de hubs (10000 vezes)
    h = clock();
    for (int i = 0; i < rep_hubs; i++) {
        solucao.hubs = selecionar_hubs(p);
    }
    h = clock() - h;
    tempo = (double)h / CLOCKS_PER_SEC;
    printf("Tempo da seleção de hubs executando %d vezes: %.5f segundos\n", rep_hubs, tempo);
     */
    
    
    // Alocar memória para as rotas
    solucao.rotas = (idRota**)malloc(numNos * sizeof(idRota*));
    for (int i = 0; i < numNos; i++) {
        solucao.rotas[i] = (idRota*)malloc(numNos * sizeof(idRota));
    }
    
    // Medir tempo da função objetivo (1 vez)
    h = clock();
    solucao.fo = calcular_custo_maximo(p, solucao.hubs, 1.0, 0.75, 1.0, solucao.rotas);
    h = clock() - h;
    tempo = (double)h / CLOCKS_PER_SEC;
    printf("Tempo do cálculo da função objetivo: %.5f segundos\n", tempo);
     
  /*
    // Medir tempo da função objetivo (1000 vezes)
    h = clock();
    for (int i = 0; i < rep_fo; i++) {
        calcular_custo_maximo(p, solucao.hubs, 1.0, 0.75, 1.0, solucao.rotas);
    }
    h = clock() - h;
    tempo = (double)h / CLOCKS_PER_SEC;
    printf("Tempo do cálculo da função objetivo executando %d vezes: %.5f segundos\n", rep_fo, tempo);
   */


    return solucao;
}



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

void gravar_solucao(idSolucao solucao, char* arquivo_saida) {
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
    fprintf (f, "OR H1 H2 DS CUSTO\n");
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

 
    clone_solucao.hubs = (int*)malloc(solucao.p * sizeof(int));
    if (clone_solucao.hubs == NULL) {
        printf("Erro ao alocar memória para hubs.\n");
        exit(1);
    }
    
    for (int i = 0; i < solucao.p; i++) {
        clone_solucao.hubs[i] = solucao.hubs[i];
    }

    
    clone_solucao.rotas = (idRota**)malloc(solucao.numNos * sizeof(idRota*));
    if (clone_solucao.rotas == NULL) {
        printf("Erro ao alocar memória para rotas.\n");
        exit(1);
    }

    for (int i = 0; i < solucao.numNos; i++) {
        clone_solucao.rotas[i] = (idRota*)malloc(solucao.numNos * sizeof(idRota));
        if (clone_solucao.rotas[i] == NULL) {
            printf("Erro ao alocar memória para rotas[%d].\n", i);
            exit(1);
        }

        for (int j = 0; j < solucao.numNos; j++) {
            clone_solucao.rotas[i][j] = solucao.rotas[i][j];
        }
    }

    return clone_solucao;
}

void liberar_solucao(idSolucao *solucao) {
    if (solucao == NULL) return;  

    free(solucao->hubs);

    if (solucao->rotas != NULL) {
        for (int i = 0; i < solucao->numNos; i++) {
            free(solucao->rotas[i]);
        }
        free(solucao->rotas);
    }

   
    solucao->hubs = NULL;
    solucao->rotas = NULL;
}

void ler_solucao(idSolucao *solucao, const char* arquivo_entrada) {
    FILE* f = fopen(arquivo_entrada, "r");
    if (f == NULL) {
        printf("Erro ao abrir o arquivo: %s\n", arquivo_entrada);
        return;
    }

    
    fscanf(f, "n: %d p: %d\n", &solucao->numNos, &solucao->p);
    
 
    solucao->hubs = (int*)malloc(solucao->p * sizeof(int));
    if (solucao->hubs == NULL) {
        printf("Erro ao alocar memória para hubs.\n");
        fclose(f);
        return;
    }

 
    fscanf(f, "HUBS: [");
    for (int i = 0; i < solucao->p; i++) {
        fscanf(f, "%d", &solucao->hubs[i]);
        if (i < solucao->p - 1) {
            fscanf(f, ", ");
        }
    }
    fscanf(f, "]\n");

    
    solucao->rotas = (idRota**)malloc(solucao->numNos * sizeof(idRota*));
    for (int i = 0; i < solucao->numNos; i++) {
        solucao->rotas[i] = (idRota*)malloc(solucao->numNos * sizeof(idRota));
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