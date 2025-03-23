#define MAX_NOS 200
#define MAX_HUBS 50
#define POP_SIZE 100
#define MAX_GEN 500
#define TAXA_MUTACAO 0.1

// Estrutura dos Nós
typedef struct nos {
    double x, y;
    int id;
} idNo;

// Estrutura de rota
typedef struct rota {
    int OR;  // Nó de origem
    int H1;  // HUB 1
    int H2;  // HUB 2
    int DS;  // Nó de destino
    double custo;  // Custo de transporte entre cada par
} idRota;

// Estrutura de solução
typedef struct solucao {
    int numNos;
    int p;
    double fo;
    int* hubs;
    idRota** rotas;
} idSolucao;

// Estrutura de indivíduo
struct Individuo {
    int* hubs;      // Vetor de hubs
    double fitness; // Valor da função objetivo
};

// Dados de entrada
int numNos;         // Número de nós
idNo* nos;          // Array de nós
double** distancias; // Matriz de distâncias

// Declarações de funções
void ler_dados(const char* arq);
void calculo_distancias();

double calcular_custo_maximo(int p, int* hubs, float beta, float alpha, float lambda, idRota** rotas);
idSolucao Construir_Solucao_inicial(int p);
idSolucao clonar_solucao(idSolucao solucao);
void liberar_solucao(idSolucao* solucao);
void gravar_solucao(idSolucao solucao, const char* arquivo_saida);
void imprimir_solucao(idSolucao solucao);
void ler_solucao(idSolucao* solucao, const char* arquivo_entrada);
void liberar_individuo(Individuo* ind);

// Funções do algoritmo genético
Individuo gerar_individuo(int p);
void gerar_populacao(Individuo* populacao, int pop_size, int p);
Individuo selecao_torneio(Individuo* pop, int pop_size);
Individuo crossover(Individuo pai1, Individuo pai2, int p);
void mutacao(Individuo* ind, int p);
Individuo algoritmo_genetico(int p, int pop_size, int max_gen, double tempo_limite);