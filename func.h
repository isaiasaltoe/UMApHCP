// Constantes
const int POP_SIZE = 50; // Tamanho da população
const int MAX_GEN = 100 ; // Número máximo de gerações
const double TAXA_MUTACAO = 0.01; // Taxa de mutação

// Estruturas de dados
struct idNo {
    double x, y;
};

struct idRota {
    int OR, H1, H2, DS;
    double custo;
};

struct Individuo {
    int* hubs; // Conjunto de hubs
    double fitness; // Valor da função objetivo (custo máximo)
};

struct idSolucao {
    int numNos, p;
    int* hubs;
    idRota** rotas;
    double fo;
};

// Variáveis globais
int numNos = 0;
idNo* nos = nullptr;
double** distancias = nullptr;

// Protótipos de funções
void ler_dados(const char* arq);
void calculo_distancias();
void liberar_individuo(Individuo* ind);
Individuo gerar_individuo(int p);
void gerar_populacao(Individuo* populacao, int pop_size, int p);
Individuo selecao_torneio(Individuo* pop, int pop_size);
Individuo crossover(Individuo pai1, Individuo pai2, int p);
void mutacao(Individuo* ind, int p);
Individuo algoritmo_genetico(int p, int pop_size, int max_gen);
double calcular_custo_maximo(int p, int* hubs, float beta, float alpha, float lambda, idRota** rotas);
void imprimir_solucao(idSolucao solucao);
void liberar_solucao(idSolucao* solucao);
