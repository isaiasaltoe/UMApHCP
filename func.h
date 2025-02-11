

#define MAX_NOS 200
#define MAX_HUBS 50 



//Estrutura dos Nós
typedef struct nos{
    double x, y;
    int id;
}idNo;



//Estrutura de rota
typedef struct rota{
int OR;  // Nó de origem
int H1;  // HUB 1
int H2;  // HUB 2
int DS; //  Nó de destino
double custo;  // Custo de transporte entre cada par 
 

}idRota;

//Estrutura de solução
typedef struct solucao{
int numNos;
int p; 
double fo;
int* hubs; 
idRota** rotas;          
}idSolucao;



//Dados de entrada
int numNos;
idNo *nos;
double **distancias;

// Declaração de funções 
void ler_dados(char* arq);
void calculo_distancias();
int* selecionar_hubs(int p);
double calcular_custo_maximo(int p, int* hubs, float beta, float alpha, float lambda, idRota** rotas);
idSolucao Construir_Solucao_inicial(int p);
idSolucao clonar_solucao(idSolucao solucao);
void liberar_solucao(idSolucao *solucao);
void gravar_solucao(idSolucao solucao, char* arquivo_saida);
void imprimir_solucao(idSolucao solucao);
void ler_solucao(idSolucao *solucao, const char* arquivo_entrada);