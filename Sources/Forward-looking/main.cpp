
#include <iostream>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <list>
#include <set>
#include <limits>
#include <stdexcept>
#include <string>
#include <cmath>
#include <sstream>
#include <sys/time.h>
#include <tuple>
#include <chrono>

#include "InstanciaPorPeriodo.h"
#include "Util.h"
#include "Objeto.h"
#include "Item.h"

#define CPLEX_GET(X) DOUBLE_TO_INT(cplex.getValue(X))
#define CPLEX_GETD(X) cplex.getValue(X)
#define EXPORTAR_MODELO
#define EXECUTABLE_PATH "./"
#define DIR "../../Instances/bromro2/"
#define INSTANCIA "3.txt"

//#define CARREGAR_SOLUCAO_CPLEX


/*
* Adiciona ao catálogo de cada período:
* 1: apenas os itens do "catálogo de entrada";
* 2: os itens do próximo período mais os itens do "catálogo de entrada";
* 3: os itens dos n próximos períodos mais os itens do "catálogo de entrada";
* 4: os itens de todos os próximos períodos mais os itens do "catálogo de entrada";
*
*/
int OPCAO_ITENS_CATALOGO = 1;

bool adicionarRestricaoArea = true;

using namespace std;

typedef std::numeric_limits< double > dbl;

typedef IloArray<IloNumVarArray> NumVarMatrix2D;
typedef IloArray<NumVarMatrix2D> NumVarMatrix3D;
typedef IloArray<IloNumArray>    FloatMatrix;
typedef IloArray<IloExprArray>   ExprMatrix;
typedef unsigned long long ullong;
typedef unsigned long ulong;


struct Ponto{
    double x, y;
    Ponto(double _x=0, double _y=0):x(_x),y(_y){}
    bool operator==(const Ponto& pt) const{
        return (DOUBLE_EQUALS(x, pt.x) && DOUBLE_EQUALS(y, pt.y));
    }
};

typedef vetor(InstanciaPorPeriodo) Instancia;
Instancia instancia;

void leArgumentos(int argc, char** argv);
void leEntrada(string arquivoEntrada);
bool resolveModeloArtigo(int periodo);
bool resolveModeloTese(int periodo);
bool resolveModelo(int ini, int fim);
bool resolveModeloRF(int ini, int fim, int numPeriodosRelaxados);

void setParametros(IloCplex &cplex, bool solucaoInteira = true);

bool atualizaAreaUtilizadaObjeto(Instancia& instancia);

inline bool cabeAlgumItemNoObjeto(Objeto& obj, vetor( vetor(Item) )& itens);
inline bool cabeAlgumItemNoObjeto(Objeto& obj, vetor(Item)& item);

void print(vetor(int)& v, int ini, int fim);
bool is_file_exist(const char *fileName);
void adicionaPontoDeCorte(vetor(Ponto)&, Ponto&);
void metodoPDA();
bool cabeItemCatalogo(vetor(int) &w_bar, vetor(int) &h_bar, Objeto& sobra);
void imprimeDeltas();
bool deltaRepetidos();

//Variáveis (dados) de entrada
int P;
int xi;


vetor( vetor(Objeto) ) objetosSobras;
vetor( vetor( vetor(Item) ) ) listaItensDoObjeto;
vetor( vetor( vetor(Item) ) ) listaItensDoObjetoCpy;
vetor( Item ) sobrasFinais;
vetor( Par ) refSobrasFinais;
vector< vector<bool> > objetoComItem;
vector< tuple<int, int, double, string> > gaps;

//Indica de qual objeto e período a sobra foi originada
//Objetos que podem ser comprados (1...m[s]) referenciam a sí próprio
vetor( vetor(Par) ) refObjPer;

int WMax = 0;
int HMax = 0;
string nomeInstancia;
string nomeArquivoInstancia;
double tempoLimiteCplex;
bool imprimirLogCplex;
bool imprimirPontosDeCorte;
bool preSolveCplex;
bool setLimiteThreadsCplex;
bool imprimirDeltasFinais;
int numThreadsCplex;
bool setLimiteMemoriaCplex;
int maxMemoriaCplex;

vector< TuplaCor > cores;
vector< pair<string, string> > padrao;
vetor(ulong) solucoesPDA;
vetor(double) tempoIteracaoPDA;

vector< vector< tuple<int, int, Delta> > > valoresDelta;

int main(int argc, char** argv) {

    std::cout.imbue(std::locale(std::cout.getloc(), new DecimalSeparator<char>(',')));
    chrono::time_point<std::chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    leArgumentos(argc, argv);

    leEntrada(nomeInstancia);

    metodoPDA();
    
    end = chrono::system_clock::now();
    chrono::duration<double> tempoGasto = end-start;

    cerr << setprecision(5) << fixed << tempoGasto.count() << endl << flush;

    cout << setprecision(5) << fixed << tempoGasto.count() << endl << endl;

    cout << "Tempo Total Gasto: " << setprecision(5) << fixed << tempoGasto.count() << endl << endl;
    cout << flush << endl;
    
    cout << "Solucoes encontradas: " << solucoesPDA.size() << endl;
    cout << "Solucao\tTempo" << endl;
    for (ullong i = 0; i < solucoesPDA.size(); ++i) {
        cout << solucoesPDA.at(i) << "\t" << tempoIteracaoPDA.at(i) << endl;
    }
    cout << endl;
    
    return 0;
}

void leArgumentos(int argc, char **argv) {
    tempoLimiteCplex = 60; //4 horas
    imprimirLogCplex = false;
    OPCAO_ITENS_CATALOGO = 1;
    preSolveCplex = true;
    setLimiteMemoriaCplex = false;
    setLimiteThreadsCplex = false;
    adicionarRestricaoArea = true;
    imprimirDeltasFinais = false;
    xi = 54;

    if (argc == 1) {
        nomeInstancia = string(DIR) + string(INSTANCIA);
        tempoLimiteCplex = 60;
        imprimirLogCplex = false;
        imprimirDeltasFinais = true;
        setLimiteMemoriaCplex = true;
        maxMemoriaCplex = 16000;
        setLimiteThreadsCplex = true;
        numThreadsCplex = 1;
        preSolveCplex = true;
        adicionarRestricaoArea = false;
        OPCAO_ITENS_CATALOGO = 1;
        xi = 4;
    }

    for (int i = 1; i < argc; i++) {
        if (string(argv[i]) == "-i") {
            nomeInstancia = string(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-tl") {
            tempoLimiteCplex = atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-logCplex") {
            imprimirLogCplex = true;
            continue;
        }
        if (string(argv[i]) == "-opc") {
            OPCAO_ITENS_CATALOGO = atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-preSolve") {
            preSolveCplex = true;
            continue;
        }
        if (string(argv[i]) == "-mem") {
            setLimiteMemoriaCplex = true;
            maxMemoriaCplex = atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-threads") {
            setLimiteThreadsCplex = true;
            numThreadsCplex = atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-lssm") {
            NUM_MAX_ITERACAO_SEM_MELHORA = atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-deltaIni") {
            valorInicialDelta = atof(argv[i + 1]);
        }
        if (string(argv[i]) == "-sigma") {
            valorSigma = atof(argv[i + 1]);
        }
        if (string(argv[i]) == "-xi") {
            xi = atoi(argv[i + 1]);
        }
    }
    xi = xi < 0 ? 0 : xi;

    //separando o nome da instancia
    string::size_type loc = nomeInstancia.find_last_of("/", nomeInstancia.size());
    string::size_type loc2 = nomeInstancia.find_last_of(".", nomeInstancia.size());
    if (!((loc != string::npos) && (loc2 != string::npos))) {
        cout << "ERRO!" << endl;
        exit(0);
    }
    nomeArquivoInstancia = string("");
    nomeArquivoInstancia.append(nomeInstancia, loc + 1, loc2 - loc - 1);

    cout << "\n\nInstancia: " << nomeInstancia << endl;
    cout << "xi: " << xi << endl;
    cout << "Metodo: Forward-Looking Matheuristic Approach" << endl;
    cout << "\t -Limite solucoes sem melhora: " << NUM_MAX_ITERACAO_SEM_MELHORA << endl;
    cout << "\t -Valor inicial dos deltas: " << valorInicialDelta << endl;
    cout << "\t -Valor sigma: " << valorSigma << endl;
    cout << "\t\t Delta_novo <- Delta_velho + (Delta_novo - Delta_velho)*" << valorSigma << "^iteracao\n";
#ifdef AREANULA_DELTAZERO
    cout << "\t -Area da sobra nula -> Delta da sobra recebe zero" << endl;
#else
    cout << "\t -Area da sobra nula -> Mantem o delta da iteracao anterior" << endl;
#endif
#ifndef PODE_ZERAR_DELTAS
    cout << "\t -Deltas dos objetos selecionados e nao utilizados nao sao alterados" << endl;
#else
    cout << "\t -Deltas dos objetos selecionados e nao utilizados sao zerados" << endl;
#endif
    
    cout << "Tempo Limite Cplex: " << tempoLimiteCplex << " (por iteracao)" << endl;
    if (setLimiteMemoriaCplex)
        cout << "Limite de memoria Cplex: " << maxMemoriaCplex << " Mb" << endl;
    if (setLimiteThreadsCplex)
        cout << "Limite de threads Cplex: " << numThreadsCplex << endl;
    cout << "Pre-solve? " << (preSolveCplex ? "Ativado" : "Desativado") << endl;
    cout << "Imprimir Log Cplex? " << std::boolalpha << imprimirLogCplex << endl;
    cout << "Opcao Itens do Catalogo: " << OPCAO_ITENS_CATALOGO << endl;
}

void metodoPDA(){
    chrono::time_point<std::chrono::system_clock> startItPDA, endItPDA;
    int s, j;
    int iteracaoMelhorSol = 0;
    int solucoesRepetidasSeguidas = 1;
    bool algumDeltaMudou = true;
    long long custoMelhorSolucao = INF;
    long long maiorValorSobrasRemanescentes = INF;
    bool encontrouSolucaoViavel;
    long long PESO = 0;
    for (s = 1; s <= P; s++) {
        for (j = 1; j <= instancia[s].m; ++j) {
            PESO += instancia[s].objeto[j].getCusto();
        }
    }
    vector< tuple<int, int, Delta> > deltas;
    for(s=1; s < P; s++) { 
        int numObj = instancia[s].m_bar;
        for (j = 1; j <= numObj/*instancia[s].m_bar*/; j++) {
            deltas.emplace_back(make_tuple(s, j, instancia[s].objeto[j].delta));
        }
    }

    valoresDelta.emplace_back(deltas);

    while(true) {
        startItPDA = chrono::system_clock::now();
        cout << "******************************************** Iteracao " << iteracao << " ********************************************" << endl;
        bool modeloViavel;
        encontrouSolucaoViavel = true;
        listaItensDoObjeto.clear();
        listaItensDoObjeto.resize((unsigned long) (P + 1));
        for (s = 1; s <= P; s++) {
            listaItensDoObjeto[s].resize((unsigned long) (instancia[s].m + 1));
            modeloViavel = resolveModeloArtigo(s);

            if(!modeloViavel) {
                encontrouSolucaoViavel = false;
            }
        }

        long long custoTotal = 0;
        long long valorSobrasRemanescentes = 0;
        for (s = 1; s <= P; s++) {
            for (j = 1; j <= instancia[s].m; ++j) {
                if (instancia[s].objeto[j].utilizado) {
                    custoTotal += instancia[s].objeto[j].getCusto();
                    valorSobrasRemanescentes += instancia[s].objeto[j].getValorSobraRemanescente();
                }
            }
        }

        if(encontrouSolucaoViavel) solucoesPDA.push_back(PESO*custoTotal - valorSobrasRemanescentes);

        if(encontrouSolucaoViavel && custoMelhorSolucao > custoTotal) {
            custoMelhorSolucao = custoTotal;
            listaItensDoObjetoCpy = listaItensDoObjeto;
            solucoesRepetidasSeguidas = 1;
            maiorValorSobrasRemanescentes = valorSobrasRemanescentes;
            iteracaoMelhorSol = iteracao;
        }
        else if(encontrouSolucaoViavel && custoMelhorSolucao == custoTotal && maiorValorSobrasRemanescentes < valorSobrasRemanescentes){
            custoMelhorSolucao = custoTotal;
            listaItensDoObjetoCpy = listaItensDoObjeto;
            solucoesRepetidasSeguidas = 1;
            maiorValorSobrasRemanescentes = valorSobrasRemanescentes;
            iteracaoMelhorSol = iteracao;
        }
        else solucoesRepetidasSeguidas++;

        algumDeltaMudou = atualizaAreaUtilizadaObjeto(instancia);

#ifdef IMPRIMIR_INFO
        cout << "Valores dos deltas" << endl;
        for (s = 1; s <= P; s++) {
            cout << "Periodo " << s << endl;
            for (j = 1; j <= instancia[s].m; ++j) {
                if(instancia[s].objeto[j].utilizado)
                    cout << setprecision(2) << fixed
                         << instancia[s].objeto[j].id << " " << instancia[s].objeto[j].codigo << " " << instancia[s].objeto[j].codigoObjetoOrigem
                         << " " << instancia[s].objeto[j].w << " x " << instancia[s].objeto[j].h << " "
                         << " O(" << instancia[s].objeto[j].getArea() << " " << (int)instancia[s].objeto[j].areaUtilizada << " "
                         << instancia[s].objeto[j].delta.objeto << ") S("
                         << (int)instancia[s].objeto[j].areaSobraDireita+(int)instancia[s].objeto[j].areaSobraSuperior
                         << " " << (int)instancia[s].objeto[j].areaSobraDireitaUtilizada+(int)instancia[s].objeto[j].areaSobraSuperiorUtilizada << " "
                         << instancia[s].objeto[j].delta.sobras << ") D("
                         << (int)instancia[s].objeto[j].areaSobraDireita << " " << (int)instancia[s].objeto[j].areaSobraDireitaUtilizada << " "
                         << instancia[s].objeto[j].delta.sobraDireita << ") T("
                         << (int)instancia[s].objeto[j].areaSobraSuperior << " " << (int)instancia[s].objeto[j].areaSobraSuperiorUtilizada << " "
                         << instancia[s].objeto[j].delta.sobraSuperior << ")" << endl;
            }
            cout << endl;
        }
#endif

        for (s = 1; s <= P; s++)
            instancia[s].resetaInstancia();

        endItPDA = chrono::system_clock::now();
        chrono::duration<double> tempoGastoItPDA = endItPDA-startItPDA;
        cout << "\tFO: " << PESO*custoTotal-valorSobrasRemanescentes << (encontrouSolucaoViavel?" (Valida)":" (Invalida)") << "\n";
        cout << "\tCusto Objetos: " << custoTotal << "\n\tValor sobras: " << valorSobrasRemanescentes << "\n";
        cout << "\tTempo gasto na iteracao: " << tempoGastoItPDA.count() << "\n" << endl;

        if(encontrouSolucaoViavel) tempoIteracaoPDA.push_back(tempoGastoItPDA.count());

        if(solucoesRepetidasSeguidas > NUM_MAX_ITERACAO_SEM_MELHORA) {
            cout << "NUMERO MAXIMO DE ITERACAO SEM MELHORA ATINGIDO" << endl;
            break;
        }

        if(!algumDeltaMudou) break;

        iteracao++;
    }

    if(imprimirDeltasFinais) imprimeDeltas();

    cout << endl;
    cerr << fixed << setprecision(5);
    if(custoMelhorSolucao != INF) {
        long long FO = PESO*custoMelhorSolucao - maiorValorSobrasRemanescentes;
        cout << endl << "Melhor Solucao: " << FO << endl;
        cout << "Custo Objetos: " << custoMelhorSolucao << "\tValor Sobras: " << maiorValorSobrasRemanescentes << endl << endl;
        cout << nomeArquivoInstancia << "\t" << xi << "\t"
             << valorSigma << "\t" << valorInicialDelta << "\t"
             << iteracao << "\t" << iteracaoMelhorSol << "\t" << FO << "\t" << custoMelhorSolucao << "\t" << maiorValorSobrasRemanescentes << "\t";

        cerr << nomeArquivoInstancia << "\t" << xi << "\t"
             << valorSigma << "\t" << valorInicialDelta << "\t"
             << iteracao << "\t" << iteracaoMelhorSol << "\t" << FO << "\t" << custoMelhorSolucao << "\t" << maiorValorSobrasRemanescentes << "\t";

        listaItensDoObjeto = listaItensDoObjetoCpy;
    }
    else {
        cout << endl << "***Solucao não encontrada***" << endl;
        cerr << "Solução não encontrada" << "\t";
    }

}

bool resolveModeloArtigo(int periodo) {
    double inicio = MyClock();
    IloEnv env;
    ulong i, j, k, l, s, j1, j2, il;
    //int m = instancia[periodo].m;
    ulong numObjetos = 0;
    ulong numSobrasValidas = 0;
    ulong m = 0;
    ulong m_bar = instancia[periodo].m_bar;
    ulong n = instancia[periodo].n;
    ulong p = instancia[periodo].p;
    ulong d = instancia[periodo].d;
    vetor(int) &n_q = instancia[periodo].n_q;
    vetor(int) &o_q = instancia[periodo].o_q;
    vetor(int) &w_bar = instancia[periodo].w_bar;
    vetor(int) &h_bar = instancia[periodo].h_bar;
    vetor(Item) &item = instancia[periodo].item;
#ifdef IMPRIMIR_INFO
    cout << "--- Periodo: " << periodo << endl;
    cout << "Tempo Limite: " << tempoLimiteCplex << endl;
#endif
    vetor(Objeto) objeto;
    objeto.push_back(Objeto());
    HMax = 0;
    WMax = 0;
#ifdef IMPRIMIR_INFO
    cout << "\tObjetos:\n";
#endif
    //Copia apenas os objetos que possuiem dimensões válidas e que podem conter ao menos um item do período
    for(j=1; j<=m_bar; ++j){
        if(instancia[periodo].objeto[j].ehObjetoValido() && periodo - instancia[periodo].objeto[j].id.periodo <= xi){
            Objeto obj = instancia[periodo].objeto[j];
#ifdef IMPRIMIR_INFO
            cout << "\t\t";
            cout << obj.id << " ";
            cout << obj.w << " x " << obj.h
                 << " " << (obj.ehSobra?"Sobra":"Objeto");
#endif
            if(cabeAlgumItemNoObjeto(obj, item)) {
#ifdef IMPRIMIR_INFO
                cout << " (Incluido)\n";
#endif
                objeto.push_back(obj);
                m++;

                WMax = (obj.w>WMax?obj.w:WMax);
                HMax = (obj.h>HMax?obj.h:HMax);
                if(obj.ehSobra)
                    numSobrasValidas++;
                else
                    numObjetos++;
            }
            else{
#ifdef IMPRIMIR_INFO
                cout << " (Nao incluido)\n";
#endif
                if(periodo < P && periodo - obj.id.periodo < xi){
                    if(obj.id.periodo != periodo) {
                        instancia[periodo + 1].adicionaObjetoSobraDireita(obj);
                        listaItensDoObjeto[obj.id.periodo][obj.id.objeto].push_back(Item(obj.w, obj.h, obj.x, obj.y, periodo, true));
                    }
                }
                else{
                    if((obj.id.periodo != P)) { //Objeto originado de sobra mas não utilizado
                        listaItensDoObjeto[obj.id.periodo][obj.id.objeto].push_back(Item(obj.w, obj.h, obj.x, obj.y, periodo, true));
                    }
                }
            }

        }
    }

#ifdef IMPRIMIR_INFO
    cout << "\tItens: " << endl;
    for(i = 1; i <= n; ++i)
        cout << "\t\t" << item[i].w << " x " << item[i].h << endl;
    cout << "\tItens catálogo: " << endl;
    for(i = 1; i <= d; ++i)
        cout << "\t\t" << w_bar[i] << " x " << h_bar[i] << endl;
    cout << endl;
    cout << "WMax = " << WMax << "\n" << "HMax = " << HMax << endl;
#endif

    try {
        IloModel modelo(env, "Cutting");
        IloCplex cplex(env);

        m_bar = 2*m + 1;

        auto L = (unsigned int) (floor(log2(WMax)) + 1);
        vector<int> pot(L+1);
        for(k=1; k<=L; ++k)
            pot[k] = (1<<(k-1));

        //Variáveis de decisão
        IloNumVarArray u(env, m + 1, 0, 1, ILOBOOL);
        NumVarMatrix2D v(env, n + 1);
        IloNumVarArray x(env, n + 1, 0.0, IloInfinity, ILOFLOAT);
        IloNumVarArray y(env, n + 1, 0.0, IloInfinity, ILOFLOAT);
        NumVarMatrix2D pi(env, n + 1);
        NumVarMatrix2D tau(env, n + 1);
        IloNumVarArray r(env, m + 1, 0, IloInfinity, ILOFLOAT);
        IloNumVarArray t(env, m + 1, 0, IloInfinity, ILOFLOAT);
        NumVarMatrix2D W_bar(env, 2);
        NumVarMatrix2D H_bar(env, 2);
        IloNumVarArray eta(env, m + 1, 0, 1, ILOBOOL);

        NumVarMatrix2D theta(env, m_bar+1);
        IloNumVarArray gamma(env, m_bar+1, 0.0, IloInfinity, ILOFLOAT);
        NumVarMatrix2D mu(env, m_bar+1);
        NumVarMatrix2D q(env, m_bar+1);
        NumVarMatrix2D omega(env, m_bar+1);
        NumVarMatrix2D zeta(env, m_bar+1);

        IloNumVarArray piso_custoDelta(env, m + 1, 0.0, IloInfinity, ILOINT);

        for (i = 1; i <= n; ++i) {
            v[i] = IloNumVarArray(env, m + 1, 0, 1, ILOBOOL);
            pi[i] = IloNumVarArray(env, n + 1, 0, 1, ILOBOOL);
            tau[i] = IloNumVarArray(env, n + 1, 0, 1, ILOBOOL);
        }
        for(j = 1; j <= m_bar; ++j){
            theta[j] = IloNumVarArray(env, L+1, 0, 1, ILOBOOL);
            omega[j] = IloNumVarArray(env, L+1, 0.0, IloInfinity, ILOFLOAT);
            zeta[j] = IloNumVarArray(env, d+1, 0, 1, ILOBOOL);
            mu[j] = IloNumVarArray(env, d+1, -IloInfinity, IloInfinity, ILOFLOAT);
            q[j] = IloNumVarArray(env, d+1, 0, 1, ILOBOOL);
        }
        W_bar[0] = IloNumVarArray(env, m+1, 0.0, IloInfinity, ILOFLOAT);
        H_bar[0] = IloNumVarArray(env, m+1, 0.0, IloInfinity, ILOFLOAT);
        W_bar[1] = IloNumVarArray(env, m_bar+1, 0.0, IloInfinity, ILOFLOAT);
        H_bar[1] = IloNumVarArray(env, m_bar+1, 0.0, IloInfinity, ILOFLOAT);

        //Número de objeto utilizados
        IloExpr objetosUsados(env);
        //Custo dos objetos utilizados
        IloExpr custoObjetos(env);
        //Valor da proporção de utilização das sobras dos objetos
        IloExpr deltasObjetos(env);
        //Custo dos objetos de entrada utilizados
        IloExpr custoObjetosEntrada(env);
        //Valor das sobras dos objetos utilizados
        IloExpr valorSobras(env);
        //Valor das sobras finais (do perído P + 1)
        IloExpr valorSobrasRemanescente(env);
        //Valor das sobras dos objetos usados mais as sobras provenientes de outros períodos não utilizadas
        IloExpr totalValorSobras(env);
        //IloExpr leftoversValue2(env);
        IloExpr numeroSobras(env);
        IloExpr objetivoOriginal(env);
        IloExpr soma_piso_custoDelta(env);
        long double PESO = 0.0;

        if (periodo != P) {
            PESO = 0.0;
            for (j = 1; j <= m; ++j) {
                PESO += objeto[j].getCustoSobra();
            }
            for (j = 1; j <= m; ++j) {
                j1 = 2 * j - 1;
                j2 = 2 * j;
                objetosUsados += u[j];
                custoObjetos += (objeto[j].getCusto() * u[j]);
                deltasObjetos += (objeto[j].getCustoSobraPorUnidadeArea() *
                                  (objeto[j].delta.sobraSuperior * gamma[j1] +    /************Verificar se eh sobra superior ***/
                                   objeto[j].delta.sobraDireita * gamma[j2]));
                soma_piso_custoDelta += piso_custoDelta[j];
                if(!objeto[j].ehSobra) {
                    custoObjetosEntrada += (objeto[j].getCusto() * u[j]);
                }

                //Iniciando os deltas com 1.0, todos os objetos são selecionados na primeira iteração, mesmo
                //que ele não contenha itens. Sem essa parte, o objeto pode ou não ser selecionado.
                valorSobras += (objeto[j].getCustoSobraPorUnidadeArea()*(gamma[j1] + gamma[j2]));

                //Calcula o valor das sobras
                totalValorSobras += (objeto[j].getCustoSobraPorUnidadeArea() * (gamma[j1] + gamma[j2]));
                if(objeto[j].ehSobra) { //Se o objeto é uma sobra de outro período e não foi utilizada é contabilizada
                    totalValorSobras += (objeto[j].getCustoSobraPorUnidadeArea()*objeto[j].getArea()*(1-u[j]));
                }
            }
            //objetivoOriginal = (PESO*(custoObjetos-deltasObjetos) - valorSobras);
            objetivoOriginal = (PESO*(custoObjetos - soma_piso_custoDelta) - valorSobras);
        }
        else {
            PESO = 0;
            for (j = 1; j <= m; ++j) {
                j1 = 2 * j - 1;
                j2 = 2 * j;
                if (!objeto[j].ehSobra) {
                    custoObjetos += (objeto[j].getCusto() * u[j]);
                }
                objetosUsados += u[j];
//                valorSobras += (objeto[j].getCustoSobraPorUnidadeArea() * (gamma[j1] + gamma[j2]));
                if (periodo - objeto[j].id.periodo < xi){
                    valorSobras += (objeto[j].getCustoSobraPorUnidadeArea() * (gamma[j1] + gamma[j2]));
                    valorSobrasRemanescente += (objeto[j].getCustoSobraPorUnidadeArea() * (gamma[j1] + gamma[j2]));
                }
                PESO += objeto[j].getCustoSobra();
            }
            objetivoOriginal = PESO * custoObjetos - valorSobras;
        }

        //=====================================================Funcao Objetivo
        modelo.add(IloMinimize(env, objetivoOriginal));

        //=====================================================Restricoes
        IloRangeArray restricoes(env);

        //
#ifdef LIMITAR_CUSTO_OBJETOS
        restricoes.add(custoObjetosEntrada - limiteCustoObjetos <= 0);
#endif
        //=====================================================
        for (j = 1; j <= m; ++j){
            j1 = 2 * j - 1;
            j2 = 2 * j;
            restricoes.add(piso_custoDelta[j] - (objeto[j].getCustoSobraPorUnidadeArea() * (objeto[j].delta.sobraSuperior * gamma[j1] + objeto[j].delta.sobraDireita * gamma[j2]) ) <= 0 );
        }

        //=====================================================
        /*
        forall(i in ITENS)
            sum(j in OBJETOS) v[i][j] == 1;
        */
        for (i = 1; i <= n; ++i) {
            IloExpr somatorio(env);
            for (j = 1; j <= m; ++j)
                somatorio += v[i][j];
            restricoes.add(somatorio == 1);
            somatorio.end();
        }

        //=====================================================
        /*
        forall(i in ITENS, j in OBJETOS)
            u[j] >= v[i][j];
        */
        for (i = 1; i <= n; ++i)
            for (j = 1; j <= m; ++j)
                restricoes.add(u[j] - v[i][j] >= 0);

        //=====================================================
        for (j = 1; j <= m; ++j) {
            IloExpr somatorio(env);
            for (i = 1; i <= n; ++i)
                somatorio += v[i][j];
            restricoes.add(u[j] - somatorio <= 0);
            somatorio.end();
        }

        //Atribui o tamanho de cada objeto
        for (j = 1; j <= m; ++j) {
            restricoes.add(W_bar[0][j] - objeto[j].w == 0);
            restricoes.add(H_bar[0][j] - objeto[j].h == 0);
        }

        //=====================================================
        for (j = 1; j <= m; ++j) {
            restricoes.add(t[j] - H_bar[0][j] <= 0);
            restricoes.add(r[j] - W_bar[0][j] <= 0);
        }

        //=====================================================
        for(j = 1; j <= m; ++j){
            j1 = 2 * j - 1;
            j2 = 2 * j;
            restricoes.add(H_bar[1][j1] - (u[j] * HMax) <= 0);
            restricoes.add(t[j] - (1 - u[j]) * HMax - H_bar[1][j1] <= 0);
            restricoes.add(H_bar[1][j1] - (t[j] + (1 - u[j]) * HMax) <= 0);

            restricoes.add(W_bar[1][j1] - (u[j] * WMax) <= 0);
            restricoes.add(W_bar[0][j] - r[j] - (1 - eta[j]) * WMax - (1 - u[j]) * WMax - (W_bar[1][j1]) <= 0);
            restricoes.add(W_bar[1][j1] - (W_bar[0][j] - r[j] + (1 - eta[j]) * WMax + (1 - u[j]) * WMax) <= 0);
            restricoes.add(W_bar[0][j] - eta[j] * WMax - (1 - u[j]) * WMax - (W_bar[1][j1]) <= 0);
            restricoes.add(W_bar[1][j1] - (W_bar[0][j] + eta[j] * WMax + (1 - u[j]) * WMax) <= 0);


            restricoes.add(W_bar[1][j2] - (u[j] * WMax) <= 0);
            restricoes.add(r[j] - (1 - u[j]) * WMax - (W_bar[1][j2]) <= 0);
            restricoes.add(W_bar[1][j2] - (r[j] + (1 - u[j]) * WMax) <= 0);

            restricoes.add(H_bar[1][j2] - (u[j] * HMax) <= 0);
            restricoes.add(H_bar[0][j] - (1 - eta[j]) * HMax - (1 - u[j]) * HMax - (H_bar[1][j2]) <= 0);
            restricoes.add(H_bar[1][j2] - (H_bar[0][j] + (1 - eta[j]) * HMax + (1 - u[j]) * HMax) <= 0);
            restricoes.add(H_bar[0][j] - t[j] - eta[j] * HMax - (1 - u[j]) * HMax - (H_bar[1][j2]) <= 0);
            restricoes.add(H_bar[1][j2] - (H_bar[0][j] - t[j] + eta[j] * HMax + (1 - u[j]) * HMax) <= 0);

        }

        //=====================================================
        for (i = 1; i <= n; ++i) {
            restricoes.add(x[i] - item[i].w / 2.0 >= 0);
            restricoes.add(y[i] - item[i].h / 2.0 >= 0);
        }
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= m; ++j) {
                restricoes.add(x[i] + item[i].w / 2.0 - (W_bar[0][j] - r[j] + (WMax) * (1 - v[i][j])) <= 0);
                restricoes.add(y[i] + item[i].h / 2.0 - (H_bar[0][j] - t[j] + (HMax) * (1 - v[i][j])) <= 0);
            }
        }

        //=====================================================
        for(j = 1; j <= m; ++j)
            for(i = 1; i <= n; ++i)
                for(il = i+1; il <= n; ++il){
                    restricoes.add( x[i] - x[il]  - (0.5*(item[i].w + item[il].w) - WMax*((1-v[i][j]) + (1-v[il][j]) + pi[i][il] + tau[i][il])) >= 0);
                    restricoes.add(-x[i] + x[il]  - (0.5*(item[i].w + item[il].w) - WMax*((1-v[i][j]) + (1-v[il][j]) + pi[i][il] + (1 - tau[i][il]))) >= 0);

                    restricoes.add( y[i] - y[il]  - (0.5*(item[i].h + item[il].h) - HMax*((1-v[i][j]) + (1-v[il][j]) + (1 - pi[i][il]) + tau[i][il])) >= 0);
                    restricoes.add(-y[i] + y[il]  - (0.5*(item[i].h + item[il].h) - HMax*((1-v[i][j]) + (1-v[il][j]) + (1 - pi[i][il]) + (1 - tau[i][il]))) >= 0);
                }

        //=====================================================
        // Restrição de simetria de itens iguais
        for (ulong q = 1; q <= p; ++q)
            for (i = o_q[q] + 1; i <= o_q[q] + n_q[q]; ++i)
                for (ulong i1 = i + 1; i1 <= o_q[q] + n_q[q]; ++i1)
                    for (j = 1; j <= m; ++j) {
                        restricoes.add(x[i1] - item[i1].w / 2.0 - (x[i] + item[i].w / 2.0 - WMax * (1 - v[i][j]) - WMax * (1 - v[i1][j]) - WMax * pi[i][i1]) >= 0);
                        restricoes.add(y[i1] - item[i1].h / 2.0 - (y[i] + item[i].h / 2.0 - HMax * (1 - v[i][j]) - HMax * (1 - v[i1][j]) - HMax * (1 - pi[i][i1])) >= 0);
                    }
        for (ulong q = 1; q <= p; ++q)
            for (i = o_q[q] + 1; i <= o_q[q] + n_q[q]; ++i)
                for (ulong i1 = o_q[q] + n_q[q] + 1; i1 <= n; ++i1)
                    for (j = 1; j <= m; ++j) {
                        restricoes.add(x[i1] + item[i1].w / 2.0 - (x[i] - item[i].w / 2.0 + WMax * (1 - v[i][j]) + WMax * (1 - v[i1][j]) + WMax * pi[i][i1] + WMax * tau[i][i1]) <= 0);
                        restricoes.add(x[i1] - item[i1].w / 2.0 - (x[i] + item[i].w / 2.0 - WMax * (1 - v[i][j]) - WMax * (1 - v[i1][j]) - WMax * pi[i][i1] - WMax * (1 - tau[i][i1])) >= 0);
                        restricoes.add(y[i1] + item[i1].h / 2.0 - (y[i] - item[i].h / 2.0 + HMax * (1 - v[i][j]) + HMax * (1 - v[i1][j]) + HMax * (1 - pi[i][i1]) + HMax * tau[i][i1]) <= 0);
                        restricoes.add(y[i1] - item[i1].h / 2.0 - (y[i] + item[i].h / 2.0 - HMax * (1 - v[i][j]) - HMax * (1 - v[i1][j]) - HMax * (1 - pi[i][i1]) - HMax * (1 - tau[i][i1])) >= 0);
                    }


        //=====================================================
        for(j = 1; j <= m_bar; ++j){
            IloExpr somatorio(env);
            for(l = 1; l <= L; ++l){
                somatorio += pow(2, l-1)*theta[j][l];
            }
            restricoes.add(W_bar[1][j] - somatorio == 0);
        }

        //=====================================================
        for(j = 1; j <= m_bar; ++j) {
            for (l = 1; l <= L; ++l) {
                restricoes.add(omega[j][l] - H_bar[1][j] <= 0);
                restricoes.add(H_bar[1][j] - (1 - theta[j][l]) * HMax - omega[j][l] <= 0);
                restricoes.add(omega[j][l] - (theta[j][l] * HMax) <= 0);
            }
        }

        //=====================================================
        for(j = 1; j <= m_bar; ++j) {
            for (i = 1; i <= d; ++i) {
                restricoes.add(w_bar[i] - (W_bar[1][j] + WMax * (1 - zeta[j][i])) <= 0);
                restricoes.add(h_bar[i] - (H_bar[1][j] + HMax * (1 - zeta[j][i])) <= 0);
            }
        }

        //=====================================================
        for(j = 1; j <= m_bar; ++j) {
            IloExpr somatorio(env);
            for (l = 1; l <= L; ++l) {
                somatorio += pow(2, l - 1) * omega[j][l];
            }
            restricoes.add(gamma[j] - somatorio <= 0);
        }
        for(j = 1; j <= m_bar; ++j) {
            IloExpr somatorio(env);
            for (i = 1; i <= d; ++i) {
                somatorio += zeta[j][i];
            }
            restricoes.add(gamma[j] - somatorio * WMax * HMax <= 0);
        }

        modelo.add(restricoes);
        cplex.extract(modelo);

        if (!imprimirLogCplex)
            cplex.setOut(env.getNullStream());

        setParametros(cplex, periodo == P);

        cplex.solve();

        bool finished;
#ifdef IMPRIMIR_INFO
        cout << "[###] Solver status: " << cplex.getStatus() << endl;
#endif
        if (cplex.getStatus() == IloAlgorithm::Optimal) {
            finished = true;
        } else if (cplex.getStatus() == IloAlgorithm::Feasible) {
            finished = true;
        } else {
#ifdef IMPRIMIR_INFO
            cerr << "[###] Solver status: " << cplex.getStatus() << endl;
#endif
            finished = false;
        }
        if(!finished) {
            env.end();
            return finished;
        }

        IloInt accumulatedIterations = cplex.getNiterations();
        IloInt accumulatedNodes = cplex.getNnodes();
        IloNum fObjBestBound = cplex.getBestObjValue();
        IloNum fObjValue = cplex.getObjValue();
        double cplexGap = cplex.getMIPRelativeGap();//(fObjValue - fObjBestBound) / (0.00000001 + fObjValue);

#ifdef IMPRIMIR_INFO
        cout.precision(dbl::digits10);
        cout << setprecision(2) << fixed;
        cout << "[###] Best objective bound: " << fObjBestBound << endl;
        cout << "[###] Cplex GAP: " << fabs(cplexGap) * 100.0 << "%" << endl;
        cout << "[###] Objective obtained: " << fObjValue << endl;
//        cout << "[###] MIP Iterations: " << accumulatedIterations << endl;
//        cout << "[###] B&B Nodes: " << accumulatedNodes << endl << endl;

//        cout << "----------------SOLUTION FOUND----------------\n";
       // cout << "Objects cost: " << CPLEX_GETD(custoObjetos) << endl;
        cout << "Custo dos objetos: " << (periodo!=P?CPLEX_GETD(custoObjetosEntrada):CPLEX_GETD(custoObjetos)) << endl;
        cout << "Valor das sobras: " << (periodo!=P?CPLEX_GETD(totalValorSobras):CPLEX_GETD(valorSobras)) << endl;
        //cout << "Leftovers value: " << CPLEX_GETD(valorSobras) << endl;
        cout << "Objects used: " << CPLEX_GETD(objetosUsados) << endl;
        cout << "Leftovers: " << CPLEX_GETD(numeroSobras) << endl;
        cout << "Cplex time: " << (MyClock() - inicio) / CLOCKS_PER_SEC << endl;
#endif
        long long cte = 0.0;
        long long FO;
        for (j = 1; j <= m; ++j) cte += objeto[j].getCustoSobra();

        if (periodo != P)
            FO = cte * CPLEX_GET(custoObjetosEntrada) - CPLEX_GET(totalValorSobras);
        else
            FO = cte * CPLEX_GET(custoObjetos) - CPLEX_GET(valorSobras);

        cout << setw(0) << setprecision(4) << fixed;
        cout << setw(3) << periodo;
        cout << setw(17) << (periodo != P ? CPLEX_GET(custoObjetosEntrada) : CPLEX_GET(custoObjetos))
             << setw(16) << (periodo != P ? CPLEX_GET(totalValorSobras) : CPLEX_GET(valorSobras))
             << setw(15) << FO << setw(15) << (MyClock() - inicio) / CLOCKS_PER_SEC
             << setw(12) << fabs(cplexGap) * 100 << "%"
             << setw(10) << numObjetos
             << setw(10) << numSobrasValidas
             << endl;
        cout << setw(0);
        if (periodo == P)
            cout << "Valor sobras remanescentes: " << CPLEX_GET(valorSobrasRemanescente) << endl;

        for(j = 1; j<=m; ++j) {
            if(CPLEX_GET(u[j]) == 1) {
                objeto[j].utilizado = true;
                ulong idx = objeto[j].codigo.objeto;
                instancia[periodo].objeto[idx].utilizado = true;
                for(i = 1; i<=n; ++i) {
                    if(CPLEX_GET(v[i][j]) == 1) {
                        ulong per = objeto[j].id.periodo;
                        ulong obj = objeto[j].id.objeto;
                        Item it = item[i];
                        it.x = CPLEX_GETD(x[i]) - item[i].w / 2.0;
                        it.y = CPLEX_GETD(y[i]) - item[i].h / 2.0;
                        it.x += objeto[j].x;
                        it.y += objeto[j].y;
//                        cout << setprecision(2) << fixed << i << ": (" << it.x << ", " << it.y  << ")\n";
                        listaItensDoObjeto[per][obj].push_back(it);
                    }
                }
            }
        }

        if(periodo < P) {
            for (j = 1; j <= m; ++j) {
                ulong per = objeto[j].id.periodo;
                ulong obj = objeto[j].id.objeto;
                /*********** Verifica se a sobra pode ser adicionada ************/
                if(periodo - objeto[j].id.periodo >= xi)
                    continue;
                if (CPLEX_GET(u[j]) == 1){
                    j1 = 2 * j - 1;
                    j2 = 2 * j;
                    ulong WsobraSuperior = CPLEX_GET(W_bar[1][j1]);
                    ulong HsobraSuperior = CPLEX_GET(H_bar[1][j1]);
                    ulong WsobraDireita  = CPLEX_GET(W_bar[1][j2]);
                    ulong HsobraDireita  = CPLEX_GET(H_bar[1][j2]);
                    ulong r_j = CPLEX_GET(r[j]);
                    ulong t_j = CPLEX_GET(t[j]);
                    double x_j, y_j;
                    if(WsobraSuperior*HsobraSuperior > 0){
                        x_j = objeto[j].x + 0;
                        y_j = objeto[j].y + (objeto[j].h - t_j);
                        Objeto novoObjeto(WsobraSuperior, HsobraSuperior, objeto[j].c, x_j, y_j, objeto[j].id, objeto[j].c_bar, true, objeto[j].codigo, objeto[j].codigoObjetoOrigem);
                        instancia[periodo + 1].adicionaObjetoSobraSuperior(novoObjeto);
                        listaItensDoObjeto[per][obj].push_back(Item(WsobraSuperior, HsobraSuperior, x_j, y_j, periodo, true));
                    }
                    if(WsobraDireita*HsobraDireita > 0){
                        x_j = objeto[j].x + (objeto[j].w - r_j);
                        y_j = objeto[j].y + 0;
                        Objeto novoObjeto(WsobraDireita, HsobraDireita, objeto[j].c, x_j, y_j, objeto[j].id, objeto[j].c_bar, true, objeto[j].codigo, objeto[j].codigoObjetoOrigem);
                        instancia[periodo + 1].adicionaObjetoSobraDireita(novoObjeto);
                        listaItensDoObjeto[per][obj].push_back(Item(WsobraDireita, HsobraDireita, x_j, y_j, periodo, true));
                    }
                }
                    // Se o objeto é originado de uma sobra aproveitável, então a sobra é passada para o próximo período
                    // Neste caso, o objeto é colocado na posição equivalente a sobra da direita
                else if (objeto[j].id.periodo != periodo) {
                    instancia[periodo + 1].adicionaObjetoSobraDireita(objeto[j]);
                    listaItensDoObjeto[per][obj].push_back(Item(objeto[j].w, objeto[j].h, objeto[j].x, objeto[j].y, periodo, true));
                }
            }
        }

        //Atualiza a área da sobra aproveitável de cada objeto (de *entrada*) do período
        for (j = 1; j <= m; ++j) {
            j1 = 2 * j - 1;
            j2 = 2 * j;
            ulong idx = objeto[j].codigo.objeto;
            if ( CPLEX_GET(u[j]) == 1){
                instancia[periodo].objeto[idx].areaSobraSuperior = CPLEX_GET(gamma[j1]);
                instancia[periodo].objeto[idx].areaSobraDireita = CPLEX_GET(gamma[j2]);
            }
            else if(objeto[j].id.periodo != periodo){
                instancia[periodo].objeto[idx].areaSobraDireita = objeto[j].w*objeto[j].h;
                instancia[periodo].objeto[idx].areaSobraSuperior = 0;
            }
        }

        //***Atualiza a área utilizada de objetos (originados de sobras?)***//
        //Contabilia a área dos itens cortados no objeto
        for(j = 1; j<=m; ++j) {
            if(CPLEX_GET(u[j]) == 1) {
                for(i = 1; i<=n; ++i) {
                    if( CPLEX_GET(v[i][j]) == 1 ){
                        ulong idx = objeto[j].codigo.objeto;
                        instancia[periodo].objeto[idx].areaUtilizada += (item[i].w*item[i].h);
                    }
                }
            }
        }

        if(periodo == P) {
            for (j = 1; j <= m; ++j) {
                ulong per = objeto[j].id.periodo;
                ulong obj = objeto[j].id.objeto;
                /*********** Verifica se a sobra pode ser adicionada ************/
                if(periodo - objeto[j].id.periodo >= xi)
                    continue;
                if ((CPLEX_GET(u[j]) == 1)) {
                    ulong r_j = CPLEX_GET(r[j]);
                    ulong t_j = CPLEX_GET(t[j]);
                    j1 = 2 * j - 1;
                    j2 = 2 * j;
                    ulong WsobraSuperior = CPLEX_GET(W_bar[1][j1]);
                    ulong HsobraSuperior = CPLEX_GET(H_bar[1][j1]);
                    ulong WsobraDireita = CPLEX_GET(W_bar[1][j2]);
                    ulong HsobraDireita = CPLEX_GET(H_bar[1][j2]);
                    double x_j, y_j;
                    ulong w, h;
                    //Contabilia a área das sobra aproveitáveis remanescentes dos objetos dos perídos anteriores à P
                    if (WsobraSuperior * HsobraSuperior > 0) {
                        x_j = objeto[j].x + 0;
                        y_j = objeto[j].y + (objeto[j].h - t_j);
                        listaItensDoObjeto[per][obj].push_back(
                                Item(WsobraSuperior, HsobraSuperior, x_j, y_j, periodo, true));
                    }
                    if (WsobraDireita * HsobraDireita > 0) {
                        x_j = objeto[j].x + (objeto[j].w - r_j);
                        y_j = objeto[j].y + 0;
                        listaItensDoObjeto[per][obj].push_back(
                                Item(WsobraDireita, HsobraDireita, x_j, y_j, periodo, true));
                    }
                    instancia[per].objeto[obj].areaSobraRemanescente +=
                            WsobraSuperior * HsobraSuperior + WsobraDireita * HsobraDireita;
                }
                else if((objeto[j].id.periodo != P)) { //Objeto originado de sobra mas não utilizado
                    listaItensDoObjeto[per][obj].push_back(Item(objeto[j].w, objeto[j].h, objeto[j].x, objeto[j].y, periodo, true));
                    instancia[per].objeto[obj].areaSobraRemanescente += objeto[j].w*objeto[j].h;
                }
            }
        }
#ifdef IMPRIMIR_INFO
        cout << endl;
#endif
        env.end();
        return finished;
    }
    catch (IloException& ex){
        cerr << "Error: " << ex << endl;
        env.end();
        exit(0);
    }
    catch (...){
        cerr << "Error" << endl;
        env.end();
        exit(0);
    }

    env.end();
    return false;
}

void print(int *v, int ini, int fim){
    for(int i = ini; i<=fim; ++i)
        cout << v[i] << "\t";
    cout << endl;
}

void print(vetor(int)& v, int ini, int fim){
    for(int i = ini; i<=fim; ++i)
        cout << v[i] << "\t";
    cout << endl;
}

/*
* P
* n1 n2 ... nP
* w[1][1] w[1][2] ... w[1][n1]
* h[1][1] h[1][2] ... h[1][n2]
* ...
* w[nP][1] w[nP][2] ... w[nP][n1]
* h[nP][1] h[nP][2] ... h[nP][n2]
*
* d
* w_bar[1] w_bar[2] ... w_bar[d]
* h_bar[1] w_bar[2] ... w_bar[d]
*
* m1 m2 ... mP
* W[1][1] W[1][2] ... W[1][m]
* H[1][1] H[1][2] ... H[1][m]
* c[1][1] c[1][2] ... c[1][m]
* ...
* W[mP][1] W[mP][2] ... W[mP][m]
* H[mP][1] H[mP][2] ... H[mP][m]
* c[mP][1] c[mP][2] ... c[mP][m]
*/

void leEntrada(string arquivoEntrada){
    int i, j, j1, s;
    int minw=INF, minh=INF;
    try{
        ifstream file(arquivoEntrada.c_str());
        if(!file){
            cerr << "Erro ao abrir o arquivo '" << arquivoEntrada << "'!" << endl;
            exit(0);
        }
        //Número de periodos
        file >> P;
        instancia.resize((unsigned int)P+2);
        //Número de itens por periodo
        //n.resize((unsigned long)(P+1));
        for(s = 1; s<=P; ++s)
            file >> instancia[s].n;

        //item.resize((unsigned long)(P+1));
        for(s = 1; s<=P; ++s)
            instancia[s].item.resize((unsigned int)(instancia[s].n+1));

        for(s = 1; s<=P; ++s){
            instancia[s].periodo = s;
            for(j = 1; j <=instancia[s].n;  ++j) {
                file >> instancia[s].item[j].w;
                minw = (instancia[s].item[j].w<minw?instancia[s].item[j].w:minw);
            }
            for(j = 1; j <=instancia[s].n; ++j) {
                file >> instancia[s].item[j].h;
                minh = (instancia[s].item[j].h<minh?instancia[s].item[j].h:minh);
            }
            for(j = 1; j <=instancia[s].n; ++j)
                instancia[s].item[j].id = s;
        }


        for(s = 1; s<=P; ++s)
            sort(instancia[s].item.begin()+1, instancia[s].item.end(),
                 [](Item i1, Item i2){return (i1.getArea() > i2.getArea()?true:(i1.getArea()==i2.getArea()?(i1.w>i2.w):false));});

        for(s = 1; s<=P; ++s){
            instancia[s].p = 1;
            for(j = 1; j <instancia[s].n; ++j){
                if(instancia[s].item[j] == instancia[s].item[j+1])
                    continue;
                instancia[s].p++;
            }
//            cout << p[s] << endl;
            instancia[s].n_q.resize((unsigned long)(instancia[s].p+1));
            instancia[s].o_q.resize((unsigned long)(instancia[s].p+1));
            int o = 1;
            int k = 1;
            for(j = 1; j <instancia[s].n; ++j){
                if(instancia[s].item[j] == instancia[s].item[j+1])
                    o++;
                else{
                    instancia[s].n_q[k++] = o;
                    o = 1;
                }
            }
//            cout << instancia[s].p << endl;
            instancia[s].n_q[k] = o;
//            print(instancia[s].n_q, 1, instancia[s].p); cout << endl;
            instancia[s].o_q[1] = 0;
//            if(2 <= p[s]) o_q[s][2] = n_q[s][1];
            for(j = 2; j <=instancia[s].p; ++j){
//                o_q[s][i] = o_q[s][i-1] + n_q[s][i];
                instancia[s].o_q[j] = instancia[s].o_q[j -1] + instancia[s].n_q[j -1];
            }
//            print(instancia[s].o_q, 1, instancia[s].p);
//            cout << endl;
        }

        //Itens do catálogo usados apenas no último período
        //Nos outros períodos, usa-se os itens do próximo período
        file >> instancia[P].d;
        instancia[P].w_bar.resize((unsigned long)(instancia[P].d+1));
        instancia[P].h_bar.resize((unsigned long)(instancia[P].d+1));
        for(j = 1; j <=instancia[P].d; ++j)
            file >> instancia[P].w_bar[j];
        for(j = 1; j <=instancia[P].d; ++j)
            file >> instancia[P].h_bar[j];

        int opcao = OPCAO_ITENS_CATALOGO;
        if(opcao == 1){
            for(s = 1; s<P; ++s){
                instancia[s].d = instancia[P].d;
                instancia[s].w_bar.resize((unsigned long)(instancia[s].d+1));
                instancia[s].h_bar.resize((unsigned long)(instancia[s].d+1));
                for(i = 1; i<=instancia[P].d; ++i){
                    instancia[s].w_bar[i] = (instancia[P].w_bar[i]);
                    instancia[s].h_bar[i] = (instancia[P].h_bar[i]);
                }
            }
        }
        else if (opcao == 2){
            for(s = 1; s<P; ++s){
                instancia[s].d = instancia[s+1].n+instancia[P].d;
                instancia[s].w_bar.resize((unsigned long)(instancia[s].d+1));
                instancia[s].h_bar.resize((unsigned long)(instancia[s].d+1));
                for(i = 1; i<=instancia[s+1].n; ++i){
                    instancia[s].w_bar[i] = instancia[s+1].item[i].w;
                    instancia[s].h_bar[i] = instancia[s+1].item[i].h;
                }
                for(i = instancia[s+1].n+1; i<=instancia[P].d+instancia[s+1].n; ++i){
                    instancia[s].w_bar[i] = (instancia[P].w_bar[i-instancia[s+1].n]);
                    instancia[s].h_bar[i] = (instancia[P].h_bar[i-instancia[s+1].n]);
                }
            }
        }
        else if(opcao == 3){
            for(s = 1; s<P; ++s){
                cout << "Periodo: " << s << endl;
                instancia[s].w_bar.push_back(0);
                instancia[s].h_bar.push_back(0);
                for(int s1 = s+1; s1<=(MIN(s+2,P)); ++s1){
                    for(i = 1; i<=instancia[s1].n; ++i) {
                        instancia[s].w_bar.push_back(instancia[s1].item[i].w);
                        instancia[s].h_bar.push_back(instancia[s1].item[i].h);
                        cout << instancia[s1].item[i].w << " x " << instancia[s1].item[i].h << endl;
                    }
                }
                for(i = 1; i<=instancia[P].d; ++i){
                    bool adicionado = false;
                    int tam = (int)instancia[s].w_bar.size();
                    for(int k=1; k<tam; ++k){
                        if(instancia[s].w_bar[k] == instancia[P].w_bar[i] && instancia[s].h_bar[k] == instancia[P].h_bar[i]) {
                            adicionado = true;
                            break;
                        }
                    }
                    if(!adicionado) {
                        instancia[s].w_bar.push_back(instancia[P].w_bar[i]);
                        instancia[s].h_bar.push_back(instancia[P].h_bar[i]);
                    }
                }
                instancia[s].d = (int)instancia[s].w_bar.size()-1;
            }
        }
        else {
            //Para cada período, adiciona como itens do catálogo os itens demandados futuros mais o itens do catálogo de entrada
            for(s = 1; s<P; ++s){
                instancia[s].w_bar.push_back(0);
                instancia[s].h_bar.push_back(0);
                for(int s1 = s+1; s1<=P/*(MIN(s+2,P))*/; ++s1){
                    for(i = 1; i<=instancia[s1].n; ++i) {
                        int tam = (int)instancia[s].w_bar.size();
                        bool adicionado = false;
                        for(int k=1; k<tam; ++k){
                            if(instancia[s].w_bar[k] == instancia[s1].item[i].w && instancia[s].h_bar[k] == instancia[s1].item[i].h) {
                                adicionado = true;
                                break;
                            }
                        }
                        if(!adicionado) {
                            instancia[s].w_bar.push_back(instancia[s1].item[i].w);
                            instancia[s].h_bar.push_back(instancia[s1].item[i].h);
                        }
                    }
                }
                //Adiciona os itens do catálogo de entrada
                for(i = 1; i<=instancia[P].d; ++i){
                    bool adicionado = false;
                    int tam = (int)instancia[s].w_bar.size();
                    for(int k=1; k<tam; ++k){
                        if(instancia[s].w_bar[k] == instancia[P].w_bar[i] && instancia[s].h_bar[k] == instancia[P].h_bar[i]) {
                            adicionado = true;
                            break;
                        }
                    }
                    if(!adicionado) {
                        instancia[s].w_bar.push_back(instancia[P].w_bar[i]);
                        instancia[s].h_bar.push_back(instancia[P].h_bar[i]);
                    }
                }
                instancia[s].d = (int)instancia[s].w_bar.size()-1;
            }
        }

        for(s = 1; s<=P; ++s) {
            file >> instancia[s].m;
            instancia[s].mHR = instancia[s].m;
        }
        instancia[P + 1].m = 0;
        instancia[1].objeto.resize((unsigned long) (instancia[1].m + 1));
        instancia[1].m_bar = instancia[1].m;
        int soma;
        for(s = 2; s <= P + 1; ++s) {
            soma = 0;
            for(int l = 1; l <= MIN(s,xi); l++){
                soma += (int)pow(2, l)*instancia[s - l].m;
            }
            instancia[s].m_bar = instancia[s].m + soma;
//            instancia[s].m_bar = instancia[s].m + 2*instancia[s-1].m_bar;
            instancia[s].objeto.resize((unsigned long) (instancia[s].m_bar + 1));
//            cerr << s << ": " << instancia[s].m << " " << instancia[s].m_bar << " " << instancia[s].objeto.size()-1 << endl;
        }

        for(s = 1; s <= P + 1; ++s){
//            cerr << "Periodo: " << s << endl;
            if (s <= P) {
                for (j = 1; j <= instancia[s].m; ++j) {
                    file >> instancia[s].objeto[j].w;
                    WMax = (instancia[s].objeto[j].w > WMax ? instancia[s].objeto[j].w : WMax);
                }
                for (j = 1; j <= instancia[s].m; ++j) {
                    file >> instancia[s].objeto[j].h;
                    HMax = (instancia[s].objeto[j].h > HMax ? instancia[s].objeto[j].h : HMax);
                }
                for (j = 1; j <= instancia[s].m; ++j) {
                    file >> instancia[s].objeto[j].c;
                    instancia[s].objeto[j].c_bar = instancia[s].objeto[j].c;
                }
            }
            for(j = 1; j <= instancia[s].m; ++j) {
                instancia[s].objeto[j].id = Par{s, j};
                instancia[s].objeto[j].codigo = instancia[s].objeto[j].codigoObjetoOrigem = Par{s, j};
                instancia[s].objeto[j].areaSobraDireita = instancia[s].objeto[j].w*instancia[s].objeto[j].h;
                instancia[s].objeto[j].areaSobraSuperior = instancia[s].objeto[j].w*instancia[s].objeto[j].h;
                instancia[s].objeto[j].delta = Delta(valorInicialDelta, valorInicialDelta, valorInicialDelta, valorInicialDelta);
            }
            if(s > 1) {
                for (j = instancia[s].m + 1, j1 = 1; j <= instancia[s].m_bar; j += 2, j1++) {
                    instancia[s].objeto[j].id = instancia[s - 1].objeto[j1].id;
                    instancia[s].objeto[j].c = instancia[s - 1].objeto[j1].c;
                    instancia[s].objeto[j].c_bar = instancia[s - 1].objeto[j1].c;
                    instancia[s].objeto[j].codigo = Par{s, j};
                    instancia[s].objeto[j].codigoObjetoOrigem = Par{s - 1, j1};
                    instancia[s].objeto[j].ehSobra = true;
                    instancia[s].objeto[j].delta = Delta(valorInicialDelta, valorInicialDelta, valorInicialDelta,
                                                         valorInicialDelta);

                    instancia[s].objeto[j + 1].id = instancia[s - 1].objeto[j1].id;
                    instancia[s].objeto[j + 1].c = instancia[s - 1].objeto[j1].c;
                    instancia[s].objeto[j + 1].c_bar = instancia[s - 1].objeto[j1].c;
                    instancia[s].objeto[j + 1].codigo = Par{s, j + 1};
                    instancia[s].objeto[j + 1].codigoObjetoOrigem = Par{s - 1, j1};
                    instancia[s].objeto[j + 1].ehSobra = true;
                    instancia[s].objeto[j + 1].delta = Delta(valorInicialDelta, valorInicialDelta, valorInicialDelta,
                                                             valorInicialDelta);

                }
            }
        }

        instancia[P+1].m = 0;
        instancia[P+1].m_bar = instancia[P+1].m + 2 * instancia[P].m_bar; 
        instancia[P+1].objeto.resize((unsigned long) (instancia[P+1].m_bar + 1));

        refObjPer.resize((unsigned long) (P + 2));
        for (s = 1; s <= P + 1; s++) {
            refObjPer[s].resize((unsigned long) (instancia[s].m_bar + 1));
            for (j = 1; j <= instancia[s].m; ++j) {
                refObjPer[s][j] = Par{s, j};
            }
        }
        for (s = 2; s <= P + 1; s++) {
            for (j = instancia[s].m + 1, j1 = 1; j <= instancia[s].m_bar; j += 2, j1++) {
                refObjPer.at(s).at(j) = refObjPer.at(s - 1).at(j1);
            }
        }

    }
    catch(exception &e){
        cout << "Error: " << e.what() << endl;
        exit(0);
    }
}

bool atualizaAreaUtilizadaObjeto(Instancia& instance){

    int s, j, per, obj;

    for(s=P; s>=1; s--){
        for(j=instance[s].m+1; j<=instance[s].m_bar; j+=2){
            //if(!instance[s].objeto[j].ehObjetoValido()) continue;

            per = instance[s].objeto[j].codigoObjetoOrigem.periodo;
            obj = instance[s].objeto[j].codigoObjetoOrigem.objeto;

            //if(!instance[per].objeto[obj].utilizado) continue;

            instance[per].objeto[obj].areaSobraDireitaUtilizada = instance[s].objeto[j].areaUtilizada;

            instance[per].objeto[obj].areaSobraSuperiorUtilizada = instance[s].objeto[j+1].areaUtilizada;

            instance[per].objeto[obj].atualizaAreaUtilizada();
        }
    }

    for(s=1; s<=P; s++) {
        for (j = 1; j <= instance[s].m_bar; j++) {

            if(!instance[s].objeto[j].ehObjetoValido()) continue;

            if(!instance[s].objeto[j].utilizado) continue;

            instance[s].objeto[j].atualizaDeltas();
        }
    }

    vector< tuple<int, int, Delta> > deltas;
    for(s=1; s <= P; s++) {
        int numObj = instance[s].m_bar;
        for (j = 1; j <= numObj; j++) {
            deltas.emplace_back(make_tuple(s, j, instance[s].objeto[j].delta));
        }
    }
    valoresDelta.emplace_back(deltas);

    return !deltaRepetidos();
}

bool cabeAlgumItemNoObjeto(Objeto& obj, vetor( vetor(Item) )& item){
    ulong numPeriodos = (ulong) item.size()-1;
    for(ulong s = 1; s<=numPeriodos; ++s) {
        ulong numItens = (ulong) item[s].size() - 1;
        for (ulong i = 1; i <= numItens; ++i) {
            if (item[s][i].w <= obj.w && item[s][i].h <= obj.h)
                return true;
        }
    }
    return false;
}

bool cabeAlgumItemNoObjeto(Objeto& obj, vetor(Item)& item){
    ulong numItens = (ulong) item.size() - 1;
    for (ulong i = 1; i <= numItens; ++i) {
        if (item[i].w <= obj.w && item[i].h <= obj.h)
            return true;
    }
    return false;
}

void setParametros(IloCplex &cplex, bool solucaoInteira){
    if(setLimiteMemoriaCplex) {
        cplex.setParam(IloCplex::WorkDir, ".");
        cplex.setParam(IloCplex::WorkMem, maxMemoriaCplex);
        cplex.setParam(IloCplex::NodeFileInd, 3);
    }
    if(setLimiteThreadsCplex)
        cplex.setParam(IloCplex::Threads, numThreadsCplex);

    cplex.setParam(IloCplex::TiLim, tempoLimiteCplex);
    cplex.setParam(IloCplex::ClockType, 2); //1 CPU time. 2 Wall clock time (total physical time elapsed).
    //cplex.setParam(IloCplex::EpOpt, 1e-9);
    if(solucaoInteira)
        cplex.setParam(IloCplex::EpAGap, 0.999999); //default: 1e-06
    else
        cplex.setParam(IloCplex::EpAGap, 1e-7);  //default: 1e-06
    cplex.setParam(IloCplex::EpGap, 0.0); //default: 1e-04
    cplex.setParam(IloCplex::PreInd, preSolveCplex);
    //cplex.setParam(IloCplex::NumericalEmphasis, true);
    //cplex.setParam(IloCplex::AdvInd, 2);
    //cplex.setParam(IloCplex::Symmetry, 0);

    /* Value	Meaning
     *   0	    Balance optimality and feasibility; default
     *   1	    Emphasize feasibility over optimality
     *   2	    Emphasize optimality over feasibility
     *   3	    Emphasize moving best bound
     *   4	    Emphasize finding hidden feasible solutions
     */
    //cplex.setParam(IloCplex::MIPEmphasis, 1);

//    cplex.setParam(IloCplex::RepeatPresolve, 3);
//    cplex.setParam(IloCplex::RelaxPreInd, 1);
//    cplex.setParam(IloCplex::HeurFreq, 20);
}


inline bool is_file_exist(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

void adicionaPontoDeCorte(vetor(Ponto)& v, Ponto& pt){
    int tam = (int) v.size();
    for(int i=0; i<tam; ++i)
        if(v[i] == pt)
            return;
    v.push_back(pt);
}

bool cabeItemCatalogo(vetor(int) &w_bar, vetor(int) &h_bar, Objeto& sobra){
    if(sobra.getArea() <= 0) return false;
    unsigned long tam = w_bar.size();
    for(unsigned long i=1; i<tam; ++i){
        if(w_bar[i] <= sobra.w && h_bar[i] <= sobra.h)
            return true;
    }
    return false;
}

void imprimeDeltas(){
    int k;
    ullong numLinhas = valoresDelta.size();
    ullong numColunas = valoresDelta[0].size();
    cout << setprecision(5) << fixed;
    vector<bool> colunaUsada(numColunas, false);
    cout << endl;
    cout << "DELTAS SOBRA DIREITA" << endl;
    k = 1;
    //vector<bool> colunaUsada(numColunas, false);
    for (ullong c = 0; c < numColunas; ++c) {
        for (ullong l = 0; l < numLinhas; ++l) {
            if (!DOUBLE_EQUALS(get<2>(valoresDelta[l][c]).sobraDireita, valorInicialDelta)) {
                colunaUsada[c] = true;
            }
        }
    }
    cout << "it\t";
    for (auto v: valoresDelta) {
        int c = 0;
        for (auto d: v) {
            if (colunaUsada[c])
                cout << "O^" << get<0>(d) << "_" << get<1>(d) << "\t";
            c++;
        }
        break;
    }
    cout << endl;
    for (auto v: valoresDelta) {
        cout << k << "\t";
        int c = 0;
        for (auto d: v) {
            if (colunaUsada[c])
                cout << get<2>(d).sobraDireita << "\t";
            c++;
        }
        cout << endl;
        k++;
    }
    cout << endl;

    cout << "DELTAS SOBRA SUPERIOR" << endl;
    k = 1;
    colunaUsada.assign(numColunas, false);
    for (ullong c = 0; c < numColunas; ++c) {
        for (ullong l = 0; l < numLinhas; ++l) {
            if (!DOUBLE_EQUALS(get<2>(valoresDelta[l][c]).sobraSuperior, valorInicialDelta)) {
                colunaUsada[c] = true;
            }
        }
    }
    cout << "it\t";
    for (auto v: valoresDelta) {
        int c = 0;
        for (auto d: v) {
            if (colunaUsada[c])
                cout << "O^" << get<0>(d) << "_" << get<1>(d) << "\t";
            c++;
        }
        break;
    }
    cout << endl;
    for (auto v: valoresDelta) {
        cout << k << "\t";
        int c = 0;
        for (auto d: v) {
            if (colunaUsada[c])
                cout << get<2>(d).sobraSuperior << "\t";
            c++;
        }
        cout << endl;
        k++;
    }
    cout << endl;
}

bool deltaRepetidos(){
    ullong numLinhas = valoresDelta.size();
    ullong numColunas = valoresDelta[0].size();
    bool deltaRepetidosSobras = false;
    bool deltaRepetidosObjeto = false;
    bool deltaRepetidosSobraDireita = false;
    bool deltaRepetidosSobraSuperior = false;

    //Verifica se o último conjunto de deltas já foi adicionado
#ifdef IMPRIMIR_INFO
    cout << "SOBRA DIREITA: LINHA " << numLinhas;
#endif
    for (ullong l2 = 0; l2 < numLinhas - 1; ++l2) {
        bool deltarepetidos = true;
        for (ullong c = 0; c < numColunas; ++c) {
            if (!COMPARA_DELTAS(get<2>(valoresDelta[numLinhas - 1][c]).sobraDireita,
                                get<2>(valoresDelta[l2][c]).sobraDireita)) {
                deltarepetidos = false;
                break;
            }
        }
        if (deltarepetidos){
#ifdef IMPRIMIR_INFO
            cout << " IGUAL A LINHA " << l2+1;
#endif
            deltaRepetidosSobraDireita = true;
            break;
        }
    }
#ifdef IMPRIMIR_INFO
    cout << endl;
    cout << "SOBRA SUPERIOR: LINHA " << numLinhas;
#endif
    for (ullong l2 = 0; l2 < numLinhas - 1; ++l2) {
        bool deltarepetidos = true;
        for (ullong c = 0; c < numColunas; ++c) {
            if (!COMPARA_DELTAS(get<2>(valoresDelta[numLinhas - 1][c]).sobraSuperior,
                                get<2>(valoresDelta[l2][c]).sobraSuperior)) {
                deltarepetidos = false;
                break;
            }
        }
        if (deltarepetidos){
#ifdef IMPRIMIR_INFO
            cout << " IGUAL A LINHA " << l2+1;
#endif
            deltaRepetidosSobraSuperior = true;
            break;
        }
    }
#ifdef IMPRIMIR_INFO
    cout << endl;
#endif

    return (deltaRepetidosSobraDireita && deltaRepetidosSobraSuperior);
}
