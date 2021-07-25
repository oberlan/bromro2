#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <iomanip>

using namespace std;

typedef unsigned long ulong;

#define randm(x)  (rand()%(ulong)(x))
#define rint(a,b) (randm((b)-(a)+1) + (a))
#define MAX(X,Y)  ((X)>(Y) ? (X) : (Y))

class Item{
public:
    ulong w, h;
    double x, y;
    ulong id;
    Item(ulong _w=0, ulong _h=0, double _x=0, double _y=0, ulong _id=0):
            w(_w),h(_h),x(_x),y(_y),id(_id)
    {
    }
    bool operator<(const Item& i) const {
        if(w < i.w) return true;
        else if(w == i.w) return (h < i.h);
        else return false;
    }

    void print() const {
        cout << w << " x " << h << "\t(" << x << ", " << y << ")\t";
    }
    bool operator==(const Item it) const{
        return (w==it.w && h==it.h);
    }
    ulong area(){
        return w*h;
    }
};

//Parâmetros de entrada
ulong semente = 1180;
ulong numMinPeriodo=4;
ulong numMaxPeriodo=8;
ulong numMinTipoObjetos=1;
ulong numMaxTipoObjetos=3;
ulong limInferiorDimensaoObjetos=10;
ulong limSuperiorDimensaoObjetos=50;
ulong numMinTipoItens=1;
ulong numMaxTipoItens=5;
ulong numMinItensPorTipo=1;
ulong numMaxItensPorTipo=3;
ulong limInferiorDimensaoItens=5;
ulong limSuperiorDimensaoItens=20;

//Valores sorteados
ulong numPeriodos;

vector<ulong> numTipoObjetos;
vector<ulong> numTotalObjetosPorPeriodo;
vector< vector<ulong> > numObjetosPorTipo;
vector< vector<ulong> > W; //Comprimento de cada objeto do período
vector< vector<ulong> > H; //Altura de cada objeto do período
vector< vector<int> > c;   //Custo por unidade de área do objeto

vector<ulong> numTipoItens;
vector<ulong> numTotalItensPorPeriodo;
vector< vector<Item> > item; //Vetor de itens do período
vector<Item> catalogo;


//Variáveis auxiliares
ulong WMax, HMax;

bool itemCabeEmUmObjeto(vector<ulong> W, vector<ulong>H, ulong w, ulong h);
void imprimeInstancia();
ulong empacotaItens(ulong W, ulong H, vector<Item> item);
vector<Item> itensDominates();
void leArgumentos(int argc, char** argv);
void geraCodigoLatex();
void imprimeParametrosUsados();

int main(int argc, char** argv) {

    leArgumentos(argc, argv);

    srand((unsigned int)semente);
    //Sortea o número de períodos
    numPeriodos = rint(numMinPeriodo, numMaxPeriodo);

    //Redimensiona o tamanho dos vetores
    W.resize(numPeriodos+1);
    H.resize(numPeriodos+1);
    c.resize(numPeriodos+1);
    numObjetosPorTipo.resize(numPeriodos+1);
    numTipoObjetos.resize(numPeriodos+1, 0);
    item.resize(numPeriodos+1);
    numTipoItens.resize(numPeriodos+1, 0);
    numTotalItensPorPeriodo.resize(numPeriodos + 1, 0);
    numTotalObjetosPorPeriodo.resize(numPeriodos + 1, 0);

    for(ulong p=1; p<=numPeriodos; ++p){

        //Define o número de tipo de objetos e suas respectivas dimensões
        numTipoObjetos[p] = rint(numMinTipoObjetos, numMaxTipoObjetos);
        W[p].resize(numTipoObjetos[p]+1, 0);
        H[p].resize(numTipoObjetos[p]+1, 0);
        c[p].resize(numTipoObjetos[p]+1, 0);
        numObjetosPorTipo[p].resize(numTipoObjetos[p]+1, 0);
        WMax = 0;
        HMax = 0;
        for(ulong i=1; i<=numTipoObjetos[p]; ++i){
            W[p][i] = rint(limInferiorDimensaoObjetos, limSuperiorDimensaoObjetos);
            H[p][i] = rint(limInferiorDimensaoObjetos, limSuperiorDimensaoObjetos);
            c[p][i] = 1;//((p<=numPeriodos/2)?1:2);
            WMax = MAX(WMax, W[p][i]);
            HMax = MAX(HMax, H[p][i]);
        }
        //Define o número de tipo de itens e suas respectivas dimensões
        numTipoItens[p] = rint(numMinTipoItens, numMaxTipoItens);
        item[p].push_back(Item());
        for(ulong i=1; i<=numTipoItens[p]; ++i){
            ulong numItemPorTipo = rint(numMinItensPorTipo, numMaxItensPorTipo);
            ulong wtmp, htmp;
            //Difine as dimensões do item e verifica se o mesmo cabe em pelo menos um objeto
            do {
                wtmp = rint(limInferiorDimensaoItens, limSuperiorDimensaoItens);
                htmp = rint(limInferiorDimensaoItens, limSuperiorDimensaoItens);
            }while(!itemCabeEmUmObjeto(W[p], H[p], wtmp, htmp));
            for(ulong j=1; j<=numItemPorTipo; ++j){
                item[p].push_back(Item(wtmp,htmp));
                numTotalItensPorPeriodo[p]++;
            }
        }

        //Empacota os itens em cada objeto para definir a quantidade de objetos de cada tipo
        for(ulong i=1; i<=numTipoObjetos[p]; ++i){
            numObjetosPorTipo[p][i] = empacotaItens(W[p][i], H[p][i], item[p]);
            numTotalObjetosPorPeriodo[p] += numObjetosPorTipo[p][i];
        }
    }

    //Define os itens dominates
    catalogo = itensDominates();

    imprimeInstancia();

    //imprimeParametrosUsados();

    //geraCodigoLatex();

    return 0;
}

void imprimeInstancia(){
    cout << numPeriodos << endl << endl;
    for(ulong p=1; p<=numPeriodos; ++p){
        cout << numTotalItensPorPeriodo[p] << "\t";
    }
    cout << endl;
    for(ulong p=1; p<=numPeriodos; ++p){
        for(ulong i=1; i<item[p].size(); ++i)
            cout << item[p][i].w << "\t";
        cout << endl;
        for(ulong i=1; i<item[p].size(); ++i)
            cout << item[p][i].h << "\t";
        cout << endl;
    }
    cout << endl;
    cout << endl;
    cout << catalogo.size() << endl;
    for(auto i: catalogo)
        cout << i.w << "\t";
    cout << endl;
    for(auto i: catalogo)
        cout << i.h << "\t";
    cout << endl;
    cout << endl;
    for(ulong p=1; p<=numPeriodos; ++p){
        cout << numTotalObjetosPorPeriodo[p] << "\t";
    }
    cout << endl;
    for(ulong p=1; p<=numPeriodos; ++p){
        for(ulong i=1; i<W[p].size(); ++i)
            for(ulong j=1; j<=numObjetosPorTipo[p][i]; ++j)
                cout << W[p][i] << "\t";
        cout << endl;
        for(ulong i=1; i<H[p].size(); ++i)
            for(ulong j=1; j<=numObjetosPorTipo[p][i]; ++j)
                cout << H[p][i] << "\t";
        cout << endl;
        for(ulong i=1; i<H[p].size(); ++i)
            for(ulong j=1; j<=numObjetosPorTipo[p][i]; ++j)
                cout << c[p][i] << "\t";
        cout << endl;
    }
}

bool itemCabeEmUmObjeto(vector<ulong> W, vector<ulong>H, ulong w, ulong h){
    ulong n = W.size();
    //if(w == 0 || h == 0) return false;
    for(ulong i=1; i<n; ++i){
        if(w <= W[i] && h <= H[i]) return true;
    }
    return false;
}

ulong empacotaItens(ulong W, ulong H, vector<Item> item){
    if(item.size() == 0){
        cout << "Erro ao empacotar itens" << endl;
        exit(0);
    }
    ulong numItens = item.size()-1;
    ulong numObjetos = 0;
    //Ordena os itens de forma decrescente pela área
    sort(item.begin()+1, item.end(), [](Item i1, Item i2){return i1.area() > i2.area();});
//    for(size_t i=1; i<=numItens; ++i){
//        cout << item[i].w << "\t" << item[i].h << "\t" << item[i].area() << endl;
//    }
    vector<bool> itemEmpacotado(numItens+1, false);
    itemEmpacotado[0] = true;
    while(true){
        numObjetos++;
        vector<ulong> objeto(W, 0);

        for(ulong i=1; i<=numItens; ++i) {
            //Verifica se o item cabe no objeto
            if(item[i].w > W || item[i].h > H) {
                itemEmpacotado[i] = true;
                continue;
            }
            if(itemEmpacotado[i]) continue;
            for (ulong j = 0; j < W; ++j) {
                bool achouEspaco = false;
                if(objeto[j] + item[i].h <= H && j+item[i].w-1 < W){
                    for (ulong k = j+1; k < item[i].w + j; ++k){
                        if(objeto[k] + item[i].h > H) {
                            break;
                        }
                    }
                    achouEspaco = true;
                }
                if(achouEspaco){
                    item[i].x = j;
                    item[i].y = objeto[j];
                    //item[i].print(); cout << endl;
                    for (ulong k = j; k < item[i].w + j; ++k){
                        objeto[k] += item[i].h;
                    }
                    itemEmpacotado[i] = true;
                    break;
                }
            }
        }

        //Verifica se empacotou todos os itens
        if (find (itemEmpacotado.begin()+1, itemEmpacotado.end(), false) == itemEmpacotado.end()) {
            break;
        }
    }

    return numObjetos;

}

vector<Item> itensDominates(){
    vector<Item> itens, todosItens;
    vector<bool> itemDominante;
    for(ulong p=1; p<=numPeriodos; ++p){
        for(ulong i=1; i<item[p].size(); ++i){
            bool adicionado = false;
            for(auto it: todosItens){
                if(it == item[p][i]){
                    adicionado = true;
                    break;
                }
            }
            if(!adicionado) {
                todosItens.push_back(item[p][i]);
                itemDominante.push_back(true);
            }
        }
    }

    for(ulong i=0; i<todosItens.size(); ++i){
        for(ulong j=0; j<todosItens.size(); ++j){
            if(i == j) continue;
            if(todosItens[i].w <= todosItens[j].w && todosItens[i].h <= todosItens[j].h){
                itemDominante[j] = false;
            }
        }
    }

    for(ulong i=0; i<todosItens.size(); ++i) {
        if (itemDominante[i]){
            itens.push_back(todosItens[i]);
        }
    }
    return itens;
}

void leArgumentos(int argc, char** argv) {
    //Valores padrões para os parâmetros
    semente = 10;                    //-s
    numMinPeriodo=4;                //-p1
    numMaxPeriodo=8;                //-p2
    numMinTipoObjetos=1;            //-p3
    numMaxTipoObjetos=3;            //-p4
    limInferiorDimensaoObjetos=10;  //-p5
    limSuperiorDimensaoObjetos=50;  //-p6
    numMinTipoItens=1;              //-p7
    numMaxTipoItens=5;              //-p8
    numMinItensPorTipo=1;           //-p9
    numMaxItensPorTipo=3;           //-p10
    limInferiorDimensaoObjetos=5;   //-p11
    limSuperiorDimensaoObjetos=20;  //-p12

    for(int i = 1 ; i < argc ; i++) {
        if (string(argv[i]) == "-s") {
            semente = (ulong)atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-p1") {
            numMinPeriodo = (ulong)atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-p2") {
            numMaxPeriodo = (ulong)atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-p3") {
            numMinTipoObjetos = (ulong)atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-p4") {
            numMaxTipoObjetos = (ulong)atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-p5") {
            limInferiorDimensaoObjetos = (ulong)atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-p6") {
            limSuperiorDimensaoObjetos = (ulong)atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-p7") {
            numMinTipoItens = (ulong)atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-p8") {
            numMaxTipoItens = (ulong)atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-p9") {
            numMinItensPorTipo = (ulong)atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-p10") {
            numMaxItensPorTipo = (ulong)atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-p11") {
            limInferiorDimensaoItens = (ulong)atoi(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-p12") {
            limSuperiorDimensaoItens = (ulong)atoi(argv[i + 1]);
            continue;
        }
    }
}

void geraCodigoLatex(){
    cout << endl << endl << endl;
    cout << "---- CODIGO LATEX DA INSTANCIA ----" << endl << endl;
    vector< vector < pair<Item, int> > > itensDistintos(numPeriodos+1);
    for(ulong p=1; p<=numPeriodos; ++p){
        vector<bool> itemAdicionado(item[p].size()+1, false);
        for(ulong i=1; i<item[p].size(); ++i){
            if(!itemAdicionado[i]) {
                itensDistintos[p].push_back(make_pair(item[p][i], 1));
                itemAdicionado[i] = true;
                for(ulong j=i+1; j<item[p].size(); ++j){
                    if(item[p][i] == item[p][j]){
                        itemAdicionado[j] = true;
                        itensDistintos[p][itensDistintos[p].size()-1].second++;
                    }
                }
            }
        }
    }
    cout << "\\multirow{" << numPeriodos << "}{*}{XXX} & " << "\\multirow{" << numPeriodos << "}{*}{" << numPeriodos << "} ";
    for(ulong p=1; p<=numPeriodos; ++p){
        if(p > 1) cout << " & ";
        cout << "& " << numTotalObjetosPorPeriodo[p] << " & ";
        for(ulong i=1; i<W[p].size(); ++i){
            if(numObjetosPorTipo[p][i] > 1) cout << numObjetosPorTipo[p][i];
            if(numObjetosPorTipo[p][i] > 1 || c[p][i] > 1) cout << "(";
            cout << W[p][i] << " $\\times$ " << H[p][i];
            if(numObjetosPorTipo[p][i] > 1 || c[p][i] > 1) cout << ")";
            if(c[p][i] > 1) cout << "[" << c[p][i] << "]";
            if(i<W[p].size()-1) cout << ", ";
        }

        cout << " & " << numTotalItensPorPeriodo[p] << " & " << numTipoItens[p] << " & ";
        if(p==1) cout << "\\multirow{" << numPeriodos << "}{*}{" << catalogo.size() << "} & ";
        else cout << " & ";
        for(ulong i=0 ; i<itensDistintos[p].size(); ++i){
            if(itensDistintos[p][i].second > 1) cout << itensDistintos[p][i].second << "(";
            bool itemDoCatalogo = false;
            for(auto it: catalogo){
                if(it == itensDistintos[p][i].first){
                    itemDoCatalogo = true;
                    break;
                }
            }
            if(itemDoCatalogo) cout << "\\underline{";
            cout << itensDistintos[p][i].first.w << " $\\times$ " << itensDistintos[p][i].first.h;
            if(itemDoCatalogo) cout << "}";
            if(itensDistintos[p][i].second > 1) cout << ")";
            if(i<itensDistintos[p].size()-1) cout << ", ";
        }
        if(p < numPeriodos) cout << "\\\\" << endl;
        else cout << "\\\\";
    }
    cout << "\\hline" << endl;



    cerr << "\\multirow{" << numPeriodos << "}{*}{XXX} & " << "\\multirow{" << numPeriodos << "}{*}{" << numPeriodos << "} ";
    for(ulong p=1; p<=numPeriodos; ++p){
        if(p > 1) cerr << " & ";
        cerr << "& " << numTotalObjetosPorPeriodo[p] << " & ";
        for(ulong i=1; i<W[p].size(); ++i){
            if(numObjetosPorTipo[p][i] > 1) cerr << numObjetosPorTipo[p][i];
            if(numObjetosPorTipo[p][i] > 1 || c[p][i] > 1) cerr << "(";
            cerr << W[p][i] << " $\\times$ " << H[p][i];
            if(numObjetosPorTipo[p][i] > 1 || c[p][i] > 1) cerr << ")";
            if(c[p][i] > 1) cerr << "[" << c[p][i] << "]";
            if(i<W[p].size()-1) cerr << ", ";
        }

        cerr << " & " << numTotalItensPorPeriodo[p] << " & " << numTipoItens[p] << " & ";
        if(p==1) cerr << "\\multirow{" << numPeriodos << "}{*}{" << catalogo.size() << "} & ";
        else cerr << " & ";
        for(ulong i=0 ; i<itensDistintos[p].size(); ++i){
            if(itensDistintos[p][i].second > 1) cerr << itensDistintos[p][i].second << "(";
            bool itemDoCatalogo = false;
            for(auto it: catalogo){
                if(it == itensDistintos[p][i].first){
                    itemDoCatalogo = true;
                    break;
                }
            }
            if(itemDoCatalogo) cerr << "\\underline{";
            cerr << itensDistintos[p][i].first.w << " $\\times$ " << itensDistintos[p][i].first.h;
            if(itemDoCatalogo) cerr << "}";
            if(itensDistintos[p][i].second > 1) cerr << ")";
            if(i<itensDistintos[p].size()-1) cerr << ", ";
        }
        if(p < numPeriodos) cerr << "\\\\" << endl;
        else cerr << "\\\\";
    }
    cerr << " \\hline" << endl << endl;
}

void imprimeParametrosUsados(){
    cout << "\n\n\n";
    cout << "+------------------------------------------------+" << endl;
    cout << "|    PARÂMETROS USADOS PARA GERAR A INSTÂNCIA    |" << endl;
    int x = 3;
    cout << "|------------------------------------------+-----|" << endl;
    cout << "| semente:                                 | " << setw(x) << semente                    << " | " << endl;
    cout << "|------------------------------------------+-----|" << endl;
    cout << "| numMinPeriodo:                           | " << setw(x) << numMinPeriodo              << " | " << endl;
    cout << "|------------------------------------------+-----|" << endl;
    cout << "| numMaxPeriodo:                           | " << setw(x) << numMaxPeriodo              << " | " << endl;
    cout << "|------------------------------------------+-----|" << endl;
    cout << "| numMinTipoObjetos:                       | " << setw(x) << numMinTipoObjetos          << " | " << endl;
    cout << "|------------------------------------------+-----|" << endl;
    cout << "| numMaxTipoObjetos:                       | " << setw(x) << numMaxTipoObjetos          << " | " << endl;
    cout << "|------------------------------------------+-----|" << endl;
    cout << "| limInferiorDimensaoObjetos:              | " << setw(x) << limInferiorDimensaoObjetos << " | " << endl;
    cout << "|------------------------------------------+-----|" << endl;
    cout << "| limSuperiorDimensaoObjetos:              | " << setw(x) << limSuperiorDimensaoObjetos << " | " << endl;
    cout << "|------------------------------------------+-----|" << endl;
    cout << "| numMinTipoItens:                         | " << setw(x) << numMinTipoItens            << " | " << endl;
    cout << "|------------------------------------------+-----|" << endl;
    cout << "| numMaxTipoItens:                         | " << setw(x) << numMaxTipoItens            << " | " << endl;
    cout << "|------------------------------------------+-----|" << endl;
    cout << "| numMinItensPorTipo:                      | " << setw(x) << numMinItensPorTipo         << " | " << endl;
    cout << "|------------------------------------------+-----|" << endl;
    cout << "| numMaxItensPorTipo:                      | " << setw(x) << numMaxItensPorTipo         << " | " << endl;
    cout << "+------------------------------------------+-----+" << endl;
    cout << endl;
}
