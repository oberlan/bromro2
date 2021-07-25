//
// Created by oberlan on 23/02/16.
//

#ifndef PROJETO_INSTANCIAPORPERIODO_H
#define PROJETO_INSTANCIAPORPERIODO_H

#include <vector>

#include "Util.h"
#include "Item.h"
#include "Objeto.h"

using namespace std;


struct InstanciaPorPeriodo {
    ulong n; //Número de itens
    ulong m; //Número de objetos de entrada
    ulong mHR;
    ulong m_bar; //Número de objetos de entrada mais número de sobras do período anterior (m[p]+2*m[p-1]+1)
    ulong d; //Quantidade de itens do catálogo
    ulong p; //Produtos distintos do período
    ulong periodo; //Periodo dos dados
    vetor(ulong) w_bar; //Comprimento dos itens do catálogo (itens do catálogo + itens do proximo período)
    vetor(ulong) h_bar; //Altura dos itens do catálogo (itens do catálogo + itens do proximo período)
    vetor(ulong) n_q;
    vetor(ulong) o_q;
    vetor(Item) item; //Itens do período
    vetor(Objeto) objeto; //ListaObjetos do perído
    vetor(Objeto) getObjetosValidos() {
        unsigned long objValidos = 0;
        for (ulong j = 1; j <= m_bar; j++)
            if (objeto[j].ehObjetoValido()) objValidos++;
        vetor(Objeto) obj(objValidos);
        for (ulong j = 1, j1 = 1; j <= m_bar; j++)
            if (objeto[j].ehObjetoValido()) {
                obj[j1] = objeto[j];
                j1++;
            }
        return obj;
    }

    unsigned int getNumeroObjetosValidos() {
        unsigned int objValidos = 0;
        for (ulong j = 1; j <= m_bar; j++)
            if (objeto[j].ehObjetoValido()) objValidos++;
        return objValidos;
    }

    void resetaInstancia() {
        //ASSERT(m_bar == objeto.size());
        for (ulong j = 1; j <= m; j++) {
            objeto[j].utilizado = false;
            objeto[j].deltaMudou = false;
            objeto[j].areaSobraDireita = 0;
            objeto[j].areaSobraSuperior = 0;
            objeto[j].areaUtilizada = 0.0;
            objeto[j].areaSobraDireitaUtilizada = 0.0;
            objeto[j].areaSobraSuperiorUtilizada = 0.0;
            objeto[j].areaSobraRemanescente = 0;
        }
        for (ulong j = m + 1; j <= m_bar; j++) {
            objeto[j].reset();
        }
    }

    ulong adicionaObjetoSobraDireita(Objeto &sobra) {
        if (sobra.codigo.periodo + 1 != periodo) {
            cerr << sobra.codigo.periodo << " " << periodo << endl;
        }
        ASSERT(sobra.codigo.periodo + 1 == periodo);
        ulong idx = ((sobra.codigo.objeto) * 2 - 1) + m;
        if (!(idx <= m_bar)) {
            cerr << " SD: " << idx << " " << m_bar << "\n";
            cerr << periodo << " " << sobra.id << " " << sobra.codigo << endl;
        }
        ASSERT(idx <= m_bar);
        ASSERT(mHR <= m_bar);
        objeto[idx].w = sobra.w;
        objeto[idx].h = sobra.h;
        objeto[idx].c = sobra.c;
        objeto[idx].c_bar = sobra.c_bar;
        objeto[idx].x = sobra.x;
        objeto[idx].y = sobra.y;
        objeto[idx].id = sobra.id;
        return idx;
    }

    ulong adicionaObjetoSobraSuperior(Objeto &sobra) {
        ASSERT(sobra.codigo.periodo + 1 == periodo);
        ulong idx = ((sobra.codigo.objeto) * 2 - 1) + m + 1;
        if (!(idx <= m_bar)) {
            cerr << "SS: " << idx << " " << m_bar << endl;
        }
        ASSERT(idx <= m_bar);
        ASSERT(mHR <= m_bar);
        objeto[idx].w = sobra.w;
        objeto[idx].h = sobra.h;
        objeto[idx].c = sobra.c;
        objeto[idx].c_bar = sobra.c_bar;
        objeto[idx].x = sobra.x;
        objeto[idx].y = sobra.y;
        objeto[idx].id = sobra.id;
        return idx;
    }

    int adicionaObjetoHR(Objeto &sobra) {
        mHR++;
        ulong idx = mHR;
        ASSERT(idx <= m_bar);
        objeto[idx].w = sobra.w;
        objeto[idx].h = sobra.h;
        objeto[idx].c = sobra.c;
        objeto[idx].c_bar = sobra.c_bar;
        objeto[idx].x = sobra.x;
        objeto[idx].y = sobra.y;
        objeto[idx].ehSobra = sobra.ehSobra;
        objeto[idx].id = sobra.id;
        return idx;
    }
};

#endif //PROJETO_INSTANCIAPORPERIODO_H
