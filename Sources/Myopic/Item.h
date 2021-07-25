//
// Created by Oberlan Romão on 05/09/15.
//

#ifndef ITEM_H
#define ITEM_H
#include <iostream>

using namespace std;

class Item{
public:
    int w, h, id;
    double x, y;
    bool ehSobra;
    Item(int _w=0, int _h=0, double _x=0, double _y=0, int _id=0, bool sobra=false):
            w(_w),h(_h),x(_x),y(_y),id(_id),ehSobra(sobra)
    {
    }
    bool operator<(const Item& i) const {
        if(w < i.w) return true;
        else if(w == i.w) return (h < i.h);
        else return false;
    }

    void print() const {
        cout << w << " x " << h << "\t(" << x << ", " << y << ")\n";
    }
    bool operator==(const Item it) const{
        return (w==it.w && h==it.h);
    }
    inline int getArea() const{
        return w*h;
    }
};
#endif //ITEM_H
