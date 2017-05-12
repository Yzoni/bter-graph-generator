//
// Created by Yorick de Boer on 11/05/2017.
//

#ifndef BTER_SEQUENTIAL_BTEREXPORT_H
#define BTER_SEQUENTIAL_BTEREXPORT_H


class BterExport {
private:
    void createMap(std::map<int, vector<int>> *m);

public:
    void exportAdjencyList();
};


#endif //BTER_SEQUENTIAL_BTEREXPORT_H
