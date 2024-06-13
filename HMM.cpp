#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>

using namespace std;

bool isMatch(const vector<char>& col) {
    int noGaps = count_if(col.begin(), col.end(), [](char c) { return c != '-'; });
    return noGaps > col.size() / 2;
}

unordered_map<int, unordered_map<char, double>> calcProbEmision(const vector<string>& secuencias, const vector<int>& colsCoincidencia) {
    string simbolos = "ACDEFGHIJKLMNOPQRSTUVWXYZ-";
    int numEstados = colsCoincidencia.size();
    unordered_map<int, unordered_map<char, int>> conteoEmision;

    for (int i = 0; i < numEstados; ++i) {
        for (char s : simbolos) {
            conteoEmision[i][s] = 0;
        }
    }

    for (const string& sec : secuencias) {
        for (int i = 0; i < numEstados; ++i) {
            char simbolo = sec[colsCoincidencia[i]];
            conteoEmision[i][simbolo]++;
        }
    }

    unordered_map<int, unordered_map<char, double>> probEmision;
    for (int i = 0; i < numEstados; ++i) {
        int total = 0;
        for (char s : simbolos) {
            total += conteoEmision[i][s] + 1;
        }
        for (char s : simbolos) {
            probEmision[i][s] = (conteoEmision[i][s] + 1.0) / total;
        }
    }

    return probEmision;
}

unordered_map<int, unordered_map<char, double>> calcProbTransicion(const vector<string>& secuencias, const vector<int>& colsCoincidencia) {
    int numEstados = colsCoincidencia.size();
    unordered_map<int, unordered_map<char, int>> conteoTransicion;

    for (int i = 0; i <= numEstados; ++i) {
        conteoTransicion[i]['M'] = 0;
        conteoTransicion[i]['X'] = 0;
        conteoTransicion[i]['Y'] = 0;
    }

    for (const string& sec : secuencias) {
        for (int i = 0; i < numEstados; ++i) {
            char simbolo = sec[colsCoincidencia[i]];
            if (simbolo == '-') {
                conteoTransicion[i]['Y']++;
            } else {
                conteoTransicion[i]['M']++;
            }
        }
    }

    unordered_map<int, unordered_map<char, double>> probTransicion;
    for (int i = 0; i <= numEstados; ++i) {
        int total = conteoTransicion[i]['M'] + conteoTransicion[i]['X'] + conteoTransicion[i]['Y'] + 3;
        probTransicion[i]['M'] = (conteoTransicion[i]['M'] + 1.0) / total;
        probTransicion[i]['X'] = (conteoTransicion[i]['X'] + 1.0) / total;
        probTransicion[i]['Y'] = (conteoTransicion[i]['Y'] + 1.0) / total;
    }

    return probTransicion;
}

vector<int> getMatch(const vector<string>& secuencias) {
    int numCols = secuencias[0].size();
    vector<int> colsCoincidencia;
    for (int col = 0; col < numCols; ++col) {
        vector<char> columna;
        for (const string& sec : secuencias) {
            columna.push_back(sec[col]);
        }
        if (isMatch(columna)) {
            colsCoincidencia.push_back(col);
        }
    }
    return colsCoincidencia;
}

void HMM(const vector<string>& secuencias) {
    vector<int> colsCoincidencia = getMatch(secuencias);
    unordered_map<int, unordered_map<char, double>> probEmision = calcProbEmision(secuencias, colsCoincidencia);
    unordered_map<int, unordered_map<char, double>> probTransicion = calcProbTransicion(secuencias, colsCoincidencia);

    cout << "Columnas de coincidencia: ";
    for (size_t i = 0; i < colsCoincidencia.size(); ++i) {
        cout << colsCoincidencia[i] << " ";
    }
    cout << endl;

    cout << "Probabilidades de emisión:" << endl;
    for (unordered_map<int, unordered_map<char, double>>::iterator it = probEmision.begin(); it != probEmision.end(); ++it) {
        cout << "Estado " << it->first << ": ";
        for (unordered_map<char, double>::iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {
            cout << jt->first << " -> " << jt->second << ", ";
        }
        cout << endl;
    }

    cout << "Probabilidades de transición:" << endl;
    for (unordered_map<int, unordered_map<char, double>>::iterator it = probTransicion.begin(); it != probTransicion.end(); ++it) {
        cout << "Estado " << it->first << ": ";
        for (unordered_map<char, double>::iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {
            cout << jt->first << " -> " << jt->second << ", ";
        }
        cout << endl;
    }
}

int main() {
    vector<string> secuencias = {
        "VGA--HAGEY",
        "V----NVDEV",
        "VEA--DVAGH",
        "VKG------D",
        "VYS--TYETS",
        "FNA--NIPKH",
        "IAGADNGAGY"
    };

    HMM(secuencias);

    return 0;
}
