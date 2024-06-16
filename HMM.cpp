#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>
#include <string>

using namespace std;

struct EmissionProbabilities {
    map<int, unordered_map<char, double>> matchState; 
};

struct TransitionProbabilities {
    unordered_map<string, double> transitions;
};

void initializeProbabilities(TransitionProbabilities &transition, int numColumns) {
    for (int i = 1; i < numColumns; ++i) {
        transition.transitions["M" + to_string(i) + "toM" + to_string(i + 1)] = 0.0;
        transition.transitions["M" + to_string(i) + "toI" + to_string(i)] = 0.0;
        transition.transitions["M" + to_string(i) + "toD" + to_string(i)] = 0.0;
        transition.transitions["I" + to_string(i) + "toM" + to_string(i + 1)] = 0.0;
        transition.transitions["I" + to_string(i) + "toI" + to_string(i)] = 0.0;
        transition.transitions["D" + to_string(i) + "toM" + to_string(i + 1)] = 0.0;
        transition.transitions["D" + to_string(i) + "toD" + to_string(i + 1)] = 0.0;
    }
}

vector<int> selectColumns(const vector<string> &sequences) {
    vector<int> selectedColumns;
    int numSequences = sequences.size();
    int sequenceLength = sequences[0].length();
    int threshold = numSequences / 2;

    for (int i = 0; i < sequenceLength; ++i) {
        int symbolCount = 0;
        for (int j = 0; j < numSequences; ++j) {
            if (sequences[j][i] != '-') {
                symbolCount++;
            }
        }
        if (symbolCount >= threshold) {
            selectedColumns.push_back(i);
        }
    }

    return selectedColumns;
}

void calculateEmissionProbabilities(const vector<string> &sequences, const vector<int> &selectedColumns, EmissionProbabilities &emission) {
    int numSequences = sequences.size();
    int alphabetSize = 20; //  20 aminoácidos
    int pseudoCount = 1;

    for (int col : selectedColumns) {
        unordered_map<char, int> counts;
        int totalCount = 0;

        for (int j = 0; j < numSequences; ++j) {
            char symbol = sequences[j][col];
            if (symbol != '-') {
                counts[symbol]++;
                totalCount++;
            }
        }

        for (char c = 'A'; c <= 'Z'; ++c) {
            int count = counts[c];
            double probability = (count + pseudoCount) / double(totalCount + pseudoCount * alphabetSize);
            emission.matchState[col][c] = probability;
        }
    }
}

void calculateTransitionProbabilities(const vector<string> &sequences, const vector<int> &selectedColumns, TransitionProbabilities &transition) {
    int numSequences = sequences.size();
    int pseudoCount = 1;
    int stateCount = 3; // M, I, D

    unordered_map<string, int> transitionCounts;
    unordered_map<string, int> stateCounts;

    for (int i = 1; i < selectedColumns.size(); ++i) {
        transitionCounts["M" + to_string(i) + "toM" + to_string(i + 1)] = 0;
        transitionCounts["M" + to_string(i) + "toI" + to_string(i)] = 0;
        transitionCounts["M" + to_string(i) + "toD" + to_string(i)] = 0;
        stateCounts["M" + to_string(i)] = 0;
    }

    for (int seqIndex = 0; seqIndex < numSequences; ++seqIndex) {
        for (int i = 1; i < selectedColumns.size(); ++i) {
            int colIndex = selectedColumns[i];
            int prevColIndex = selectedColumns[i - 1];

            if (sequences[seqIndex][prevColIndex] != '-') {
                stateCounts["M" + to_string(i)]++;
                if (sequences[seqIndex][colIndex] != '-') {
                    transitionCounts["M" + to_string(i) + "toM" + to_string(i + 1)]++;
                } else {
                    transitionCounts["M" + to_string(i) + "toD" + to_string(i)]++;
                }
            } else {
                transitionCounts["M" + to_string(i) + "toI" + to_string(i)]++;
            }
        }
    }

    for (int i = 1; i < selectedColumns.size(); ++i) {
        string mToM = "M" + to_string(i) + "toM" + to_string(i + 1);
        string mToI = "M" + to_string(i) + "toI" + to_string(i);
        string mToD = "M" + to_string(i) + "toD" + to_string(i);

        transition.transitions[mToM] = (transitionCounts[mToM] + pseudoCount) / double(stateCounts["M" + to_string(i)] + pseudoCount * stateCount);
        transition.transitions[mToI] = (transitionCounts[mToI] + pseudoCount) / double(stateCounts["M" + to_string(i)] + pseudoCount * stateCount);
        transition.transitions[mToD] = (transitionCounts[mToD] + pseudoCount) / double(stateCounts["M" + to_string(i)] + pseudoCount * stateCount);
    }
}


void HMM(const vector<string> &sequences) {
    EmissionProbabilities emission;
    TransitionProbabilities transition;

    vector<int> selectedColumns = selectColumns(sequences);
    initializeProbabilities(transition, selectedColumns.size());
    calculateEmissionProbabilities(sequences, selectedColumns, emission);
    calculateTransitionProbabilities(sequences, selectedColumns, transition);


    cout << "Probabilidades de emisión (estado de coincidencia):" << endl;
    for (auto &colEntry : emission.matchState) {
        cout << "Columna " << colEntry.first << ":" << endl;
        for (auto &entry : colEntry.second) {
            cout << "bm" << colEntry.first + 1 << entry.first << " = " << entry.second << endl;
        }
    }

    cout << "Probabilidades de transición:" << endl;
    for (int i = 1; i < selectedColumns.size(); ++i) {
        string mToM = "M" + to_string(i) + "toM" + to_string(i + 1);
        string mToI = "M" + to_string(i) + "toI" + to_string(i);
        string mToD = "M" + to_string(i) + "toD" + to_string(i);

        cout << mToM << " = " << transition.transitions[mToM] << endl;
        cout << mToI << " = " << transition.transitions[mToI] << endl;
        cout << mToD << " = " << transition.transitions[mToD] << endl;
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

}
