#ifndef CUSTOM_VECTOR_H
#define CUSTOM_VECTOR_H

#include <stdexcept>
#include <vector>
using namespace std;

template <class T>
class CustomVector {
   public:
    int size;
    vector<int> sizes;
    vector<vector<T>> vectors;

    CustomVector() : size(0) {}

    void extend(vector<T>& v) {
        int s = v.size();
        sizes.push_back(s);
        size += s;
        vectors.push_back(v);
    }

    T get(int index) const {
        if (index >= size || index < 0) {
            throw out_of_range("Index out of range");
        }

        int offset = 0;
        int currentTotal = 0;

        for (int i = 0; i < sizes.size(); i++) {
            currentTotal += sizes[i];
            if (currentTotal > index) {
                offset = index - (currentTotal - sizes[i]);
                return vectors[i][offset];
            }
        }

        throw out_of_range("Index out of range");
    }

    typename std::vector<T>::iterator getPointer(int index) {
        if (index >= size) {
            throw out_of_range("Index out of range");
        }

        int offset = 0;
        int currentTotal = 0;

        for (int i = 0; i < sizes.size(); i++) {
            currentTotal += sizes[i];
            if (currentTotal > index) {
                offset = index - (currentTotal - sizes[i]);
                return (vectors[i].begin() + offset);
            }
        }

        throw out_of_range("Index out of range");
    }
};

#endif  // CUSTOM_VECTOR_H