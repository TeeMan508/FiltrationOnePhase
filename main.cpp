#include <iostream>
#include "utils/params/params.h"

int main() {
    Parameters params;
    params.load("../data/example.txt");

    return 0;
}
