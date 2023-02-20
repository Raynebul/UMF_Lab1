#pragma once
#include <iostream>
#include <windows.h>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

using namespace std;

typedef double type;
typedef vector<type> Vector;


class Matrix
{
public:
   vector<double> di, u1, u2, l1, l2;
   int k = -1;
};

