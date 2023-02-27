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

class Area
{
public:
	type x1, x2, y1, y2; // ??? спорно очень
	int lambda, gamma, fId; // - функция 
};

