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
typedef vector<vector<type>> Matrix;


class Area
{
public:
	type x1, x2, y1, y2;
	int lambda, gamma, fId; // - функция 
};

