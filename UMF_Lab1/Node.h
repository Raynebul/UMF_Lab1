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

class Node
{
public:
   type x, y; // Координаты точек
   int area;
};

