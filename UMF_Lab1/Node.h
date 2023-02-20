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

class Node
{
public:
   type x, y; // Координаты точек
 //  type kx, ky; // Коэффициент разрядки
  // type hx, hy; // Коэффициент приращения
  // int Nx, Ny; // Кол-во разбиений
};

