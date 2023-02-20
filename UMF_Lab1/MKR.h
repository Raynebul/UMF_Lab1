#pragma once
#include "Grid.h"
#include "Area.h"
#include "Boundary.h"
#include <cmath>

using namespace std;

typedef double type;
typedef vector<type> Vector;
typedef vector<vector<type>> Matrix;


class MKR
{
public:
   Vector hx, hy; // ����������� ����������
   Vector kx, ky; // ����������� ��������
   int countX,countY, countArea;
   Vector vectorX, vectorY; //������� x � y;
   vector<Area> vectorAreas;
   vector<Node> grid;
   vector<Boundary> vectorB1, vectorB2x, vectorB2y;
   MKR();
   void Boundaries();
   void Areas();
   void Grid();

   int input();
   void GaussZeidel();
   void BlockRelaxation();

};

