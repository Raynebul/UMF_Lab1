#pragma once
#include "Matrix.h"
#include "Area.h"
#include "Boundary.h"
#include "Node.h"
#include <cmath>

//typedef vector<vector<type>> Matrix;


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
   Matrix Aij;
   MKR();
   void Boundaries();
   void Areas();
   void Grid();
   int makeMatrix();
   int input();
   void GaussZeidel();
   void BlockRelaxation();

};

