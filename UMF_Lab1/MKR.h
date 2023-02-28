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
   const int count_of_diag = 5;
   Vector hx, hy; // коэффиценты приращения
   Vector kx, ky; // коэффициент разрядик
   int countX,countY, countArea;
   Vector vectorX, vectorY; //векторы x и y;
   vector<Area> vectorAreas;
   vector<Node> grid;
   vector<Boundary> vectorB1, vectorB2x, vectorB2y;
   Matrix Aij;
   Vector f;
   type func(Node point);
   MKR();
   int sub(Vector& a, Vector& b);
   type norm(Vector& a);
   void Boundaries();
   void Areas();
   void Grid();
   int makeMatrix();
   int input();
   type step(Vector l1, Vector l2, Vector u1, Vector u2, Vector& di, Vector& b, Vector& x0, Vector& x, type w,
      int& N, int& m, Vector& z, int i);
   // Блочная релаксация с блоком 1x1
   type GaussZeidel(Vector l1, Vector l2, Vector u1, Vector u2, Vector& di, Vector& f, Vector& x0, type& eps, type& w, int& N, int& m, int& max_iter);
   type teta(Node point, int idBound);
   int solve(string filename);

   /*
   int BlockRelaxation();
   int BlockSize();
   void Factorize();
   void BlockRelazationIteration(int blocksize);
   type ResidualBlock(int i);
   Vector Ri_(int blocksize, int i);
   type FindAElement(int a1, int a2);
   Vector LU_SOL(int i, int blocksize, Vector& Ri);
   */
};

