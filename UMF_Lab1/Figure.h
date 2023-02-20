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

class Figure
{
public:
   virtual bool isInnerNode(type x, type y, type eps) = 0;
   virtual bool isEdgeNode(type x, type y, type eps) = 0;
   virtual bool IsLeft(type x, type y, type eps) = 0;
   virtual bool IsRight(type x, type y, type eps) = 0;
   virtual bool IsBottom(type x, type y, type eps) = 0;
   virtual bool IsTop(type x, type y, type eps) = 0;
};

class LshapedFigure : public Figure {
private:
   type scale;
   type ratio;
public:
   bool isInnerNode(type x, type y, type eps) override;
   bool isEdgeNode(type x, type y, type eps) override;
   bool IsLeft(type x, type y, type eps) override;
   bool IsRight(type x, type y, type eps) override;
   bool IsBottom(type x, type y, type eps) override;
   bool IsTop(type x, type y, type eps) override;
};

class Square : public Figure {
private:
   type scale;
public:
   Square(type scale);
   bool isEdgeNode(type x, type y, type eps) override;
   bool isInnerNode(type x, type y, type eps) override;
   bool IsLeft(type x, type y, type eps) override;
   bool IsRight(type x, type y, type eps) override;
   bool IsBottom(type x, type y, type eps) override;
   bool IsTop(type x, type y, type eps) override;
};


