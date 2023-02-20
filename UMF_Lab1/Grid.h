#pragma once
#include "Figure.h"
#include "Node.h"

class setka
{
public:
   enum class tests {Test1, Test2, Test3, Test4, Test5, Test6};
   type x, y;
   Vector hx, hy; // коэффиценты приращения
   Vector kx, ky; // коэффициенты разрядки
   Matrix u; //матрица
   //// выбор функции
   type function(tests test);
   type diff_function(tests test);
};

enum class ThirdConditionsSide {
   LEFT,
   RIGHT,
   BOTTOM,
   TOP,
   NONE,
};
// Равномерная сетка
class RegularGrid {
private:
   double _scale;
   double _stepX;
   double _stepY;
   Figure* _figure;
   function2D _u;
   function2D _f;
   ThirdConditionsSide _side;
   double _epsBase = 1e-5;
   double _beta = 1.0;
   double _lambda = 1.0;
   double _gamma = 0.0;
public:
   RegularGrid(double scale, double stepX, double stepY, Figure* figure,
      function2D f, function2D u, ThirdConditionsSide side);
   template<typename T>
   T convertToSystemOfEquations(SystemOfEquationsFactory<T>* factory);
private:
   vector<Node> nodes();
   double leftDerivativeByX(double x, double y);
   double leftDerivativeByY(double x, double y);
   double rightDerivativeByX(double x, double y);
   double rightDerivativeByY(double x, double y);
};
RegularGrid::RegularGrid(double scale, double stepX, double stepY, Figure*
   figure, function2D f, function2D u, ThirdConditionsSide side) {
   _scale = scale;
   _stepX = stepX;
   _stepY = stepY;
   _figure = figure;
   _u = u;
   _f = f;
   _side = side;
}
vector<Node> RegularGrid::nodes() {
   vector<Node> nodes = vector<Node>();
   int numberOfNodesByX = int(_scale / _stepX);
   int numberOfNodesByY = int(_scale / _stepY);
   double eps = _stepX < _stepY ? _stepX * _epsBase : _stepY * _epsBase;
   for (int i = 0; i <= numberOfNodesByX; i++) {
      for (int j = 0; j <= numberOfNodesByY; j++) {
         double x = (double)i * _stepX;
         double y = (double)j * _stepY;
         if (_figure->isEdgeNode(x, y, eps)) {
            nodes.push_back(Node(NodeType::EDGE, x, y));
         }
         else {
            if (_figure->isInnerNode(x, y, eps)) {
               nodes.push_back(Node(NodeType::INNER, x, y));
            }
            else {
               nodes.push_back(Node(NodeType::DUMMY, x, y));
            }
         }
      }
   }
   return nodes;
}
double RegularGrid::leftDerivativeByX(double x, double y) {
   double h = 1e-9;
   return (_u(x, y) - _u(x - h, y)) / h;
}
double RegularGrid::leftDerivativeByY(double x, double y) {
   double h = 1e-9;
   return (_u(x, y) - _u(x, y - h)) / h;
}
double RegularGrid::rightDerivativeByX(double x, double y) {
   double h = 1e-9;
   return (_u(x + h, y) - _u(x, y)) / h;
}
double RegularGrid::rightDerivativeByY(double x, double y) {
   double h = 1e-9;
   return (_u(x, y + h) - _u(x, y)) / h;
}
template<typename T>
T RegularGrid::convertToSystemOfEquations(SystemOfEquationsFactory<T>* factory) {
   int width = int(_scale / _stepX) + 1;
   double eps = _stepX < _stepY ? _stepX * _epsBase : _stepY * _epsBase;
   vector<Node> nodes = this->nodes();
   int countOfNodes = nodes.size();
   vector<double> di = vector<double>(countOfNodes);
   vector<double> al1 = vector<double>(countOfNodes - 1);
   vector<double> al2 = vector<double>(countOfNodes - width);
   vector<double> au1 = vector<double>(countOfNodes - 1);
   vector<double> au2 = vector<double>(countOfNodes - width);
   vector<double> b = vector<double>(countOfNodes);
   for (int i = 0; i < countOfNodes; i++) {
      double x = nodes[i].x;
      double y = nodes[i].y;
      switch (nodes[i].type) {
      case NodeType::EDGE: {
         switch (_side) {
         case ThirdConditionsSide::NONE: {
            di[i] = 1.0;
            b[i] = _u(x, y);
            break;
         }
         case ThirdConditionsSide::LEFT: {
            if (_figure->isLeft(x, y, eps)) {
               auto uInPoint = _u(x, y);
               auto rightByX = rightDerivativeByX(x, y);
               auto ubeta = -_lambda * rightByX / _beta +
                  uInPoint;
               di[i] = _lambda / _stepY + _beta;
               au2[i] = -_lambda / _stepY;
               b[i] = -_lambda * rightByX + _beta * (uInPoint -
                  ubeta) + _beta * ubeta;
            }
            else {
               di[i] = 1.0;
               b[i] = _u(x, y);
            }
            break;
         }
         case ThirdConditionsSide::RIGHT: {
            if (_figure->isRight(x, y, eps)) {
               auto uInPoint = _u(x, y);
               auto leftByX = leftDerivativeByX(x, y);
               auto ubeta = _lambda * leftByX / _beta +
                  uInPoint;
               di[i] = _lambda / _stepY + _beta;
               al2[i - width] = -_lambda / _stepY;
               b[i] = _lambda * leftByX + _beta * (uInPoint -
                  ubeta) + _beta * ubeta;
            }
            else {
               di[i] = 1.0;
               b[i] = _u(x, y);
            }
            break;
         }
         case ThirdConditionsSide::BOTTOM: {
            if (_figure->isBottom(x, y, eps)) {
               auto uInPoint = _u(x, y);
               auto rightByY = rightDerivativeByY(x, y);
               auto ubeta = -_lambda * rightByY / _beta +
                  uInPoint;
               di[i] = _lambda / _stepY + _beta;
               au1[i] = -_lambda / _stepY;
               b[i] = -_lambda * rightByY + _beta * (uInPoint -
                  ubeta) + _beta * ubeta;
            }
            else {
               di[i] = 1.0;
               b[i] = _u(x, y);
            }
            break;
         }
         case ThirdConditionsSide::TOP: {
            if (_figure->isTop(x, y, eps)) {
               auto uInPoint = _u(x, y);
               auto leftByY = leftDerivativeByY(x, y);
               auto ubeta = _lambda * leftByY / _beta +
                  uInPoint;
               di[i] = _lambda / _stepY + _beta;
               al1[i - 1] = -_lambda / _stepY;
               b[i] = _lambda * leftByY + _beta * (uInPoint -
                  ubeta) + _beta * ubeta;
            }
            else {
               di[i] = 1.0;
               b[i] = _u(x, y);
            }
            break;
         }
         }
         break;
      }
      case NodeType::INNER: {
         di[i] = _lambda * (2.0 / (_stepX * _stepX) + 2.0 / (_stepY *
            _stepY)) + _gamma;
         al1[i - 1] = -_lambda / (_stepY * _stepY);
         au1[i] = -_lambda / (_stepY * _stepY);
         al2[i - width] = -_lambda / (_stepX * _stepX);
         au2[i] = -_lambda / (_stepX * _stepX);
         b[i] = _f(x, y);
         break;
      }
      case NodeType::DUMMY: {
         di[i] = 1.0;
         b[i] = 0.0;
         break;
      }
      }
   }
   return factory->createSystem(di, au1, au2, al1, al2, b);
}

//Неравномерная сетка
class IrregularGrid {
private:
   double _scale;
   double _stepX;
   double _stepY;
   Figure* _figure;
   function2D _u;
   function2D _f;
   ThirdConditionsSide _side;
   double _epsBase = 1e-5;
   double _beta = 1.0;
   double _lambda = 1.0;
   double _gamma = 0.0;
public:
   IrregularGrid(double scale, double stepX, double stepY, Figure* figure,
      function2D f, function2D u, ThirdConditionsSide side);
   template<typename T>
   T convertToSystemOfEquations(SystemOfEquationsFactory<T>* factory);
private:
   vector<Node> nodes();
   int width();
   //Левая разность
   double leftDerivativeByX(double x, double y);
   double leftDerivativeByY(double x, double y);
   //Правая разность
   double rightDerivativeByX(double x, double y);
   double rightDerivativeByY(double x, double y);
};
IrregularGrid::IrregularGrid(double scale, double stepX, double stepY, Figure*
   figure, function2D f, function2D u, ThirdConditionsSide side) {
   _scale = scale;
   _stepX = stepX;
   _stepY = stepY;
   _figure = figure;
   _f = f;
   _u = u;
   _side = side;
}
template<typename T>
T IrregularGrid::convertToSystemOfEquations(SystemOfEquationsFactory<T>* factory)
{
   vector<Node> nodes = this->nodes();
   int countOfNodes = nodes.size();
   int w = width();
   double eps = _stepX < _stepY ? _stepX * _epsBase : _stepY * _epsBase;
   vector<double> di = vector<double>(countOfNodes);
   vector<double> al1 = vector<double>(countOfNodes - 1);
   vector<double> al2 = vector<double>(countOfNodes - w);
   vector<double> au1 = vector<double>(countOfNodes - 1);
   vector<double> au2 = vector<double>(countOfNodes - w);
   vector<double> b = vector<double>(countOfNodes);
   for (int i = 0; i < countOfNodes; i++) {
      double x = nodes[i].x;
      double y = nodes[i].y;
      switch (nodes[i].type) {
      case NodeType::EDGE: {
         switch (_side) {
         case ThirdConditionsSide::NONE: {
            di[i] = 1.0;
            b[i] = _u(x, y);
            break;
         }
         case ThirdConditionsSide::BOTTOM: {
            if (_figure->isBottom(x, y, eps)) {
               auto uInPoint = _u(x, y);
               auto rightByY = rightDerivativeByY(x, y);
               auto ubeta = -_lambda * rightByY / _beta +
                  uInPoint;
               di[i] = _lambda / _stepY + _beta;
               au1[i] = -_lambda / _stepY;
               b[i] = -_lambda * rightByY + _beta * (uInPoint -
                  ubeta) + _beta * ubeta;
            }
            else {
               di[i] = 1.0;
               b[i] = _u(x, y);
            }
            break;
         }
         case ThirdConditionsSide::TOP: {
            if (_figure->isTop(x, y, eps)) {
               auto uInPoint = _u(x, y);
               auto leftByY = leftDerivativeByY(x, y);
               auto ubeta = _lambda * leftByY / _beta +
                  uInPoint;
               di[i] = _lambda / _stepY + _beta;
               al1[i - 1] = -_lambda / _stepY;
               b[i] = _lambda * leftByY + _beta * (uInPoint -
                  ubeta) + _beta * ubeta;
            }
            else {
               di[i] = 1.0;
               b[i] = _u(x, y);
            }
            break;
         }
         }
         break;
      }
      case NodeType::INNER: {
         double hx = nodes[i + w].x - nodes[i].x;
         double hy = nodes[i + 1].y - nodes[i].y;
         double hxPrev = nodes[i].x - nodes[i - w].x;
         double hyPrev = nodes[i].y - nodes[i - 1].y;
         //cout << y << "," << endl;
         di[i] = _lambda * (2.0 / (hxPrev * hx) + 2.0 / (hyPrev * hy))
            + _gamma;
         al1[i - 1] = -2.0 * _lambda / (hyPrev * (hy +
            hyPrev));//(hxPrev * (hx + hxPrev));
         au1[i] = -2.0 * _lambda / (hy * (hy + hyPrev));//(hx * (hx +
         hxPrev));
         al2[i - w] = -2.0 * _lambda / (hxPrev * (hx + hxPrev));
         au2[i] = -2.0 * _lambda / (hx * (hx + hxPrev));
         b[i] = _f(x, y);
         break;
      }
      case NodeType::DUMMY: {
         di[i] = 1.0;
         b[i] = 0.0;
         break;
      }
      }
   }
   return factory->createSystem(di, au1, au2, al1, al2, b);
}
vector<Node> IrregularGrid::nodes() {
   vector<Node> nodes = vector<Node>();
   int sumByX = int(_scale / _stepX);
   int sumByY = int(_scale / _stepY);
   int stepX = 1;
   int stepY;
   double eps = _stepX < _stepY ? _stepX * _epsBase : _stepY * _epsBase;
   for (int i = 0; i <= sumByX; i += stepX) {
      stepY = 1;
      for (int j = 0; j <= sumByY; j += stepY) {
         double x = _stepX * (double)i;
         double y = _stepY * (double)j;
         if (_figure->isEdgeNode(x, y, eps)) {
            nodes.push_back(Node(NodeType::EDGE, x, y));
         }
         else {
            if (_figure->isInnerNode(x, y, eps)) {
               nodes.push_back(Node(NodeType::INNER, x, y));
            }
            else {
               nodes.push_back(Node(NodeType::DUMMY, x, y));
            }
         }
         if (j >= sumByY / 2) {
            stepY = 2;
         }
      }
      if (i >= sumByX / 2) {
         stepX = 2;
      }
   }
   return nodes;
}
double IrregularGrid::leftDerivativeByX(double x, double y) {
   double h = 1e-9;
   return (_u(x, y) - _u(x - h, y)) / h;
}
double IrregularGrid::leftDerivativeByY(double x, double y) {
   double h = 1e-9;
   return (_u(x, y) - _u(x, y - h)) / h;
}
double IrregularGrid::rightDerivativeByX(double x, double y) {
   double h = 1e-9;
   return (_u(x + h, y) - _u(x, y)) / h;
}
double IrregularGrid::rightDerivativeByY(double x, double y) {
   double h = 1e-9;
   return (_u(x, y + h) - _u(x, y)) / h;
}
int IrregularGrid::width() {
   int sumByY = int(_scale / _stepY);
   int stepY = 1;
   int w = 0;
   for (int j = 0; j <= sumByY; j += stepY, w++) {
      if (j >= sumByY / 2) {
         stepY = 2;
      }
   }
   return w;
}

