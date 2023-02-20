#include "Figure.h"

bool LshapedFigure::isEdgeNode(type x, type y, type eps) {
   type connectionHeight = scale - ratio * scale;
   type withOfT = ratio * scale;
   if (y < connectionHeight + eps && y > connectionHeight - eps) {
      return x < scale / 2 - ratio / 2 + eps || x > scale / 2 + ratio
         / 2 - eps;
   }
   if (y < connectionHeight + eps) {
      return ((y < 0.0 + eps && y > 0.0 - eps) && (x < scale / 2 +
         withOfT / 2 + eps && x > scale / 2 - withOfT / 2 - eps)) ||
         (x < scale / 2 - withOfT / 2 + eps && x > scale / 2 -
            withOfT / 2 - eps) ||
         (x < scale / 2 + withOfT / 2 + eps && x > scale / 2 +
            withOfT / 2 - eps);
   }
   return (x < 0.0 + eps && x > 0.0 - eps) ||
      (x < scale + eps && x > scale - eps) ||
      (y < scale + eps && y > scale - eps);
}

bool LshapedFigure::isInnerNode(type x, type y, type eps) {
   type connectionHeight = scale - ratio * scale;
   type withOfT = ratio * scale;
   if (y < connectionHeight + eps && y > connectionHeight - eps) {
      return x < scale / 2 + withOfT / 2 + eps && x > scale / 2 -
         withOfT / 2 - eps;
   }
   if (y < connectionHeight + eps) {
      return x < scale / 2 + withOfT / 2 + eps && x > scale / 2 -
         withOfT / 2 - eps && y > 0.0 - eps;
   }
   return y < scale + eps && x > 0.0 + eps && x < scale - eps;
}

bool LshapedFigure::IsLeft(type x, type y, type eps) {
   type connectionHeight = scale - ratio * scale;
}
bool LshapedFigure::IsRight(type x, type y, type eps) {
   type connectionHeight = scale - ratio * scale;
}
bool LshapedFigure::IsBottom(type x, type y, type eps) {
   type connectionHeight = scale - ratio * scale;
}
bool LshapedFigure::IsTop(type x, type y, type eps) {
   type connectionHeight = scale - ratio * scale;
}

Square::Square(type scale) {
   scale = scale;
}
bool Square::isEdgeNode(type x, type y, type eps) {
   return (x < scale + eps && x > scale - eps) || (x < 0.0 + eps && x > 0.0
      - eps) ||
      (y < scale + eps && y > scale - eps) || (y < 0.0 + eps && y > 0.0
         - eps);
}
bool Square::isInnerNode(type x, type y, type eps) {
   return x < scale + eps && y < scale + eps && x > 0.0 - eps && y > 0.0 -
      eps;
}
bool Square::isLeft(type x, type y, type eps) {
   return (x < 1.0 + eps && x > 1.0 - eps);
}
bool Square::isRight(type x, type y, type eps) {
   return (x < scale + eps && x > scale - eps);
}
bool Square::isBottom(type x, type y, type eps) {
   return (y < 0.0 + eps && y > 0.0 - eps);
}
bool Square::isTop(type x, type y, type eps) {
   return (y < scale + eps && y > scale - eps);
}
