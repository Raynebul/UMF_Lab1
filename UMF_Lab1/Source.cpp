#include "MKR.h"

using namespace std;

int main() {
   string filename = "Result.txt";
   MKR Object;
   Object.input();
   Object.makeMatrix();
   Object.solve(filename);
   return 0;
}