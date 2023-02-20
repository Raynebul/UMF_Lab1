#include "MKR.h"

MKR::MKR()
{

};

void MKR::GaussZeidel() 
{

};

void MKR::BlockRelaxation()
{

};

void MKR::Areas()
{
	ifstream in("Areas.txt");
	in >> countArea;
	vectorAreas.resize(countArea);

	// x1 x2 y1 y2 
	for (int currArea = 0; currArea < countArea; currArea++)
	{
		in >> vectorAreas[currArea].x1 >> vectorAreas[currArea].x2
			>> vectorAreas[currArea].y1 >> vectorAreas[currArea].y2
			>> vectorAreas[currArea].lambda >> vectorAreas[currArea].gamma >> vectorAreas[currArea].fId;

		for (int b = 0; b < countY && vectorAreas[currArea].y2 >= vectorY[b]; )
			if (vectorAreas[currArea].y1 <= vectorY[b])
				for (int a = 0; a < countX && vectorAreas[currArea].x2 >= vectorX[a]; )
					if (vectorAreas[currArea].x1 <= vectorX[a])
						grid[b * countX + a].area = currArea;

	}
	in.close();
}

void MKR::Grid()
{
	ifstream in("Grid.txt");
	vector<Area> vectorN;
	type X, Y; //точки x и y, коэффициент разр€дки?
	int Nx, Ny; //количество разбиений

	in >> countX >> countY;

	//vectorX.resize(countX); // все x
	//vectorY.resize(countY); // все 


	vectorX.resize(countX); // все x
	vectorY.resize(countY); // все 

	//in >> vectorX[0] >> vectorY[0];
	in >> vectorX[0] >> vectorY[0];

	for (int currCountX = 0; currCountX < countX - 1; )
	{
		//in >> X >> Nx >> kx;
		in >> X >> Nx >> kx[currCountX];
		double hx;
		if (kx[currCountX] == 1) // равномерна€ сетка
		{
			hx = (X - vectorX[currCountX]) / Nx; // коэффициент приращени€, получаем по формуле (1+1+1+...)
			for (int p = 1; p < Nx;)
				vectorX[currCountX + p] = vectorX[currCountX] + hx * p;
			currCountX += Nx;
		}
		else // неравномерна€ сетка
		{
			hx = (X - vectorX[currCountX]) * (kx[currCountX] - 1) / (pow(kx[currCountX], Nx + 1) - 1); //геометричеса€прогресси€ (1+ kx + k2x + k3x...)
			for (int p = 0; p < Nx - 1; currCountX++, p++)
				vectorX[currCountX + 1] = vectorX[currCountX] + hx * pow(kx[currCountX], p);
			currCountX++;
		}

		vectorX[currCountX] = X;
	}

	// строим равномерную или неравномерную сетку
	for (int currCountY = 0; currCountY < countY - 1; )
	{
		in >> Y >> Ny >> ky[currCountY];
		double hy;
		if (ky[currCountY] == 1)
		{
			hy = (Y - vectorY[currCountY]) / Ny; // коэффициент приращени€, получаем по формуле (1+1+1+...)
			for (int p = 1; p < Ny; p++)
				vectorY[currCountY + p] = vectorY[currCountY] + hy * p;
			currCountY += Ny;
		}
		else
		{
			hy = (Y - vectorY[currCountY]) * (ky[currCountY] - 1) / (pow(ky[currCountY], Ny + 1) - 1); //геометричеса€прогресси€ (1+ kx + k2x + k3x...)
			for (int p = 0; p < Ny - 1; currCountY++, p++)
				vectorY[currCountY + 1] = vectorY[currCountY] + hy * pow(ky[currCountY], p);
			currCountY++;
		}
		vectorY[currCountY] = Y;
	}
	in.close();

	// сетка класса Node();
	grid.resize(countY * countX);

	for (int i = 0; i < countY; i++)
		for (int j = 0; j < countX; j++)
		{
			grid[i * countX + j].x = vectorX[j];
			grid[i * countX + j].y = vectorY[i];
		}

}

void MKR::Boundaries()
{
	ifstream in("B1.txt");
	int boundId;
	in >> countArea;
	vectorB1.resize(countArea);
	type x1, x2, y1, y2;
	for (int currArea = 0; currArea < countArea; currArea++)
	{
		in >> x1 >> x2 >> y1 >> y2 >> boundId;
		vectorB1[currArea].teta = boundId;
		vectorB1[currArea].nodes.reserve(countY * countX);

		// если лини€ горизонтальна€
		if (y1 == y2)
		{
			int a, b;
			for (b = 0; b < countY && y1 != vectorY[b]; );
			for (a = 0; a < countX && x2 >= vectorX[a]; )
				if (x1 <= vectorX[a])
					vectorB1[currArea].nodes.push_back(b * countX + a);
		}
		else
			//если лини€ вертикальна€
			if (x1 == x2)
			{
				int a, b;
				for (a = 0; a < countX && x1 != vectorX[a]; );
				for (b = 0; b < countY && y2 >= vectorY[b]; )
					if (y1 <= vectorY[b])
						vectorB1[currArea].nodes.push_back(b * countX + a);
			}
			else
				printf_s("Not bound");
	}
	in.close();

	// вертикальна€ (х - одинаковый, у - разный)
	in.open("B2x.txt");
	in >> countArea;
	vectorB2x.resize(countArea);

	for (int currArea = 0; currArea < countArea; )
	{
		in >> x1 >> y1 >> y2 >> boundId;
		vectorB2x[currArea].teta = boundId;
		vectorB2x[currArea].nodes.reserve(countY * countX);

		int a, b;
		for (a = 0; a < countX && x1 != vectorX[a]; );
		for (b = 0; b < countY && y2 >= vectorY[b]; )
			if (y1 <= vectorY[b])
				vectorB2x[currArea].nodes.push_back(b * countX + a);
	}
	in.close();

	// горизонтальна€ (у - одинаковый, х - разный)
	in.open("B2y.txt");
	in >> countArea;
	vectorB2y.resize(countArea);

	for (int currArea = 0; currArea < countArea; )
	{
		in >> x1 >> x2 >> y1 >> boundId;
		vectorB2y[currArea].teta = boundId;
		vectorB2y[currArea].nodes.reserve(countY * countX);
		int a, b;
		for (b = 0; b < countY && y1 != vectorY[b]; );
		for (a = 0; a < countX && x2 >= vectorX[a]; )
			if (x1 <= vectorX[a])
				vectorB2y[currArea].nodes.push_back(b * countX + a);
	}

	in.close();

	Aij.k = countX - 2;

}

int MKR::input()
{
	// ввод данных по сетке    
	//Vector vectorX, vectorY; // вектор x и y
	Grid();
	Areas();
	Boundaries();
	return 0;
}
