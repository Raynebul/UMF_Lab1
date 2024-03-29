#include "MKR.h"

MKR::MKR()
{

};

int MKR::sub(Vector& a, Vector& b) // b = a - b
{
	for (int i = 0; i < a.size(); )
		b[i] = a[i] - b[i];
	return 0;
}

type MKR::norm(Vector& a) // ����� ������ ���������� �� ����� ���������
{
	type res = 0;
	for (int i = 0; i < a.size(); i++)
		res += a[i] * a[i];
	res = sqrt(res);
	return res;
}

type MKR::teta(Node point, int idBound)
{
	switch (idBound)
	{
	case 0:
		return 0;
	case 1:
		return 0; // ������ �������
	case 2:
		return 1;
	case 3:
		return 2 * point.x;
	case 4:
		return 3 * pow(point.x, 2);
	case 5:
		return pow(point.x, 4) + pow(point.y, 4);
	case 6:
		return point.x * point.y;
	case 7:
		return point.x * point.x * point.y;
	case 8:
		return cos(point.x);
	default:
		cout << "Something in teta get wrong" << endl;
		break;
	}
	return 0;

}

type MKR::tetaX(Node point, int idBound)
{
	switch (idBound)
	{
	case 0:
		return 0;
	case 1:
		return 0; // ������ �������
	case 2:
		return 1;
	case 3:
		return 2 * point.x;
	case 4:
		return 3* pow(point.x, 2);
	case 5:
		return 4*pow(point.x, 3);
	case 6:
		return point.y;
	case 7:
		return 2* point.x * point.y;
	case 8:
		return -sin(point.x);
	default:
		cout << "Something in teta get wrong" << endl;
		break;
	}
	return 0;

}

type MKR::tetaY(Node point, int idBound)
{
	switch (idBound)
	{
	case 0:
		return 0;
	case 1:
		return 0; // ������ �������
	case 2:
		return 1;
	case 3:
		return 2 * point.y;
	case 4:
		return 3*pow(point.y, 2);
	case 5:
		return 4* pow(point.y, 3);
	case 6:
		return point.x;
	case 7:
		return point.x * point.x;
	case 8:
		return 0;
	default:
		cout << "Something in teta get wrong" << endl;
		break;
	}
	return 0;

}

type MKR::func(Node point)
{
	switch (vectorAreas[point.area].fId)
	{
	case 0:
		return 0;
	case 1:
		return 10; // ������ �������
	case 2:
		return point.x + point.y;
	case 3:
		return point.x*point.x+point.y*point.y-4;
	case 4:
		return pow(point.x, 3) + pow(point.y, 3) -6*point.x - 6*point.y;
	case 5:
		return pow(point.x, 4) + pow(point.y, 4) - 12 * point.x*point.x - 12 * point.y * point.y;
	case 6:
		return point.x * point.y;
	case 7:
		return point.x* point.x * point.y - 2*point.y;
	case 8:
		return 2*cos(point.x);
	default:
		cout << "Something in func get wrong" << endl;
		break;
	}
	return 0;
}

type MKR::u(Node point, int idBound)
{
	switch (idBound)
	{
	case 0:
		return 0;
	case 1:
		return 10; // ������ �������
	case 2:
		return point.x + point.y;
	case 3:
		return point.x * point.x + point.y * point.y;
	case 4:
		return pow(point.x, 3) + pow(point.y, 3);
	case 5:
		return pow(point.x, 4) + pow(point.y, 4);
	case 6:
		return point.x*point.y;
	case 7:
		return point.x* point.x * point.y;
	case 8:
		return cos(point.x);
	default:
		cout << "Something in func get wrong" << endl;
		break;
	}
	return 0;
}

type MKR::step(Vector l1, Vector l2, Vector u1, Vector u2, Vector& di, Vector& b, Vector& x0, Vector& x, type w, int& N, int& m, Vector& z, int i)
{
	double t = 0;
	t += di[i] * x0[i];// ������� ������� ���������

	// ����
	if (i > 0)
		t += l1[i-1] * z[i-1];
	if (i > m + 1)
		t += l2[i-m-2] * z[i-m-2];

	// ����
	if (i < N - 1)
		t += u1[i] * x0[i+1];
	if (i < N - m - 2)
		t += u2[i] * x0[i+m+2];

	t = b[i] - t;
	return t;
}

// �������� ��������� (�������� � ���������� - ��� ������ x0[6] ���������� 0, � ����� ��������� x ���������� nan(ind)
type MKR::GaussZeidel(Vector l1, Vector l2, Vector u1, Vector u2, Vector& di, Vector& f, Vector& x0,type& eps, type& w, int& N, int& m, int& max_iter)
{
	type delta = 1;  // ������������� �������
	for (int k = 0; k < max_iter || delta >= eps; k++)   // ������� ��������
	{
		delta = 0;
		for (int i = 0; i < N; i++)
		{
			type t = step(l1, l2, u1, u2, di, f, x0, x0, w, N, m, x0, i);
			x0[i] = x0[i] + w * t / di[i];
		}
		for (int i = 0; i < N; i++)
		{
			type t = step(l1, l2, u1, u2, di, f, x0, x0, w, N, m, x0, i);
			delta += t * t;
		}

		delta = sqrt(delta) / norm(f);   // ������������� �������      
	}
	return delta;
}

void MKR::Areas()
{
	ifstream in("Areas.txt");
	in >> countArea;
	vectorAreas.resize(countArea);

	// x1 x2 y1 y2 
	// ������� ������� �� ���� ����� - (x1, y1) � (x2, y2) - ����� ������ � ������ �������;
	for (int currArea = 0; currArea < countArea; currArea++)
	{
		//���� �����, ������, ����� � �������
		in >> vectorAreas[currArea].x1 >> vectorAreas[currArea].x2
			>> vectorAreas[currArea].y1 >> vectorAreas[currArea].y2
			>> vectorAreas[currArea].lambda >> vectorAreas[currArea].gamma >> vectorAreas[currArea].fId;

		// ����� ����� ������ � ������� ����� - ����������� �� < � >
		for (int b = 0; b < countY && vectorAreas[currArea].y2 >= vectorY[b]; b++)
			if (vectorAreas[currArea].y1 <= vectorY[b])
				for (int a = 0; a < countX && vectorAreas[currArea].x2 >= vectorX[a]; a++)
					if (vectorAreas[currArea].x1 <= vectorX[a])
						grid[b * countX + a].area = currArea; //������� ����� (x, y) � �����

	}
	in.close();
}

void MKR::Grid()
{
	ifstream in("Grid.txt");
	vector<Area> vectorN;
	type X, Y; //����� x � y
	int Nx, Ny; //���������� ���������
	type kx, ky;
	in >> countX >> countY;

	//vectorX.resize(countX); // ��� x
	//vectorY.resize(countY); // ��� 


	vectorX.resize(countX); // ��� x
	vectorY.resize(countY); // ��� 

	//in >> vectorX[0] >> vectorY[0];
	in >> vectorX[0] >> vectorY[0];

	// ������� ���� ���� x
	for (int currCountX = 0; currCountX < countX-1; )
	{
		//in >> X >> Nx >> kx;
		// x, ���-�� ���������, ����. ��������
		in >> X >> Nx >> kx;
		type hx;
		if (kx == 1) // ����������� �����
		{
			hx = (X - vectorX[currCountX]) / Nx; // ����������� ����������, �������� �� ������� (1+1+1+...)
			for (int p = 1; p < Nx; p++)
				vectorX[currCountX + p] = vectorX[currCountX] + hx * p; //��� ���������� ����������� h
			currCountX += Nx; //���� ����� ��������� ���������
		}
		else // ������������� �����
		{
			// = sqrt(kx);
			//kx
			hx = (X - vectorX[currCountX]) * (kx - 1) / (pow(kx, Nx) - 1); //����������������������� (1+ kx + k2x + k3x...)
			for (int p = 0; p < Nx; p++)
				vectorX[currCountX + 1 + p] = vectorX[currCountX] + hx * pow(kx, p); //��� ���������� ����������� h
			//vectorX[currCountX + 1] = (vectorX[currCountX] + vectorX[currCountX + 2]) / 2;
			//vectorX[currCountX + 3] = (vectorX[currCountX+2] + vectorX[currCountX + 4]) / 2;
			//vectorX[currCountX + 5] = (vectorX[currCountX+4] + X) / 2;
			currCountX+=Nx+1;
			//currCountX = 6;
		}

		vectorX[currCountX] = X;
	}

	// ������ ����������� ��� ������������� �����
	for (int currCountY = 0; currCountY < countY-1; )
	{
		in >> Y >> Ny >> ky;
		type hy;
		if (ky == 1)
		{
			hy = (Y - vectorY[currCountY]) / Ny; // ����������� ����������, �������� �� ������� (1+1+1+...)
			for (int p = 1; p < Ny; p++)
				vectorY[currCountY + p] = vectorY[currCountY] + hy * p; //��� ���������� ����������� h
			currCountY += Ny;
		}
		else
		{
			hy = (Y - vectorY[currCountY]) * (ky - 1) / (pow(ky, Ny + 1) - 1); //�������������� ���������� (1+ kx + k2x + k3x...)
			for (int p = 0; p < Ny - 1; currCountY++, p++)
				vectorY[currCountY + 1] = vectorY[currCountY] + hy * pow(ky, p); //��� ���������� ����������� h
			currCountY++;
		}
		vectorY[currCountY] = Y;
	}
	in.close();

	// ����� ������ Node();
	grid.resize(countY * countX);

	//��������� ������� ����� ���������� ? ����� ����� ������ ��������� ???
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
	// ���� ������ ������� - �� ��������
	for (int currArea = 0; currArea < countArea; currArea++)
	{
		in >> x1 >> x2 >> y1 >> y2 >> boundId;
		vectorB1[currArea].teta = boundId;
		vectorB1[currArea].nodes.reserve(countY * countX);

		// ���� ����� ��������������
		if (y1 == y2)
		{
			int a, b;
			for (b = 0; b < countY && y1 != vectorY[b]; b++);
			for (a = 0; a < countX && x2 >= vectorX[a]; a++)
				if (x1 <= vectorX[a])
					vectorB1[currArea].nodes.push_back(b * countX + a); //������� �� �������� ������� �����
		}
		else
			//���� ����� ������������
			if (x1 == x2)
			{
				int a, b;
				for (a = 0; a < countX && x1 != vectorX[a]; a++);
				for (b = 0; b < countY && y2 >= vectorY[b]; b++)
					if (y1 <= vectorY[b])
						vectorB1[currArea].nodes.push_back(b * countX + a); //������� �� �������� ������� �����
			}
			else
				printf_s("Not bound");
	}
	in.close();

	// ������������ (� - ����������, � - ������)
	in.open("B2x.txt");
	in >> countArea;
	vectorB2x.resize(countArea);

	for (int currArea = 0; currArea < countArea; currArea++)
	{
		in >> x1 >> y1 >> y2 >> boundId;
		vectorB2x[currArea].teta = boundId;
		vectorB2x[currArea].nodes.reserve(countY * countX);

		int a, b;
		for (a = 0; a < countX && x1 != vectorX[a]; a++);
		for (b = 0; b < countY && y2 >= vectorY[b]; b++)
			if (y1 <= vectorY[b])
				vectorB2x[currArea].nodes.push_back(b * countX + a);
	}
	in.close();

	// �������������� (� - ����������, � - ������)
	in.open("B2y.txt");
	in >> countArea;
	vectorB2y.resize(countArea);

	for (int currArea = 0; currArea < countArea; currArea++)
	{
		in >> x1 >> x2 >> y1 >> boundId;
		vectorB2y[currArea].teta = boundId;
		vectorB2y[currArea].nodes.reserve(countY * countX); //��������� ������
		int a, b;
		for (b = 0; b < countY && y1 != vectorY[b]; b++);
		for (a = 0; a < countX && x2 >= vectorX[a]; a++)
			if (x1 <= vectorX[a])
				vectorB2y[currArea].nodes.push_back(b * countX + a);
	}

	in.close();

	Aij.k = countX - 2;

}

// ���� ���� ������ �� ������
int MKR::input()
{
	// ���� ������ �� �����    
	//Vector vectorX, vectorY; // ������ x � y
	// �����
	Grid();
	// �������
	Areas();
	// �������
	Boundaries();
	return 0;
}

int MKR::makeMatrix()
{

	// ������������ ������ - 5 ����������, ������ ��������� ������
	Aij.di.resize(grid.size());
	Aij.u1.resize(grid.size() - 1);
	Aij.l1.resize(grid.size() - 1);
	Aij.u2.resize(grid.size() - countX); // k = CountX-2
	Aij.l2.resize(grid.size() - countX);
	f.resize(grid.size()); // ��� b

	f[0] = func(grid[0]);

	// ��������� ��������� ����, ���� ��� �� �������/������ �������
	for (int i = 0; i < countX+1; i++)
		//if (!vectorAreas[grid[i].area].lambda)
		{
			f[i] = 0;
			Aij.di[i] = 1;
		}

	for (int i = grid.size() - countX-1; i < grid.size(); i++)
		//if (!vectorAreas[grid[i].area].lambda)
		{
			f[i] = 0;
			Aij.di[i] = 1;
		}

	// ���� � 5 ���� (����/��� ������ �� �������, ������� ������� �������� � ��������� ��������)
	for (int s = countX +1; s < grid.size() - 1 - countX; s++)
	{
		// hx=x(i)-x(i-1), hy=y(i)-y(i-1)
		type hx1 = grid[s].x - grid[s - 1].x;
		type hx2 = grid[s + 1].x - grid[s].x;
		type hy1 = grid[s].y - grid[s - Aij.k - 3].y;
		type hy2 = grid[s + countX].y - grid[s].y;
		type curLambda = vectorAreas[grid[s].area].lambda;

		// ���� �������� ������������, ���� ������� ��������� ���� ��� �������
		//if (curLambda && !(s % Aij.k+1 == 0 || Aij.k + 1 == 1))
		if (grid[s].area != 0)
		{
			if (hx1 > 0 && hx2 > 0)
			{
				//������� � ���������� �������
				f[s] = func(grid[s]);
				Aij.di[s] += curLambda * 2 / (hx1 * hx2) + curLambda * 2 / (hy1 * hy2) + vectorAreas[grid[s].area].gamma;
				Aij.u1[s] += curLambda * -2 / (hx2 * (hx1 + hx2));
				Aij.l1[s - 1] += curLambda * -2 / (hx1 * (hx1 + hx2));
				Aij.u2[s] += curLambda * -2 / (hy2 * (hy1 + hy2));
				Aij.l2[s - countX] += curLambda * -2 / (hy1 * (hy1 + hy2));
			}
			else
			{
				f[s] = 0;
				Aij.di[s] = 1;
				Aij.u1[s] = 0;
				Aij.l1[s-1] = 0;
				Aij.u2[s] = 0;
				Aij.l2[s - countX] = 0;
			}
		}
		else
		{
			f[s] = 0;
			Aij.di[s] = 1;
			Aij.u1[s] = 0;
			Aij.l1[s - 1] = 0;
			Aij.u2[s] = 0;
			Aij.l2[s - countX] = 0;
		}
	}
		//f[grid.size() - 1] = func(grid[grid.size() - 1]);

		// ������� 2
			//   // ��
	for (int i = 0; i < vectorB2y.size(); i++)
	{
		for (int j = 0; j < vectorB2y[i].nodes.size(); j++)
		{
			int v = vectorB2y[i].nodes[j]; // ���������� ����� �������� ����

			f[v] = tetaY(grid[v], vectorB2y[i].teta);

			// �������� ��������� �������� ������
			// ���������, ���� �� ��� �������� �� ������
			if ((v + 1) % countX != 0)
			{
				Aij.u1[v] = 0; 
			}
			if ((v + 1) % countX != 1)
			{
				Aij.l1[v - 1] = 0; 
			}
			if (v < grid.size() - countX && vectorAreas[grid[v + countX].area].lambda != 0) // 																			���� ���� ���� �� ���������
			{
				type h = grid[v + countX].y - grid[v].y;
			//	f[v] = -tetaY(grid[v], vectorB2y[i].teta);
				Aij.u2[v] = 1 / h;
				Aij.di[v] = -1 / h;
				if (v >= countX)
				{
					//type h = grid[v].y - grid[v - countX].y;
					//Aij.l2[v - countX] = -1 / h;
					Aij.l2[v - countX] = 0;
				}

			}
			else
				if (v >= countX)
				{
				//	f[v] = tetaY(grid[v], vectorB2y[i].teta);
					type h = grid[v].y - grid[v - countX].y;
					Aij.l2[v - countX] = -1 / h;
					Aij.di[v] = 1 / h;
				}
			if (v + countX < grid.size())
			{
				if (grid[v + countX].area == 0 && grid[v + countX].y != countY)
				{
					type h = grid[v+countX].y - grid[v].y;
					Aij.u2[v] = -1 / h;
					Aij.di[v] = 1 / h;
				}
			}
			if (v >= countX)
			{
				if (grid[v - countX].area == 0 && grid[v - countX].y != grid[0].y)
				{
					type h = grid[v].y - grid[v - countX].y;
					Aij.l2[v - countX] = -1 / h;
					Aij.di[v] = 1 / h;
				}
			}
		}
	}

		//   // ��
		for (int i = 0; i < vectorB2x.size(); i++)
			for (int j = 0; j < vectorB2x[i].nodes.size(); j++)
			{
				int v = vectorB2x[i].nodes[j]; // ���������� ����� �������� ����

				// �������� ��������� �������� ������
				// ���������, ���� �� ��� �������� �� ������
				if ((v+1) % countX == 1 && vectorAreas[grid[v + 1].area].lambda != 0) // ���� 																				���� ������ �� ���������
				{
					f[v] = -tetaX(grid[v], vectorB2x[i].teta);
					double h = grid[v+1].x - grid[v].x;
					Aij.u1[v] = 1 / h; 
					Aij.di[v] = -1 / h;
					
					if ((v + 1) % countX != 0)
					{  // ���� � ����� ����
						//type h = grid[v].x - grid[v - 1].x;
						//Aij.l1[v - 1] = -1 / h; 
						Aij.l1[v - 1] = 0;
					}
				}
				else
					if ((v+1) % countX  == 0) // ���� ��� ���� ������, �� ���� �����
					{
						f[v] = tetaX(grid[v], vectorB2x[i].teta);
						type h = grid[v].x - grid[v-1].x;
						Aij.di[v] = 1 / h;
						Aij.l1[v - 1] = -1 / h; 
						//Aij.u1[v] = 0;
					}

				if (v < grid.size() - countX)
				{
					Aij.u2[v] = 0;
				}
					
				if (v >= countX)
				{
					Aij.l2[v - countX] = 0;
				}
				if (grid[v].x < grid[countX-1].x)
				{
					if (grid[v + 1].area == 0 && grid[v].x != grid[0].x)
					{
						type h = grid[v].y - grid[v - 1].y;
						Aij.l1[v - 1] = -1 / h;
						Aij.di[v] = 1 / h;
					}
				}
				if (grid[v].x > grid[0].x)
				{
					if (grid[v - 1].area == 0 && grid[v].x != grid[0].x)
					{
						type h = grid[v+1].y - grid[v].y;
						Aij.u1[v] = -1 / h;
						Aij.di[v] = 1 / h;
					}
				}

			}
	//	Aij.di[13] = 1;
	//	Aij.l2[8] = -1;
	  // f[13] = tetaY(grid[13], grid[13].area);

	//	Aij.di[41] = 2;
	//	Aij.l2[36] = -2;
	//	f[41] = tetaY(grid[41], 8);

	//	Aij.di[42] = 2;
	//	Aij.l2[37] = -2;
	//	f[42] = tetaY(grid[42], 8);

		//Aij.di[43] = 2;
	//	Aij.l2[38] = -2;
	//	f[43] = tetaY(grid[43], 8);

	 //  Aij.di[19] = 1;
	//	Aij.l2[12] = -1;
		//f[19] = tetaX(grid[19], grid[19].area);

		// ������� 1
		for (int i = 0; i < vectorB1.size(); i++)
		{
			for (int j = 0; j < vectorB1[i].nodes.size(); j++)
			{
				int v = vectorB1[i].nodes[j]; // ���������� ����� �������� ����
				Aij.di[v] = 1; // 1 �� �������
				f[v] = u(grid[v], vectorB1[i].teta);

				// �������� ��������� �������� ������
				// ���������, ���� �� ��� �������� �� ������
				if (v % (countX) != (Aij.k + 1))
					Aij.u1[v] = 0;

				if (v % countX != 0)
					Aij.l1[v - 1] = 0;

				if (v < grid.size() - countX)
					Aij.u2[v] = 0;

				if (v >= countX)
					Aij.l2[v - countX] = 0;
			}
		}
		return 0;
}

int MKR::solve(string filename)
{
	int N = grid.size(), m = Aij.k, max_iter = 100;
	type W = 1, eps = 1e-15;
	Vector u(N);

	type temp = GaussZeidel(Aij.l1, Aij.l2, Aij.u1, Aij.u2, Aij.di, f, u, eps, W, N, m, max_iter);
	ofstream out(filename);
	for (int i = 0; i < u.size(); i++)
	{
		out << setprecision(21) << u[i] << endl;
	}

	return temp;
}


//////////////////////////////
///////������� ����������/////
//////////////////////////////

/*
int MKR::BlockRelaxation()
{
	omega = 0.01;
	ofstream output("output.txt"); //����, ���� ����� ��� ������������.
	//ReadFile();
	int blocksize = blocksize_;
	Factorize(blocksize); // ������������� �������
	while (omega < 2.01) //����������� �������� ����������
	{
		int k;
		residual = epsilon + 1;
		x_next = x_start;
		clock_t begin = clock();
		for (k = 0; k<maxiter && residual > epsilon; k++)
		{
			BlockRelazationIteration(blocksize);
		}
		clock_t end = clock();
		clock_t time = end - begin;
		if (k < k_min)
		{
			k_min = k;
			residual_min = residual;
			opt_omega = omega;
			clockf = time;
		}
		WriteFile(k, output);
		omega += 0.01;
	}
	output.close();
	return 0;
}

int  MKR::BlockSize()
{
	// int block = m + 2;
	if (blocksize_ > m + 2 || n % blocksize_ != 0)
		return 0;
	if (blocksize_ == 1)
		return 1;
	return 2;
	// while (block > 1 && n % block != 0)
		// block--;
	 //return block;
}

void  MKR::Factorize(int blocksize)
{
	org_di = A[2]; //��������� ����������� �������� ���������� ��� ���������� residual
	org_u = A[3];
	org_l = A[1];
	int num = n / blocksize;
	int u = 3;
	int l = 1;

	for (int i = 0; i < num; i++)
	{
		for (int j = i * blocksize; j < (i + 1) * blocksize; j++)
		{
			if (j - 1 >= i * blocksize)
			{
				A[2][j] -= A[l][j] * A[u][j - 1];
			}
			if (j + 1 < (i + 1) * blocksize)
			{
				A[l][j + 1] /= A[2][j];
			}
		}
	}
};

void MKR::BlockRelazationIteration(int blocksize)
{
	type nev = 0.0;
	//type sumnev = 0.0;
	residual = 0.0;
	int num = n / blocksize;
	for (int i = 0; i < num; i++)
	{
		Vector Ri = Ri_(blocksize, i);
		for (size_t i = 0; i < blocksize; i++)
		{
			Ri[i] *= omega;
		}

		Vector solution = LU_SOL(i, blocksize, Ri);
		for (int j = 0; j < blocksize; j++)
		{
			x_next[i * blocksize + j] += solution[j];
		}
	}
	for (int i = 0; i < n; i++)
	{
		residual += ResidualBlock(i) * ResidualBlock(i);
	}
	residual = sqrt(residual) / sqrt(norm);
	int k = 0;

}

type MKR::ResidualBlock(int i)
{
	int index[5] = { -m - 2, -1, 0, 1, m + 2 };
	type sum = 0.0;
	sum += org_di[i] * x_next[i];
	if (i + index[0] >= 0 && i + index[0] < n)
		sum += A[0][i] * x_next[i + index[0]];
	if (i + index[1] >= 0 && i + index[1] < n)
		sum += org_l[i] * x_next[i + index[1]];
	if (i + index[3] >= 0 && i + index[3] < n)
		sum += org_u[i] * x_next[i + index[3]];
	if (i + index[4] >= 0 && i + index[4] < n)
		sum += A[4][i] * x_next[i + index[4]];
	return F[i] - sum;
}

Vector MKR::Ri_(int blocksize, int i)
{
	int num = n / blocksize;
	Vector Ri(blocksize);
	for (int j = 0; j < blocksize; j++)
	{
		Ri[j] = F[i * blocksize + j];
	}
	for (int j = 0; j < num; j++)
	{
		int start_col = j * blocksize;
		int end_col = (j + 1) * blocksize;
		int start_row = i * blocksize;
		int end_row = (i + 1) * blocksize;
		for (int a1 = start_row; a1 < end_row; a1++)
		{
			sum sum = 0.0;
			for (int a2 = start_col; a2 < end_col; a2++)
			{
				sum += FindAElement(a1, a2) * x_next[a2];
			}
			Ri[a1 % blocksize] -= sum;
			// residual += (F[a1]-sum) * (F[a1] - sum);
		}
	}
	return Ri;
};

type MKR::FindAElement(int a1, int a2)
{
	int index[5] = { -m - 2, -1, 0, 1, m + 2 };

	if (a2 == a1 + index[0])
		return A[0][a1];
	if (a2 == a1 + index[1])
		return org_l[a1];
	if (a2 == a1)
		return org_di[a1];
	if (a2 == a1 + index[3])
		return org_u[a1];
	if (a2 == a1 + index[4])
		return A[4][a1];
	return 0;
};

Vector MKR::LU_SOL(int i, int blocksize, Vector& Ri)
{
	Vector y(blocksize);
	int u = 3;
	int l = 1;
	y[0] = Ri[0];
	for (int counter = 1; counter < blocksize; counter++) //Ly=Ri
	{
		int new_offset = i * blocksize + counter;
		y[counter] = (Ri[counter] - A[l][i * blocksize + counter] * y[counter - 1]); // / A[2][i * blocksize + counter];
	}


	Vector x(blocksize);
	x[blocksize - 1] = y[blocksize - 1] / A[2][(i + 1) * blocksize - 1];
	for (int counter = blocksize - 2; counter >= 0; counter--) //UYi=y
	{
		x[counter] = (y[counter] - A[u][i * blocksize + counter] * x[counter + 1]) / A[2][i * blocksize + counter];
	}
	return x;
}
*/

//////////////////////////////
///////������� ����������/////
//////////////////////////////
