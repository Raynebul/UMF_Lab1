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

type MKR::norm(Vector& a) // норма корень квадратный из суммы квадратов
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
	case 1: // в соответствии с тестом
		return 10;

	default:
		cout << "Something in teta get wrong" << endl;
		break;
	}
	return 0;

}

type MKR::step(Vector l1, Vector l2, Vector u1, Vector u2, Vector& di, Vector& b, Vector& x0, Vector& x, type w, int& N, int& m, Vector& z, int i)
{
	double t = 0;
	t += di[i] * x0[i];// элемент главной диагонали

	// ниже
	if (i > 0)
		t += l1[i - 1] * z[i - 1];
	if (i > m + 1)
		t += l2[i - m - 2] * z[i - m - 2];

	// выше
	if (i < N - 1)
		t += u1[i] * x0[i + 1];
	if (i < N - m - 2)
		t += u2[i] * x0[i + m + 2];

	t = b[i] - t;
	return t;
}

// Работает нормально (проблема с диагональю - там короче x0[6] получается 0, и затем следующие x получаются nan(ind)
type MKR::GaussZeidel(Vector l1, Vector l2, Vector u1, Vector u2, Vector& di, Vector& f, Vector& x0,type& eps, type& w, int& N, int& m, int& max_iter)
{
	type delta = 1;  // относительная невязка
	for (int k = 0; k < max_iter || delta >= eps; k++)   // счетчик итераций
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

		delta = sqrt(delta) / norm(f);   // относительная невязка      
	}
	return delta;
}

type MKR::func(Node point)
{
	/*if (vectorAreas[point.area].fId == 1)
	{
		return 10;
	}
	else
	{
		return 0;
	}
	return 0;*/
	switch (vectorAreas[point.area].fId)
	{
	case 0:
		return 0;
	case 1:
		return 10; // меняем функцию
	case 2:
		return 10;
	default:
		cout << "Something in func get wrong" << endl;
		break;
	}
	return 0;
}

void MKR::Areas()
{
	ifstream in("Areas.txt");
	in >> countArea;
	vectorAreas.resize(countArea);

	// x1 x2 y1 y2 
	// область состоит из двух точек - (x1, y1) и (x2, y2) - левая нижняя и правая верхняя;
	for (int currArea = 0; currArea < countArea; currArea++)
	{
		//ввод точке, лямбды, гаммы и функции
		in >> vectorAreas[currArea].x1 >> vectorAreas[currArea].x2
			>> vectorAreas[currArea].y1 >> vectorAreas[currArea].y2
			>> vectorAreas[currArea].lambda >> vectorAreas[currArea].gamma >> vectorAreas[currArea].fId;

		// какие точки входят в область сетки - проверяется по < и >
		for (int b = 0; b < countY && vectorAreas[currArea].y2 >= vectorY[b]; b++)
			if (vectorAreas[currArea].y1 <= vectorY[b])
				for (int a = 0; a < countX && vectorAreas[currArea].x2 >= vectorX[a]; a++)
					if (vectorAreas[currArea].x1 <= vectorX[a])
						grid[b * countX + a].area = currArea; //область точки (x, y) в сетке

	}
	in.close();
}

void MKR::Grid()
{
	ifstream in("Grid.txt");
	vector<Area> vectorN;
	type X, Y; //точки x и y
	int Nx, Ny, kx, ky; //количество разбиений

	in >> countX >> countY;

	//vectorX.resize(countX); // все x
	//vectorY.resize(countY); // все 


	vectorX.resize(countX); // все x
	vectorY.resize(countY); // все 

	//in >> vectorX[0] >> vectorY[0];
	in >> vectorX[0] >> vectorY[0];

	// Сначала идут узлы x
	for (int currCountX = 0; currCountX < countX-1; )
	{
		//in >> X >> Nx >> kx;
		// x, кол-во разбиений, коэф. разрядки
		in >> X >> Nx >> kx;
		type hx;
		if (kx == 1) // равномерная сетка
		{
			hx = (X - vectorX[currCountX]) / Nx; // коэффициент приращения, получаем по формуле (1+1+1+...)
			for (int p = 1; p < Nx; p++)
				vectorX[currCountX + p] = vectorX[currCountX] + hx * p; //для нахождения последующих h
			currCountX += Nx; //идут через следующее разбиение
		}
		else // неравномерная сетка
		{
			hx = (X - vectorX[currCountX]) * (kx - 1) / (pow(kx, Nx + 1) - 1); //геометричесаяпрогрессия (1+ kx + k2x + k3x...)
			for (int p = 0; p < Nx - 1; currCountX++, p++)
				vectorX[currCountX + 1] = vectorX[currCountX] + hx * pow(kx, p); //для нахождения последующих h
			currCountX++;
		}

		vectorX[currCountX] = X;
	}

	// строим равномерную или неравномерную сетку
	for (int currCountY = 0; currCountY < countY-1; )
	{
		in >> Y >> Ny >> ky;
		type hy;
		if (ky == 1)
		{
			hy = (Y - vectorY[currCountY]) / Ny; // коэффициент приращения, получаем по формуле (1+1+1+...)
			for (int p = 1; p < Ny; p++)
				vectorY[currCountY + p] = vectorY[currCountY] + hy * p; //для нахождения последующих h
			currCountY += Ny;
		}
		else
		{
			hy = (Y - vectorY[currCountY]) * (ky - 1) / (pow(ky, Ny + 1) - 1); //геометрическая прогрессия (1+ kx + k2x + k3x...)
			for (int p = 0; p < Ny - 1; currCountY++, p++)
				vectorY[currCountY + 1] = vectorY[currCountY] + hy * pow(ky, p); //для нахождения последующих h
			currCountY++;
		}
		vectorY[currCountY] = Y;
	}
	in.close();

	// сетка класса Node();
	grid.resize(countY * countX);

	//двумерная матрица через одномерную ? может лучше просто двумерную ???
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
	// Учёт первых краевых - по границам
	for (int currArea = 0; currArea < countArea; currArea++)
	{
		in >> x1 >> x2 >> y1 >> y2 >> boundId;
		vectorB1[currArea].teta = boundId;
		vectorB1[currArea].nodes.reserve(countY * countX);

		// если линия горизонтальная
		if (y1 == y2)
		{
			int a, b;
			for (b = 0; b < countY && y1 != vectorY[b]; b++);
			for (a = 0; a < countX && x2 >= vectorX[a]; a++)
				if (x1 <= vectorX[a])
					vectorB1[currArea].nodes.push_back(b * countX + a); //убираем из областей краевые точки
		}
		else
			//если линия вертикальная
			if (x1 == x2)
			{
				int a, b;
				for (a = 0; a < countX && x1 != vectorX[a]; a++);
				for (b = 0; b < countY && y2 >= vectorY[b]; b++)
					if (y1 <= vectorY[b])
						vectorB1[currArea].nodes.push_back(b * countX + a); //убираем из областей краевые точки
			}
			else
				printf_s("Not bound");
	}
	in.close();

	// вертикальная (х - одинаковый, у - разный)
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

	// горизонтальная (у - одинаковый, х - разный)
	in.open("B2y.txt");
	in >> countArea;
	vectorB2y.resize(countArea);

	for (int currArea = 0; currArea < countArea; )
	{
		in >> x1 >> x2 >> y1 >> boundId;
		vectorB2y[currArea].teta = boundId;
		vectorB2y[currArea].nodes.reserve(countY * countX); //выделение памяти
		int a, b;
		for (b = 0; b < countY && y1 != vectorY[b]; );
		for (a = 0; a < countX && x2 >= vectorX[a]; )
			if (x1 <= vectorX[a])
				vectorB2y[currArea].nodes.push_back(b * countX + a);
	}

	in.close();

	Aij.k = countX - 2;

}

// Ввод всех данных из файлов
int MKR::input()
{
	// ввод данных по сетке    
	//Vector vectorX, vectorY; // вектор x и y
	// Сетка
	Grid();
	// Области
	Areas();
	// Краевые
	Boundaries();
	return 0;
}

int MKR::makeMatrix()
{

	// Диагональный формат - 5 диагоналей, каждый отдельный вектор
	Aij.di.resize(grid.size());
	Aij.u1.resize(grid.size() - 1);
	Aij.l1.resize(grid.size() - 1);
	Aij.u2.resize(grid.size() - 2 - Aij.k);
	Aij.l2.resize(grid.size() - 2 - Aij.k);
	f.resize(grid.size()); // наш b

	f[0] = func(grid[0]);

	// заполнить фиктивные узлы, если они по верхней/нижней границе
	for (int i = 0; i < Aij.k + 2; i++)
		//if (!vectorAreas[grid[i].area].lambda)
		{
			f[i] = 0;
			Aij.di[i] = 1;
		}

	for (int i = grid.size() - 2 - Aij.k; i < grid.size(); i++)
		//if (!vectorAreas[grid[i].area].lambda)
		{
			f[i] = 0;
			Aij.di[i] = 1;
		}

	// блок с 5 диаг (верх/ниж границ не смотрим, боковые границы поправим в граничных условиях)
	for (int s = Aij.k + 3; s < grid.size() - 2 - Aij.k; s++)
	{
		// hx=x(i)-x(i-1), hy=y(i)-y(i-1)
		type hx1 = grid[s].x - grid[s - 1].x;
		type hx2 = grid[s + 1].x - grid[s].x;
		type hy1 = grid[s].y - grid[s - Aij.k - 2].y;
		type hy2 = grid[s + Aij.k + 2].y - grid[s].y;
		type curLambda = vectorAreas[grid[s].area].lambda;

		// если разность отрицательна, взят боковой граничный узел как текущий
		//if (curLambda && !(s % Aij.k+1 == 0 || Aij.k + 1 == 1))
		if (!(s % Aij.k + 1 == 0 || Aij.k + 1 == 1))
		{
			if (hx1 > 0 && hx2 > 0)
			{
				//формула с оператором лапласа
				f[s] = func(grid[s]);
				Aij.di[s] += curLambda * 2 / (hx1 * hx2) + curLambda * 2 / (hy1 * hy2) + vectorAreas[grid[s].area].gamma;
				Aij.u1[s] += curLambda * -2 / (hx2 * (hx1 + hx2));
				Aij.l1[s - 1] += curLambda * -2 / (hx1 * (hx1 + hx2));
				Aij.u2[s] += curLambda * -2 / (hy2 * (hy1 + hy2));
				Aij.l2[s - Aij.k - 2] += curLambda * -2 / (hy1 * (hy1 + hy2));
			}
			else
			{
				f[s] = 0;
				Aij.di[s] = 1;
			}
		}
		else
		{
			f[s] = 0;
			Aij.di[s] = 1;
		}
	}
		//f[grid.size() - 1] = func(grid[grid.size() - 1]);

		// краевые 2
		//   // ОХ
		for (int i = 0; i < vectorB2x.size(); i++)
			for (int j = 0; j < vectorB2x[i].nodes.size(); j++)
			{
				int v = vectorB2x[i].nodes[j]; // глобальный номер текущего узла

				// занулить остальные элементы строки
				// проверить, есть ли эти элементы на строке
				if (v % (Aij.k + 2) != (Aij.k + 1) && vectorAreas[grid[v + 1].area].lambda != 0) // есть 																				узел правее не фиктивный
				{
					f[v] = -teta(grid[v], vectorB2x[i].teta);
					double h = grid[v + 1].x - grid[v].x;
					Aij.u1[v] = 1 / h;
					Aij.di[v] = -1 / h;

					if (v % (Aij.k + 2) != 0) // если и левее узел
						Aij.l1[v - 1] = 0;
				}
				else
					if (v % (Aij.k + 2) != 0) // если нет узла правее, но есть левее
					{
						f[v] = teta(grid[v], vectorB2x[i].teta);
						type h = grid[v].x - grid[v - 1].x;
						Aij.di[v] = 1 / h;
						Aij.l1[v - 1] = -1 / h;
					}

				if (v < grid.size() - Aij.k - 2)
					Aij.u2[v] = 0;
				if (v >= Aij.k + 2)
					Aij.l2[v - Aij.k - 2] = 0;
			}

		//   // ОУ
		for (int i = 0; i < vectorB2y.size(); i++)
		{
			for (int j = 0; j < vectorB2y[i].nodes.size(); j++)
			{
				int v = vectorB2y[i].nodes[j]; // глобальный номер текущего узла

				f[v] = teta(grid[v], vectorB2y[i].teta);

				// занулить остальные элементы строки
				// проверить, есть ли эти элементы на строке
				if (v % (Aij.k + 2) != (Aij.k + 1))
					Aij.u1[v] = 0;

				if (v % (Aij.k + 2) != 0)
					Aij.l1[v - 1] = 0;

				if (v < grid.size() - Aij.k - 2 && vectorAreas[grid[v + Aij.k + 2].area].lambda != 0) // 																			есть узел выше не фиктивный
				{
					type h = grid[v + Aij.k + 2].y - grid[v].y;
					Aij.u2[v] = 1 / h;
					Aij.di[v] = -1 / h;
					if (v >= Aij.k + 2)
						Aij.l2[v - Aij.k - 2] = 0;
				}
				else
					if (v >= Aij.k + 2)
					{
						type h = grid[v].y - grid[v - Aij.k - 2].y;
						Aij.l2[v - Aij.k - 2] = -1 / h;
						Aij.di[v] = 1 / h;
					}
			}
		}

		// краевые 1
		for (int i = 0; i < vectorB1.size(); i++)
		{
			for (int j = 0; j < vectorB1[i].nodes.size(); j++)
			{
				int v = vectorB1[i].nodes[j]; // глобальный номер текущего узла
				Aij.di[v] = 1; // 1 на главной
				f[v] = teta(grid[v], vectorB1[i].teta);

				// занулить остальные элементы строки
				// проверить, есть ли эти элементы на строке
				if (v % (Aij.k + 2) != (Aij.k + 1))
					Aij.u1[v] = 0;

				if (v % (Aij.k + 2) != 0)
					Aij.l1[v - 1] = 0;

				if (v < grid.size() - Aij.k - 2)
					Aij.u2[v] = 0;

				if (v >= Aij.k + 2)
					Aij.l2[v - Aij.k - 2] = 0;
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
		out << u[i] << endl;
	}

	return temp;
}


//////////////////////////////
///////БЛОЧНАЯ РЕЛАКСАЦИЯ/////
//////////////////////////////

/*
int MKR::BlockRelaxation()
{
	omega = 0.01;
	ofstream output("output.txt"); //файл, куда будет все записываться.
	//ReadFile();
	int blocksize = blocksize_;
	Factorize(blocksize); // факториизация матрицы
	while (omega < 2.01) //оптимальный параметр релаксации
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
	org_di = A[2]; //оставляем изначальные значения диагоналей для нахождение residual
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
///////БЛОЧНАЯ РЕЛАКСАЦИЯ/////
//////////////////////////////
