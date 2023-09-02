#include <iostream>
#include <fstream>
#include <ctime>
#include "mpi.h"

using namespace std;

/*Вариант 1. Классический вариант КА:
Поле – двумерная матрица, замкнутая в тор; 4 соседа у клетки;
состояние автомата – бинарное (0-мертв, 1-жив);
S(t) – исходное состояние; N(t) – количество «живых» соседей.
Описание переходов: S(t)=1, N(t)<2 –> S(t+1)=0;
					N(t)<4 –> S(t+1)=1;
					N(t)=4 –> S(t+1)=0;
					S(t)=0, N(t)=2 –> S(t+1)=1.*/

void WithOutput(int** array_first, int** array_second, int& size, int& iterations) {
	int neighbors = 0;
	int** tmp = NULL;
	bool noChanges = false;
	ofstream out;
	out.open("out.txt");

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			out << array_first[i][j] << ' ';
		}
		out << endl;
	}
	out << "\n\n\n\n\n";
	do {
		out << "Iterations: " << iterations << endl;
		noChanges = true;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				neighbors += array_first[i][(j + 1) % size];//правый
				neighbors += array_first[(i + 1) % size][j];//нижний
				neighbors += array_first[i][(j - 1 + size) % size];//левый
				neighbors += array_first[(i - 1 + size) % size][j];//верхний
				/*S(t) – исходное состояние; N(t) – количество «живых» соседей.
				Описание переходов: S(t) = 1, N(t) < 2 – > S(t + 1) = 0;
									N(t) < 4 – > S(t + 1) = 1;
									N(t) = 4 – > S(t + 1) = 0;
									S(t) = 0, N(t) = 2 – > S(t + 1) = 1.*/
				array_second[i][j] = (neighbors == 2) + array_first[i][j] * (neighbors == 3);
				neighbors = 0;
				noChanges = array_second[i][j] == array_first[i][j] && noChanges;
				out << array_second[i][j] << ' ';
			}
			out << endl;
		}
		out << "\n\n\n\n\n";
		tmp = array_first;
		array_first = array_second;
		array_second = tmp;
		iterations--;
	} while (iterations != 0 && !noChanges);
	out.close();
}

void WithoutOutputParal(int** array_first, int**& array_second, const int& size, int& iterations, const int& myrank, const int& ranksize) {
	clock_t delay;
	delay = clock();

	int max_size = ranksize > size ? size : ranksize;//ограничивает максимальное число процессов
	int size_r_tmp = size / max_size + (size % max_size != 0);//размер блока для процесса
	int workers = size / size_r_tmp + (size % size_r_tmp != 0);//определяем количество работающих процессов
	int be_work = myrank < workers;//переменная состояния процесса (рабочий/нерабочий)
	MPI_Comm comm;
	MPI_Comm_split(MPI_COMM_WORLD, be_work, myrank, &comm);//для каждого процесса создаётся коммутатор группы

	if (be_work) {
		//очищение второго массива
		if (array_second) {
			for (int i = 0; i < size; ++i) {
				delete array_second[i];
			}
			delete array_second;
		}

		int noChanges;
		int neighbors = 0;
		int* tmp = NULL, * r_tmp = NULL;
		MPI_Request* req_l = NULL, * req_r = NULL;
		MPI_Status* status = NULL;

		if (!myrank) {
			tmp = new int[workers * size_r_tmp * size];//буфер для отправки блоков каждому процессу

			for (int i = 0; i < workers; ++i)//цикл для каждого процесса
				for (int j = 0; j < size_r_tmp; ++j)//цикл по строкам текущего процесса
					for (int k = 0; k < size; ++k)//цикл по столбцам текущего процесса
						//заполнение буфера матрицей, уложенной в строку. Число -1 вносится в случае недостатка столбцов для последнего процесса (например, при разбиении все процессы получили по 3 строки, а последний процесс - 1 строку)
						tmp[i * size_r_tmp * size + j * size + k] = (i * size_r_tmp + j) < size ? array_first[i * size_r_tmp + j][k] : -1;
		}

		r_tmp = new int[size_r_tmp * size];//буфер в каждом процессе для принятия блока данных
		MPI_Scatter(tmp, size_r_tmp * size, MPI_INT, r_tmp, size_r_tmp * size, MPI_INT, 0, comm);//рассылка блоков матрицы каждому процессу
		if (tmp) delete tmp;

		//очищение первого массива
		if (array_first) {
			for (int i = 0; i < size; ++i) delete array_first[i];
			delete array_first;
		}

		//формирование вычисляемой матрицы в каждом процессе
		array_second = new int* [size_r_tmp];
		array_first = new int* [size_r_tmp + 2];//к блоку данных добавляются 2 строки для записи граничных данных из других процессов
		for (int i = 0; i < size_r_tmp + 2; ++i) {
			array_first[i] = new int[size];
			if (i < size_r_tmp) array_second[i] = new int[size];
			if (i > 0 && i < (size_r_tmp + 1))
				for (int j = 0; j < size; ++j) {
					array_first[i][j] = array_second[i - 1][j] = r_tmp[(i - 1) * size + j];
				}
		}
		delete r_tmp;

		//вычисление поледней необрабатываемой строки
		//необходимо, чтобы на место этой строки записывать граничные данные из таблицы другого процесса
		int q = size_r_tmp - 1;
		while (q > -1 && array_second[q][0] == -1) --q;
		if (q == -1) q = 0;

		while (1) {
			req_r = new MPI_Request;
			req_l = new MPI_Request;
			status = new MPI_Status;

			//отправка граничной верхней строки предыдущему процессу
			MPI_Isend(array_first[1], size, MPI_INT, (myrank - 1 + workers) % workers, 0, comm, req_l);
			//отправка граничной нижней строки следующему процессу
			MPI_Isend(array_first[q + 1], size, MPI_INT, (myrank + 1) % workers, 1, comm, req_r);

			//приём граничной верхней строки следующего процесса
			MPI_Irecv(array_first[q + 2], size, MPI_INT, (myrank + 1) % workers, 0, comm, req_r);
			//приём граничной нижней строки предыдущего процесса
			MPI_Irecv(array_first[0], size, MPI_INT, (myrank - 1 + workers) % workers, 1, comm, req_l);

			//ожидание всех строк
			MPI_Wait(req_r, status);
			MPI_Wait(req_l, status);
			delete req_r;
			delete req_l;
			delete status;

			//блок вычислений
			noChanges = true;
			for (int i = 0; i < size_r_tmp; i++) {
				if (array_first[(i + 1) % size_r_tmp + 1][0] == -1) break;
				for (int j = 0; j < size; j++) {
					neighbors += array_first[i + 1][(j + 1) % size];//правый
					neighbors += array_first[i + 2][j];//нижний
					neighbors += array_first[i + 1][(j - 1 + size) % size];//левый
					neighbors += array_first[i][j];//верхний
					/*S(t) – исходное состояние; N(t) – количество «живых» соседей.
					Описание переходов: S(t) = 1, N(t) < 2 – > S(t + 1) = 0;
										N(t) < 4 – > S(t + 1) = 1;
										N(t) = 4 – > S(t + 1) = 0;
										S(t) = 0, N(t) = 2 – > S(t + 1) = 1.*/
					array_second[i][j] = (neighbors == 2) + array_first[i + 1][j] * (neighbors == 3);
					neighbors = 0;
					noChanges = (array_second[i][j] == array_first[i + 1][j]) && noChanges;
				}
			}

			r_tmp = new int[workers];//буфер в мастер-процессе для сбора переменной noChanches
			MPI_Gather(&noChanges, 1, MPI_INT, r_tmp, 1, MPI_INT, 0, comm);//сбор переменной noChanches со всех процессов в мастер-процесс

			--iterations;
			//вычисление общей переменной наличия изменений
			if (!myrank) {
				noChanges = 0;
				for (int i = 0; i < workers; ++i)
					noChanges += r_tmp[i];
				noChanges = noChanges != workers && iterations;//если хоть один процесс возвращает 0 и итерации не закончились, то работа продолжается
			}
			delete r_tmp;
			//рассылка всем процессам принятого решения о продолжении работы
			MPI_Bcast(&noChanges, 1, MPI_INT, 0, comm);

			if (noChanges) {//продолжать работу
				for (int i = 0; i < size_r_tmp; ++i) {
					int* c = array_first[i + 1];
					array_first[i + 1] = array_second[i];
					array_second[i] = c;
				}
			}
			else {//завершаем работу
				for (int i = 0; i < size_r_tmp + 2; ++i) delete array_first[i];
				delete array_first;

				r_tmp = new int[size_r_tmp * size];//буфер для пересылки блока матрицы
				tmp = new int[workers * size_r_tmp * size];//буфер для сбора всей матрицы в мастер-процессе

				for (int i = 0; i < size_r_tmp; ++i)
					for (int j = 0; j < size; ++j)
						r_tmp[i * size + j] = array_second[i][j];

				//сбор блоков матрицы в буфер мастера-процесса
				MPI_Gather(r_tmp, size_r_tmp * size, MPI_INT, tmp, size_r_tmp * size, MPI_INT, 0, comm);

				delete r_tmp;
				for (int i = 0; i < size_r_tmp; ++i) delete array_second[i];
				delete array_second;

				//формирование матрицы из буфера
				if (!myrank) {
					array_second = new int* [size];
					for (int i = 0; i < size; ++i) {
						array_second[i] = new int[size];
						for (int j = 0; j < size; ++j) array_second[i][j] = tmp[i * size + j];
					}
					delete tmp;
					delay = clock() - delay;
					cout << "Calculation time = " << delay << endl;
				}
				break;
			}
		}
	}
}

int main(int argc, char** argv) {
	int myrank, ranksize;
	MPI_Init(&argc, &argv);//Инициализация MPI
	//Определяем свой номер в группе:
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	//Определяем размер группы:
	MPI_Comm_size(MPI_COMM_WORLD, &ranksize);
	int iter, mode, size;
	int** array_first = NULL, ** array_second = NULL;
	int* buf = new int[3];

	if (!myrank) {
		int seed;
		cout << "Enter seed of random: ";
		cin >> seed;
		srand(seed);
		cout << "Enter number of iterations: ";
		cin >> iter;
		cout << "Enter size of array: ";
		cin >> size;
		cout << "Enter number of mode (0. With output / 1. Without output): ";
		cin >> mode;

		buf[0] = iter; buf[1] = size; buf[2] = mode;
		array_first = new int* [size];
		array_second = new int* [size];

		for (int i = 0; i < size; i++) {
			array_first[i] = new int[size];
			array_second[i] = new int[size];
			for (int j = 0; j < size; j++) {
				array_first[i][j] = rand() % 2;
			}
		}

		if (!mode) WithOutput(array_first, array_second, size, iter);
	}

	//рассылка всем процессам информации об количестве итераций, размере матрицы, режиме работы
	MPI_Bcast(buf, 3, MPI_INT, 0, MPI_COMM_WORLD);
	iter = buf[0]; size = buf[1]; mode = buf[2];
	delete[] buf;

	if (mode) {
		WithoutOutputParal(array_first, array_second, size, iter, myrank, ranksize);
	}

	MPI_Finalize();
	return 0;
}