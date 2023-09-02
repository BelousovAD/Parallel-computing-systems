#include <iostream>
#include <fstream>
#include "omp.h"

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

void WithoutOutput(int** array_first, int** array_second, int& size, int& iterations) {
	int neighbors = 0;
	int** tmp = NULL;
	double start, stop;
	bool noChanges;
	start = omp_get_wtime();
	do {
		noChanges = true;
		#pragma omp parallel shared(array_first, array_second, size) firstprivate(neighbors)
		{
			#pragma omp for reduction(&&:noChanges) nowait 
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
					noChanges = (array_second[i][j] == array_first[i][j]) && noChanges;
				}
			}
		}
		tmp = array_first;
		array_first = array_second;
		array_second = tmp;
		iterations--;
	} while (iterations != 0 && !noChanges);
	stop = omp_get_wtime();
	cout << "Calculation time = " << stop - start << endl
		<< "Size = " << size << " X " << size << endl;
}

int main() {
	int iter, mode, seed, size;
	cout << "Enter seed of random: ";
	cin >> seed;
	srand(seed);
	cout << "Enter number of iterations: ";
	cin >> iter;
	cout << "Enter size of array: ";
	cin >> size;
	cout << "Enter number of mode (0. With output / 1. Without output): ";
	cin >> mode;
	int** array_first = new int* [size];
	int** array_second = new int* [size];
	for (int i = 0; i < size; i++) {
		array_first[i] = new int[size];
		array_second[i] = new int[size];
		for (int j = 0; j < size; j++) {
			array_first[i][j] = rand() % 2;
		}
	}
	if (mode) {
		WithoutOutput(array_first, array_second, size, iter);
	}
	else {
		WithOutput(array_first, array_second, size, iter);
	}
}