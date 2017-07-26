#include "stdafx.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include <fstream>
#include <math.h>

using namespace std;

double** create2DArray(int size);

void fillValuesRandomly(double **Array2D, int size);

double multiplySerial(int size, double **numArray1, double **numArray2, double **resultArray);

void free2DArray(double **Array2D, int size);

bool checkConfidenceInterval(int stepCount, double sum, double sqrSum);

int main()
{
	srand((double)time(0));
	int multi;

	ofstream outputFile;
	outputFile.open("serial.txt");

	cout << "...................... " << endl;

	outputFile << "\n" << "Serial Results: " << endl;

	for (int i = 200; i < 2001; i = i + 200) {
		cout << "size: " << i << endl;
		int size = i;

		double sum = 0;
		double sqrSum = 0;
		double seqTime = 0;

		// Initialize arrays
		double** initial2DArray1 = create2DArray(size);
		double** initial2DArray2 = create2DArray(size);
		double** resultArray = create2DArray(size);

		// Assume doing 20 times is sufficient enough (proved later...)
		int stepCount = 20;
		for (int step = 0; step<stepCount; step += 1) {

			fillValuesRandomly(initial2DArray1, size);
			fillValuesRandomly(initial2DArray2, size);

			double serialTime = multiplySerial(size, initial2DArray1, initial2DArray2, resultArray);

			sum += serialTime;
			sqrSum += serialTime*serialTime;

		}

		free2DArray(initial2DArray1, size);
		free2DArray(initial2DArray2, size);
		free2DArray(resultArray, size);

		bool isConfidentEnough = checkConfidenceInterval(stepCount, sum, sqrSum);
		double avgTime = sum / stepCount;

		cout << "Is within 95% confident: " << isConfidentEnough << endl;
		cout << "Average Time => " << avgTime << endl << endl;

		outputFile << size << " : " << avgTime << endl;
	}
	cout << "...................... " << endl;

	outputFile.close();

	cin >> multi;
	return 0;
}

double** create2DArray(int size)
{
	double** array2D;

	array2D = new double*[size];

	for (int h = 0; h < size; h++)
	{
		array2D[h] = new double[size];
	}

	return array2D;
}

void fillValuesRandomly(double **Array2D, int size)
{
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			double r = rand();
			Array2D[i][j] = r;
		}
	}
	return;
}

double multiplySerial(int size, double **numArray1, double **numArray2, double **resultArray)
{
	double sum;

	double start = omp_get_wtime();

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			double sum = 0;
			for (int k = 0; k < size; k++) {
				sum += numArray1[i][k] * numArray2[k][j];
			}
			resultArray[i][j] = sum;
		}
	}

	double end = omp_get_wtime();
	double totalTime = end - start;

	//cout << endl << "Time for Sequential Multiplication: " << totalTime << endl;

	return totalTime;
	//return resultArray;
}

bool checkConfidenceInterval(int stepCount, double sum, double sqrSum)
{
	double avg = sum / stepCount;
	double sqrAvg = sqrSum / stepCount;

	double stdDeviation = sqrt(sqrAvg - avg*avg);

	// For 95% confidence interval
	double expectedSamples = pow((100 * 1.96*stdDeviation) / (5 * avg), 2);

	if (expectedSamples < stepCount) {
		return true;
	}
	return false;
}

void free2DArray(double **Array2D, int size)
{
	for (int i = 0; i < size; i++) {
		delete[] Array2D[i];
	}
	//Free the array of pointers
	delete[] Array2D;
}
