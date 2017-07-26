#include "stdafx.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <omp.h>
#include <fstream>
#include <math.h>

using namespace std;

bool checkConfidenceInterval(int stepCount, double sum, double sqrSum);

int main()
{
	srand((double)time(0));
	int multi;

	ofstream outputFile;
	outputFile.open("serial.txt");

	cout << "...................... " << endl;

	outputFile << "\n" << "Parallel run optimised Results: " << endl;

	for (int i = 200; i < 2001; i = i + 200) {
		cout << "size: " << i << endl;
		int size = i;

		double sum = 0;
		double sqrSum = 0;

		// Assume doing 20 times is sufficient enough (proved later...)
		int stepCount = 20;
		for (int step = 0; step<stepCount; step += 1) {

            vector< vector<double> > a(size, vector<double>(size)), b(size, vector<double>(size)), c(size, vector<double>(size));

            // Initialize vectors to random numbers
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    a[i][j] = (double)rand() / RAND_MAX * 10;
                    b[i][j] = (double)rand() / RAND_MAX * 10;
                    c[i].push_back(0);
                }
            }

            double start = omp_get_wtime();

            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    for (int k = 0; k < size; k++)
                    {
                            c[i][j] += a[i][k] * b[k][j];
                    }
                }
            }

            double end = omp_get_wtime();
            double totalTime = end - start;

			sum += totalTime;
			sqrSum += totalTime*totalTime;
		}

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
