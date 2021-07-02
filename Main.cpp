#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <omp.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <sysinfoapi.h>
#include <windows.h>

using namespace std;

//Note that this is a serial implementation with a periodic grid
vector<vector<bool>> grid, new_grid;
int imax, jmax;
int max_steps;
int todo;
string files[12]; // CHANGE VALUE TO YOUR NUMBER OF CORES!!!
int freq;
int ppm;
int nthreads = omp_get_max_threads();

int step;

int num_neighbours(int ii, int jj)
{
	int ix, jx;
	int cnt = 0;
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
			if (i != 0 || j != 0)
			{
				ix = (i + ii + imax) % imax;
				jx = (j + jj + jmax) % jmax;
				if (grid[ix][jx]) cnt++;
			}
	return cnt;
}

void grid_to_file(int it)
{
	stringstream fname;
	fstream f1;
	fname << "output" << "_" << it << ".dat";
	f1.open(fname.str().c_str(), ios_base::out);
	for (int i = 0; i < imax; i++)
	{
		for (int j = 0; j < jmax; j++)
			f1 << grid[i][j] << "\t";
		f1 << endl;
	}
	f1.close();
}

void grid_to_pbm(int it)
{	
	//How many threads can we run at the same time?

	//Creating a char array as a buffer: writing large files to buffer is much faster than writing them char by char!
	//I have opted to create a vector instead of an array because of easy resizing.
	vector<char> lines;

	//I know the number of bytes I will need to allocate for my buffer. 2 for each pixel ('0' or '1' and '\t') and 1 for each newline
	lines.reserve(imax*jmax*2+imax);

	//We start the parallel region
	#pragma omp parallel
	{
		//We start a single thread to create all tasks
		#pragma omp single
		{
			//Looping through the number of threads we can create
			for (int num = 0; num < nthreads; num++)
			{
				//We create a task for each possible thread
				#pragma omp task
				{
					//Each task loops over an i/nthread sized chunk of lines
					for (int i = num*imax/nthreads; i < (num+1)*imax/nthreads; i++)
					{
						for (int j = 0; j < jmax; j++)
						{
							//Apologies for the complicated indices.
							//Adding a 1 or a 0 depending on grid[i][j]
                            lines[((i*(jmax))+j)*2+i] = grid[i][j] ? '1' : '0';
							//Adding a '\t'
							lines[((i*(jmax))+j)*2+i+1] = '\t';
						}
						//Adding a newline
						lines[(i+1)*jmax*2+i] = '\n';
					}
				}
			}
		}
	}

	string s(&lines[0], &lines[0]+imax*jmax*2+imax);

	files[(it/freq)%nthreads] = s;


	if ((it+1)%nthreads==0 || it == max_steps-1)
	{
		#pragma omp parallel
		{
			//We start a single thread to create all tasks
			#pragma omp single
			{
				for(int i = 0; i<min(nthreads,todo);i++)
				{
					#pragma omp task
					{			
						//Creating the filestream
						ofstream f1;

						//Name of the file
						string fname;
						fname += "output_" + to_string((it/nthreads)*nthreads+i*freq) + ".pbm";
						//Opening the file
						f1.open(fname, ios_base::out);

						//Creating the file's header
						string header;
						header += "P1 " + to_string(imax) + " " + to_string(jmax) + "\n";
						f1 << header;
						f1 << files[i];
					}
				}
			}
		}
		todo-=nthreads;
	}
	//And finally we write our buffer! Unfortunately this operation cannot be parallelised (only one thread can access the disk at a time)
}


void grid_to_ppm(int it, int pixels)
{
	string fname;
	fname = fname + "images/output_" + to_string(it) + ".ppm";
	ofstream ofs(fname, ios_base::out | ios_base::binary);
    
	cout << fname << endl;
	{
		ofs << "P6" << endl << imax*(pixels+1) << ' ' << jmax*(pixels+1) << endl << "255" << endl;
		for (auto i = 0; i < imax*(pixels+1); i++)
		{
			for (auto j = 0; j < jmax*(pixels+1); j++)
			if (i%(pixels+1) == 0 || j%(pixels+1) == 0)
			{
				ofs << (char) (192) << (char) (192) << (char) (192);
			}
			else
			{
				ofs << (char) ((!grid[i/(pixels+1)][j/(pixels+1)])*255) << (char) (!grid[i/(pixels+1)][j/(pixels+1)]*255) << (char) (!grid[i/(pixels+1)][j/(pixels+1)]*255);       // red, green, blue
			}
		}
	}
    ofs.close();
}

void do_iteration(void)
{
	//We start the parallel region
	#pragma omp parallel
	{
		//We start a single thread to create all tasks
		#pragma omp single
		{
			//Looping through the number of threads we can create
			for (int num = 0; num < nthreads; num++)
			{
				//We create a task for each possible thread
				#pragma omp task
				{
					//Each task loops over an i/nthread sized chunk of lines
					for (int i = num*imax/nthreads; i < (num+1)*imax/nthreads; i++)
					{
						int num_n;
						for (int j = 0; j < jmax; j++)
						{
                            new_grid[i][j] = grid[i][j];
							num_n = num_neighbours(i, j);
							
							if (grid[i][j])
							{
								if (num_n != 2 && num_n != 3)
									new_grid[i][j] = false;
							}
							else if (num_n == 3) new_grid[i][j] = true;
						}
					}
				}
			}
		}
	}
	grid.swap(new_grid);
}

int main(int argc, char *argv[])
{

	MEMORYSTATUSEX statex;

	statex.dwLength = sizeof (statex);

	GlobalMemoryStatusEx (&statex);

	//cout << "TOTAL RAM: " << statex.ullTotalPhys/1024/1024 << " AVAILABLE RAM: "<<statex.ullAvailPhys/1024/1024 << endl;

	int TOTMEM = statex.ullAvailPhys/1024/1024;
	int step = 0;

	string command;
	cout << "Welcome to our little parallel code for Conway's Game of Life!" <<endl;

	cout << "With what frequency would you like us to print updates?" << endl;
	cin >> command;
	freq = stoi(command);

	ppm = 2;

	while(ppm !=0 && ppm!=1)
	{
		cout << "Would you like us to print pbm files or ppm files? Input 0 for pbm files or 1 for ppm files" << endl;
		cin >> command;
		ppm = stoi(command);
	}

	cout << "Input imax: ";
	cin >> command;
	imax = stoi(command);

	cout << "\nInput jmax: ";
	cin >> command;
	jmax = stoi(command);

	cout << "\nInput iter: ";
	cin >> command;
	max_steps = stoi(command);
	todo = stoi(command);


	srand(time(NULL));
	grid.resize(imax, vector<bool>(jmax));
	new_grid.resize(imax, vector<bool>(jmax));
	
	//How many threads can we run at the same time?
	nthreads = omp_get_max_threads();

	#pragma for parallel
	for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++) grid[i][j] = (rand() % 2);

	double start = omp_get_wtime();
	for (int n = 0; n < max_steps; n++)
	{
		//TIMING IS VERY EXPENSIVE!
		//200,200,400 is 5.7 secs without timing and 8.5 with!

		do_iteration();
		
		if(n%freq ==0)
		{
			if (ppm==1)
			{
				grid_to_ppm(n,10);
			}
			else
			{
				grid_to_pbm(n);
			}
		}		
	};
	double end = omp_get_wtime();
	cout << "Finding step " << max_steps << " for a grid of size " <<imax << "x" << jmax << " took " << end-start <<" seconds \n";
	return 0;
}
