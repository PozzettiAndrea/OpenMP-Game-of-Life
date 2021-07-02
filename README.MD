# Parallisation for Conway's game of life

## Introduction
We aimed to parallelise a serial code implementation of ConwayGame of Life using OpenMP. This project involves attempted optimisations at GOL algorithm as well as the printing functions. 

## Structure
This programme includes four functions to implement Gonway's Game of Life. `num_neighbours(int ii, int jj)` is used to find every cell which in one of two possible states and the states of it eight neighnours. `grid_to_file(int it)`, the function which prints matrices output into different `.dat` files with its iteration number.`grid_to_pbm(int it)`, the main function to print data into pbm format at output. `do_iteration(void)` is used to loop over cells in matrices. 

## Installation

1) The software is installed by cloning the GitHub repository to your local machine. Fully functional program lives on Main branch.
2) For more advanced usage run the program using an IDE such as Visual Studio Community or a text editor such as Visual Studio Code with a compiler.

To compile make sure to turn on optimisation and openmp flags

g++ main.cpp -fopenmp -O3


## PLEASE NOTE:

UPDATE THE SIZE OF THE STRING ARRAY FILES[] to the number of cores of your computer or it might not work.

Unfortunately the frequency functionality is broken for pbm!

## Functionality

Running the project navigates the user through a user interface that allows them to choose what type of image to create from the GoL code. The user also could time the code if they wish with the parts of the code that are commented out. 

 - Grid tests were done for the parallel and serial case on a 12 core system, as well as for the parallel case on a 2 core system from a 100x100 grid to a 10000x10000 grid.

