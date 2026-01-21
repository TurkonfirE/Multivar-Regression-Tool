# Multivariate Regression Tool

A C++ tool that automatically performs multivariate least squares regression analysis on multiple models; testing linear, quadratic, exponential, and logarithmic models to find the best fit

## Features
- Matrix operations implemented from scratch (multiplication, transposition, inversion) in C++ to ensure efficiency when performed on big data
- Multivariate least squares regression
- R-squared calculations
- Automatic model selection based on best R-squared

## Usage
chmod +x modeler
./modeler -vf data.txt

## Input Format
Space-separated values, one row per line. Last column is the output variable.

## Example
1 1 10
2 1 13 
1 2 15
2 2 18

## Author
Samuel Jennings - University of Michigan
