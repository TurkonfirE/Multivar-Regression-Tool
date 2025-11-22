#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <string>

using namespace std;

class Model {
public:
    // Constructor
    Model(bool verbose, const std::string& filename);
    
    // Destructor
    ~Model();
    
    // Create model
    void make_model();
    
private:
    // Data information
    struct Matrix {
        vector<vector<double>> data;
        vector<string> colNames;
        vector<string> rowNames;
        vector<double> output;
    };

    Matrix dataMatrix;
    
    // Verbose
    bool verbose;

    // Best model + info
    string bestModel;
    double bestRSquared;
    Matrix bestCoefficients;
    // *** Model operations ***
    // Parse the data
    void parseData(const string& filename);

    // Run least squares and r-squared analyses
    void runAnalyses();

    // Display best model results + info
    void showResults();

    // Find the R Squared value
    double calculateRSquared(const Matrix& y_actual, const Matrix& y_predicted);

    // *** Matrix operations ***
    // Perform the least squares analysis
    Matrix leastSquares(const Matrix& input, const Matrix& output, bool exp);

    // Transform a matrix (lin, quad, exp, log)
    Matrix transform(const Matrix& m, string type);

    // Matrix multiplication
    Matrix multiply(const Matrix& a, const Matrix& b);

    // Transpose a matrix
    Matrix transpose(const Matrix& m);

    // Invert a matrix
    Matrix inverse(const Matrix& m);
};

#endif