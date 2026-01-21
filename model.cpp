#include "model.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

Model::Model(bool verbose, const std::string& filename) 
    : verbose(verbose), bestRSquared(0.0) {
    // Load data in constructor
    parseData(filename);
}

Model::~Model() {
    // Destructor
}

void Model::make_model() {
    // Run model comparisons
    runAnalyses();
    
    // Display results
    showResults();
}

void Model::parseData(const std::string& filename) {
    ifstream file(filename);

    if (!file.is_open()) {
        throw runtime_error("Could not open file!");
    }

    string line;

    // Add to data line by line
    while (getline(file, line)) {
        vector<double> row; 
        stringstream ss(line);
        double val;

        while (ss >> val) {
            row.push_back(val);
        }
        dataMatrix.output.push_back(row[row.size() - 1]);
        row.pop_back();
    
        dataMatrix.data.push_back(row);
    }

    if (verbose) {
        cout << "Loaded " << dataMatrix.data.size() << " rows with " 
        << dataMatrix.data[0].size() << " columns" << endl;
    
        int rowsToShow = min(5, (int)dataMatrix.data.size());
        cout << "The first " << rowsToShow << " rows (inputs -> outputs):\n";
        for (int i = 0; i < rowsToShow; ++i) {
            for (int j = 0; j < dataMatrix.data[0].size(); ++j) {
                cout << dataMatrix.data[i][j] << " ";
            }
            cout << "-> " << dataMatrix.output[i] << endl;
        }
    }
}

void Model::runAnalyses() {
    // Run least squares on different models
    // Calculate R-squared for each
    // Keep track of best model
    // (X^T*X)^-1*X^T*y
    Matrix y;
    for (double val : dataMatrix.output) {
        vector<double> row;
        row.push_back(val);
        y.data.push_back(row);
    }

    Matrix coefficients_lin, coefficients_quad, coefficients_exp, coefficients_log;
    vector<double> rValues;

    Matrix lin = transform(dataMatrix, "linear");
    coefficients_lin = leastSquares(lin, y, false);
    Matrix y_predicted_lin = multiply(lin, coefficients_lin);
    double linRSquared = calculateRSquared(y, y_predicted_lin);
    rValues.push_back(linRSquared);
    if (verbose) {
        cout << "Linear R^2: " << linRSquared << endl;
    }
    
    if (y.data.size() >= 1 + 2 * dataMatrix.data[0].size()) {
        Matrix quad = transform(dataMatrix, "quadratic");
        coefficients_quad = leastSquares(quad, y, false);
        Matrix y_predicted_quad = multiply(quad, coefficients_quad);
        double quadRSquared = calculateRSquared(y, y_predicted_quad);
        rValues.push_back(quadRSquared);
        if (verbose) {
            cout << "Quadratic R^2: " << quadRSquared << endl;
        }
    }
    
    Matrix exp = transform(dataMatrix, "exponential");
    coefficients_exp = leastSquares(exp, y, true);
    Matrix y_predicted_exp = multiply(exp, coefficients_exp);
    double expRSquared = calculateRSquared(y, y_predicted_exp);
    rValues.push_back(expRSquared);
    if (verbose) {
        cout << "Exponential R^2: " << expRSquared << endl;
    }

    Matrix log = transform(dataMatrix, "logarithmic");
    coefficients_log = leastSquares(log, y, false);
    Matrix y_predicted_log = multiply(log, coefficients_log);
    double logRSquared = calculateRSquared(y, y_predicted_log);
    rValues.push_back(logRSquared);
    if (verbose) {
        cout << "Logarithmic R^2: " << logRSquared << endl;
    }

    double maxR = 0.0;
    int index = 0;
    for (int i = 0; i < rValues.size(); ++i) {
        if (abs(rValues[i]) > maxR) {
            maxR = rValues[i];
            index = i;
        }
    }
    bestRSquared = maxR;
    
    if (index == 0) {
        bestModel = "linear";
        bestCoefficients = coefficients_lin;
    }
    else if (index == 1) {
        bestModel = "quadratic";
        bestCoefficients = coefficients_quad;
    }
    else if (index == 2) {
        bestModel = "exponential";
        bestCoefficients = coefficients_exp;
    }
    else if (index == 3) {
        bestModel = "logarithmic";
        bestCoefficients = coefficients_log;
    }
}

// Equation specific output formatting
void Model::showResults() {
    cout << "\n---------------" << endl;
    cout << "  BEST MODEL" << endl;
    cout << "---------------" << endl;
    cout << "Model type: " << bestModel << endl;
    cout << "R-squared: " << bestRSquared << endl;
    cout << "Coefficients: ";
    int quadCount = 0;

    if (bestModel == "exponential") {
        cout << bestCoefficients.data[0][0];
        cout << " + e^(";
    }
    for (int i = 0; i < bestCoefficients.data.size(); ++i) {
        if (bestModel == "exponential" && i == 0) {
            // Nothing
        }
        else {
            cout << bestCoefficients.data[i][0];
        }
        if (bestModel == "linear") {
            if (i != 0) {
                cout << "*x" << i;
            }

            if (i != bestCoefficients.data.size() - 1) {
                cout << " + ";
            }
        }
        if (bestModel == "quadratic") {
            if (i == 0) {
                // Nothing
            }
            else if (i <= (bestCoefficients.data.size() - 1) / 2) {
                cout << "*x" << i;
            }
            else {
                ++quadCount;
                cout << "*x" << quadCount << "^2";
            }
            if (i != bestCoefficients.data.size() - 1) {
                cout << " + ";
            }
        }
        if (bestModel == "exponential") {
            if (i != 0) {
                cout << "*x" << i;
                if (i != bestCoefficients.data.size() - 1) {
                    cout << " + ";
                }
                else {
                    cout << ")";
                }
            }
        }
        if (bestModel == "logarithmic") {
            if (i != 0) {
                cout << "*ln(x" << i << ")";
            }
            
            if (i != bestCoefficients.data.size() - 1) {
                cout << " + ";
            }
        }
    }

    if (verbose) {
        cout << "\n\nThis model explains " << (bestRSquared * 100) 
        << "% of the variance in the output.\n";
    }
}

double Model::calculateRSquared(const Matrix& y_actual, const Matrix& y_predicted) {
    // r^2 = [(Ya * Yp) / (||Ya||*||Yp||)]^2.
    int n = y_actual.data.size(); // Num of data points
    // Find the mean of y
    double yMean = 0.0;
    double yPMean = 0.0;
    for (int i = 0; i < n; ++i) {
        yMean += y_actual.data[i][0];
        yPMean += y_predicted.data[i][0];
    }
    yMean /= n;
    yPMean /= n;

    double yDevMag = 0.0;
    double yPredDevMag = 0.0;
    double numerator = 0.0;
    
    for (int i = 0; i < n; ++i) {
        double yDev = y_actual.data[i][0] - yMean;
        double yPDev = y_predicted.data[i][0] - yPMean;

        yDevMag += yDev * yDev;
        yPredDevMag += yPDev * yPDev;
        numerator += yDev * yPDev;
    }

    yDevMag = sqrt(yDevMag);
    yPredDevMag = sqrt(yPredDevMag);

    double denominator = yDevMag * yPredDevMag;
    return (numerator / denominator) * (numerator / denominator);
}

Model::Matrix Model::leastSquares(const Matrix& input, const Matrix& output, bool exp) {
    // Perform matrix calculations
    auto inputTransposed = transpose(input);
    auto step1 = multiply(inputTransposed, input);
    auto step2 = inverse(step1);
    auto step3 = multiply(step2, inputTransposed);
    // Resulting coefficients (x*)

    if (exp) {
        // When doing least squares on exponential, you transform the output
        auto outputCopy = output;
        for (int i = 0; i < outputCopy.data.size(); ++i) {
            outputCopy.data[i][0] = log(outputCopy.data[i][0]);
        }
        return multiply(step3, outputCopy);
    }
    else {
        return multiply(step3, output); 
    } 
}

Model::Matrix Model::transform(const Matrix& m, string type) {
    Matrix result = m;
    
    // Nothing for linear model

    // Nothing for exponential (you take log of output)

    // Quadratic
    if (type == "quadratic") {
        for (int i = 0; i < m.data.size(); ++i) {
            for (int j = 0; j < m.data[0].size(); ++j) {
                result.data[i].push_back(m.data[i][j] * m.data[i][j]);
            }
        }
    }

    // Logarithmic
    else if (type == "logarithmic") {
        for (int i = 0; i < m.data.size(); ++i) {
            for (int j = 0; j < m.data[0].size(); ++j) {
                result.data[i][j] = log(m.data[i][j]);
            }
        }
    }
    
    // Add ones for every model!
    for (int i = 0; i < m.data.size(); ++i) {
            result.data[i].insert(result.data[i].begin(), 1.0);
    }

    return result;
}

Model::Matrix Model::multiply(const Matrix& a, const Matrix& b) {
    Matrix result;

    for (int i = 0; i < a.data.size(); ++i) {
        vector<double> row;

        for (int j = 0; j < b.data[0].size(); ++j) {
            double sum = 0;

            for (int k = 0; k < b.data.size(); ++k) {
                sum += a.data[i][k] * b.data[k][j];
            }

            row.push_back(sum);
        }

        result.data.push_back(row);
    }

    return result;
}

Model::Matrix Model::transpose(const Matrix& m) {
    Matrix result;

    for (int i = 0; i < m.data[0].size(); ++i) {
        vector<double> row;

        for (int j = 0; j < m.data.size(); ++j) {
            row.push_back(m.data[j][i]);
        }

        result.data.push_back(row);
    }
    return result;
}

Model::Matrix Model::inverse(const Matrix& m) {
    int n = m.data.size();
    
    // Check if square matrix
    if (n != m.data[0].size()) {
        throw std::runtime_error("Matrix must be square to invert");
    }
    
    // Create augmented matrix [A | I]
    vector<vector<double>> aug(n, vector<double>(2 * n));
    
    // Fill left side with m, right side with identity
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            aug[i][j] = m.data[i][j];
            aug[i][j + n] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    // Forward elimination
    for (int i = 0; i < n; ++i) {
        // Find pivot
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(aug[k][i]) > abs(aug[maxRow][i])) {
                maxRow = k;
            }
        }
        
        // Swap rows
        swap(aug[i], aug[maxRow]);
        
        // Check for singular matrix
        if (abs(aug[i][i]) < 1e-10) {
            throw std::runtime_error("Matrix is singular, cannot invert");
        }
        
        // Scale pivot row
        double pivot = aug[i][i];
        for (int j = 0; j < 2 * n; ++j) {
            aug[i][j] /= pivot;
        }
        
        // Eliminate column
        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = aug[k][i];
                for (int j = 0; j < 2 * n; ++j) {
                    aug[k][j] -= factor * aug[i][j];
                }
            }
        }
    }
    
    // Extract right half (the inverse)
    Matrix result;
    for (int i = 0; i < n; ++i) {
        vector<double> row;
        for (int j = n; j < 2 * n; ++j) {
            row.push_back(aug[i][j]);
        }
        result.data.push_back(row);
    }
    
    return result;
}

