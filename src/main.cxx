#define _USE_MATH_DEFINES
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

double pseudoVoigt(double x, double h, double x0, double sigma, double gamma, double eta); /* x is a vector, needs to be rewritten */
double scherrerEquation(double fwhm, double wavelength = 1.5406, double theta = 0.0, double shapeFactor = 0.9);
std::pair<std::vector<double>, std::vector<double>> readXrdData(const std::string& filePath);
std::vector<double> polyfit(const std::vector<double>& x, const std::vector<double>& y, double degree);
std::vector<double> evaluatePolynomial(const std::vector<double>& x, const std::vector<double>& coeffs);
std::vector<double> baselineCorrection(const std::vector<double>& x, const std::vector<double>& y, double degree);
std::pair<std::vector<double>, std::vector<double>> detectPeaks(const std::vector<double>& x, const std::vector<double>& y);

int main() {
	std::string filePath = "./CU_B_OBS.txt";
	std::vector<double> [x, y] = readXrdData(filePath);
	std::vector<double> yCorrected = baselineCorrection(x, y, 2);
	std::vector<double> [peaks, ySmooth] = detectPeaks(x, yCorrected);

	if (peaks.size() == 0) {
		std::cout << "No peaks detected." << std::endl;
		return 1;
	}

	std::vector<double> peakParams = fitPeaks(x, yCorrected, peaks);

	/* Code to plot the graphs here */

	for ()

	return 0;
}

double pseudoVoigt(double x, double h, double x0, double sigma, double gamma, double eta) {
	double gaussian = std::exp(-((x - x0) * (x - x0)) / (2 * sigma * sigma));
	double lorentzian = ((gamma / 2) * (gamma / 2)) / ((x - x0) * (x - x0) + (gamma / 2) * (gamma / 2));
	return h * (eta * gaussian + (1 - eta) * lorentzian);
}

double scherrerEquation(double fwhm, double wavelength=1.5406, double theta=0.0, double shapeFactor=0.9) {
	double thetaRad = theta * (M_PI / 180.0);
	return (shapeFactor * wavelength) / (fwhm * std::cos(thetaRad));
}

std::pair<std::vector<double>, std::vector<double>> readXrdData(const std::string& filePath) {
	std::ifstream file(filePath);
	std::vector<double> x, y;

	if (!file.is_open()) {
		std::cerr << "Error: Could not open file: " << filePath << std::endl;
		return {x, y};
	}

	std::string line;
	while(std::getline(file, line)) {
		std::istringstrean iss(line);
		double xVal, yVal;
		if (iss >> xVal >> yVal) {
			x.push_back(xVal);
			y.push_back(yVal);
		}
	}

	file.close();
	return {x, y};
}

std::vector<double> polyfit(const std::vector<double>& x, const std::vector<double>& y, double degree) {
	int n = x.size();
	Eigen::MatrixXd X(n, degree + 1);
	Eigen::VectorXd Y(n);

	// Construct Vandermonde matrix X and vector Y;
	for (int i = 0; i < n; ++i) {
		Y(i) = y[i];
		double val = 1.0;
		for (int j = 0; j <= degree; ++j) {
			X(i, j) = val;
			val *= x[i];
		}
	}

	// Solve for coeffs using normal equations: (XᵀX)c = XᵀY
	Eigen::VectorXd coeffs = (X.transpose() * X).ldlt().solve(X.transpose() * Y);

	return std::vector<double>(coeffs.data(), coeffs.data() + coeffs.size());
}

std::vector<double> evaluatePolynomial(const std::vector<double>& x, const std::vector<double>& coeffs) {
	std::vector<double> result(x.size(), 0.0);
	for (size_t i = 0; i < x.size(); ++i) {
		double val = 0.0;
		double xi = 1.0;
		for (size_t j = 0; j < coeffs.size(); ++j) {
			val += coeffs[j] * xi;
			xi *= x[i];
		}
		result[i] = val;
	}
	return result;
}

std::vector<double> baselineCorrection(const std::vector<double>& x, const std::vector<double>& y, double degree) {
	// 1. Fit polynomial to data
	std::vector<double> coeffs = polyfit(x, y, degree);

	// 2. Evaluate polynomial(baseline) at the given x values
	std::vector<double> baseline = evaluatePolynomial(x, coeffs);

	// 3. Subtract baseline from original y values to get corrected data
	std::vector<double> correctedData(y.size(), 0.0);
	for (size_t i = 0; i < y.size(); ++i) {
		correctedData[i] = y[i] - baseline[i];
	}

	return correctedData; /* This + last 2 functions were just 4 lines in python (T_T) */
}

std::pair<std::vector<double>, std::vector<double>> detectPeaks(const std::vector<double>& x, const std::vector<double>& y) {

}
