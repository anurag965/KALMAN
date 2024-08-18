#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include<ctime>


struct Vector3d {
    double v[3];

    Vector3d() : v{0.0, 0.0, 0.0} {}

    Vector3d(double x, double y, double z) : v{x, y, z} {}

    double& operator[](int i) { return v[i]; }
    const double& operator[](int i) const { return v[i]; }

    // Vector addition
    Vector3d operator+(const Vector3d& other) const {
        return Vector3d(v[0] + other.v[0], v[1] + other.v[1], v[2] + other.v[2]);
    }

    // Vector subtraction
    Vector3d operator-(const Vector3d& other) const {
        return Vector3d(v[0] - other.v[0], v[1] - other.v[1], v[2] - other.v[2]);
    }

    // Scalar multiplication
    Vector3d operator*(double scalar) const {
        return Vector3d(v[0] * scalar, v[1] * scalar, v[2] * scalar);
    }

    // Dot product
    double dot(const Vector3d& other) const {
        return v[0] * other.v[0] + v[1] * other.v[1] + v[2] * other.v[2];
    }
};

struct Matrix3x3 {
    double m[3][3];

    Matrix3x3() {
        // Initialize to zero
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                m[i][j] = 0.0;
    }

    Matrix3x3(double a00, double a01, double a02,
              double a10, double a11, double a12,
              double a20, double a21, double a22) {
        m[0][0] = a00; m[0][1] = a01; m[0][2] = a02;
        m[1][0] = a10; m[1][1] = a11; m[1][2] = a12;
        m[2][0] = a20; m[2][1] = a21; m[2][2] = a22;
    }

    // Identity matrix
    static Matrix3x3 identity() {
        return Matrix3x3(1, 0, 0, 0, 1, 0, 0, 0, 1);
    }

    // Matrix multiplication
    Matrix3x3 operator*(const Matrix3x3& other) const {
        Matrix3x3 result;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                result.m[i][j] = m[i][0] * other.m[0][j] +
                                 m[i][1] * other.m[1][j] +
                                 m[i][2] * other.m[2][j];
        return result;
    }

    // Matrix-vector multiplication
    Vector3d operator*(const Vector3d& v) const {
        return Vector3d{
            m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
        };
    }

    // Matrix transpose
    Matrix3x3 transpose() const {
        return Matrix3x3(
            m[0][0], m[1][0], m[2][0],
            m[0][1], m[1][1], m[2][1],
            m[0][2], m[1][2], m[2][2]
        );
    }

    // Matrix inversion using Gaussian elimination
    Matrix3x3 inverse() const {
        Matrix3x3 result = Matrix3x3::identity();
        Matrix3x3 m = *this;
        
        for (int i = 0; i < 3; ++i) {
            double maxVal = fabs(m.m[i][i]);
            int maxRow = i;
            for (int k = i + 1; k < 3; ++k) {
                if (fabs(m.m[k][i]) > maxVal) {
                    maxVal = fabs(m.m[k][i]);
                    maxRow = k;
                }
            }

            if (fabs(maxVal) < 1e-10) {
                throw std::runtime_error("Matrix is singular and cannot be inverted");
            }

            if (maxRow != i) {
                std::swap(m.m[i], m.m[maxRow]);
                std::swap(result.m[i], result.m[maxRow]);
            }

            double pivot = m.m[i][i];
            for (int j = 0; j < 3; ++j) {
                m.m[i][j] /= pivot;
                result.m[i][j] /= pivot;
            }

            for (int k = 0; k < 3; ++k) {
                if (k != i) {
                    double factor = m.m[k][i];
                    for (int j = 0; j < 3; ++j) {
                        m.m[k][j] -= factor * m.m[i][j];
                        result.m[k][j] -= factor * result.m[i][j];
                    }
                }
            }
        }
        return result;
    }

    // Scalar multiplication
    Matrix3x3 operator*(double scalar) const {
        Matrix3x3 result;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                result.m[i][j] = m[i][j] * scalar;
        return result;
    }

    // Matrix subtraction
    Matrix3x3 operator-(const Matrix3x3& other) const {
        Matrix3x3 result;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                result.m[i][j] = m[i][j] - other.m[i][j];
        return result;
    }

    // Matrix addition
    Matrix3x3 operator+(const Matrix3x3& other) const {
        Matrix3x3 result;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                result.m[i][j] = m[i][j] + other.m[i][j];
        return result;
    }
};

// Kalman Filter class
class KalmanFilter {
public:
    KalmanFilter(const Matrix3x3& F, const Matrix3x3& H, const Matrix3x3& Q, const Matrix3x3& R)
        : F(F), H(H), Q(Q), R(R) {
        P = Matrix3x3::identity();
        x = Vector3d();
    }

    void predict() {
        x = F * x;
        P = F * P * F.transpose() + Q;
    }

    void update(const Vector3d& z) {
        Matrix3x3 S = H * P * H.transpose() + R;
        Matrix3x3 K = P * H.transpose() * S.inverse();
        x = x + K * (z - H * x);
        P = (Matrix3x3::identity() - K * H) * P;
    }

    Vector3d getState() const { return x; }

private:
    Matrix3x3 F; // State transition matrix
    Matrix3x3 H; // Measurement matrix
    Matrix3x3 Q; // Process noise covariance
    Matrix3x3 R; // Measurement noise covariance
    Matrix3x3 P; // Estimate covariance
    Vector3d x;  // State estimate
};

// Function to read CSV file into a vector of Vector3d
std::vector<Vector3d> readCSV(const std::string& filename) {
    std::vector<Vector3d> data;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        double x, y, z;

        std::getline(ss, item, ',');
        x = std::stod(item);
        std::getline(ss, item, ',');
        y = std::stod(item);
        std::getline(ss, item, ',');
        z = std::stod(item);

        data.emplace_back(x, y, z);
    }

    file.close();
    return data;
}

// Function to write a vector of Vector3d to a CSV file

void writeCSV(const std::string& filename, const std::vector<Vector3d>& data) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    for (const auto& v : data) {
        file << v[0] << "," << v[1] << "," << v[2] << "\n";
    }

    file.close();
}

int main() {
    
    try {
        
        // Example CSV file paths
        std::string inputFile = "D:\\Code\\data.csv";
        std::string outputFile = "D:\\Code\\output_data.csv";

        // Read data from CSV file
        std::vector<Vector3d> data = readCSV(inputFile);
        
        // Kalman filter setup
        Matrix3x3 F = Matrix3x3::identity(); // State transition matrix
        Matrix3x3 H = Matrix3x3::identity(); // Measurement matrix
        Matrix3x3 Q = Matrix3x3(0.005, 0, 0, 0, 0.005, 0, 0, 0, 0.005); // Process noise covariance
        Matrix3x3 R = Matrix3x3(0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05); // Measurement noise covariance

        KalmanFilter kf(F, H, Q, R);

        std::vector<Vector3d> filteredData;

        
        for (const auto& v : data) {
            kf.predict();
            kf.update(v);
            filteredData.push_back(kf.getState());
        }
        
        // Write processed data to a new CSV file
        writeCSV(outputFile, filteredData);
       
        std::cout << "Data processing complete. Results saved to " << outputFile << std::endl;
    } catch (const std::runtime_error& e){
        std::cerr << "Error: " << e.what() << std::endl;
    }
    
    
    

    return 0;
}
