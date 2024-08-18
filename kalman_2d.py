
import matplotlib.pyplot as plt
import numpy as np
import csv

def read_csv(file_path):
    data = []
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            try:
                data.append([float(value) for value in row])
            except ValueError:
                print("Error parsing CSV data")
    return np.array(data)

def kalman_filter(z, x_est, P, Q, R):
    x_pred = x_est  # State prediction
    P_pred = P + Q  # Covariance prediction

    # Measurement update
    K = P_pred / (P_pred + R)  # Kalman gain
    x_est = x_pred + K * (z - x_pred)  # State estimate
    P = (1 - K) * P_pred  # Covariance estimate

    return x_est, P

file_path = "D:\Code\mpu_blah_data.csv"
x = read_csv(file_path)

Q = 1e-3  # Process noise covariance
R = 1e-2  # Measurement noise covariance
x_est = 0  # Initial estimate
P = 1  # Initial covariance estimate

filtered = []

for z in x:
    x_est, P = kalman_filter(z, x_est, P, Q, R)
    filtered.append(x_est)


plt.plot(x,color='r')
plt.plot(filtered,color='b')
#plt.xlim(0, 60)
plt.show()