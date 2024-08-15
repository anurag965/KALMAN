import serial
import csv
import time

# Set the serial port and baud rate
serial_port = 'COM4'  # Change this to your Arduino's COM port
baud_rate = 9600
file_name = 'mpu_blah_data.csv'

# Create a serial connection
ser = serial.Serial(serial_port, baud_rate)
time.sleep(2)  # Wait for the connection to be established

# Open a CSV file for writing
with open(file_name, mode='a', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)

    

    try:
        while True:
            # Read a line from the serial port
            line = ser.readline().decode('utf-8').rstrip()
            if line:
                # Split the line into a list
                data = line.split(',')
                # Write the data to the CSV file
                csv_writer.writerow(data)
                csv_file.flush() 
                #print(data)  # Optional: Print the data to console
    except KeyboardInterrupt:
        print("Logging stopped.")
    finally:
        ser.close()