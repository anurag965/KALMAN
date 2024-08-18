import pandas as pd
import matplotlib.pyplot as plt

def read_csv(filename):
    """Reads a CSV file and returns a DataFrame."""
    return pd.read_csv(filename, header=None, names=['ax', 'ay', 'az'])

def plot_data(df1, df2):
    """Plots data from two DataFrames on combined subplots."""
    fig, axs = plt.subplots(3, 1, figsize=(12, 15), sharex=True)

    # Plot 'ax' from both DataFrames in the first subplot
    axs[0].plot(df1.index, df1['ax'], label='File 1 - ax', color='r', linestyle='-', marker='')
    axs[0].plot(df2.index, df2['ax'], label='File 2 - ax', color='b', linestyle='-', marker='')
    axs[0].set_title('Comparison of ax')
    axs[0].set_ylabel('ax')
    axs[0].legend()
    axs[0].grid(True)

    # Plot 'ay' from both DataFrames in the second subplot
    axs[1].plot(df1.index, df1['ay'], label='File 1 - ay', color='r', linestyle='-', marker='')
    axs[1].plot(df2.index, df2['ay'], label='File 2 - ay', color='b', linestyle='-', marker='')
    axs[1].set_title('Comparison of ay')
    axs[1].set_ylabel('ay')
    axs[1].legend()
    axs[1].grid(True)

    # Plot 'az' from both DataFrames in the third subplot
    axs[2].plot(df1.index, df1['az'], label='File 1 - az', color='r', linestyle='-', marker='')
    axs[2].plot(df2.index, df2['az'], label='File 2 - az', color='b', linestyle='-', marker='')
    axs[2].set_title('Comparison of az')
    axs[2].set_xlabel('Index')
    axs[2].set_ylabel('az')
    axs[2].legend()
    axs[2].grid(True)

    # Adjust layout
    plt.tight_layout()
    plt.show()

def main():
    # File paths to the CSV files
    # file1 = 'D:\\Code\\mpu_blah_data.csv'
    file1 = 'D:\\Code\\data.csv'
    file2 = 'D:\\Code\\output_data.csv'

    # Read the CSV files
    df1 = read_csv(file1)
    df2 = read_csv(file2)

    # Plot the data
    plot_data(df1, df2)

if __name__ == '__main__':
    main()
