import pandas as pd
import argparse
import matplotlib.pyplot as plt

# Step 1: Set up command-line argument parsing
def parse_args():
    parser = argparse.ArgumentParser(description="Calculate mean depth of coverage in 1kb bins.")
    parser.add_argument('-i', '--input', required=True, help="Path to input depth file (tab-separated)")
    parser.add_argument('-o', '--output', required=True, help="Path to output BED file with 1kb bin coverage")
    parser.add_argument('-p', '--plot', required=False, help="Path to save the plot image (e.g., plot.png)")
    return parser.parse_args()

# Step 2: Main function to compute the coverage
def compute_coverage(input_file, output_file):
    # Load depth file into a DataFrame
    depth_data = pd.read_csv(input_file, sep='\t', header=None, names=['chr', 'pos', 'depth'])

    # Create 1kb bins
    depth_data['bin'] = (depth_data['pos'] - 1) // 1000  # Divide by 1000 and floor to get the bin number

        # Calculate the mean depth for each bin
    bin_means = depth_data.groupby(['chr', 'bin'])['depth'].mean().reset_index()

    # Create the output in BED format with 1kb windows
    bin_means['start'] = bin_means['bin'] * 1000 + 1
    bin_means['end'] = (bin_means['bin'] + 1) * 1000

    # Save the results to a new file
    bin_means[['chr', 'start', 'end', 'depth']].to_csv(output_file, sep='\t', header=False, index=False)

    print(f"Coverage in 1kb bins has been saved to '{output_file}'.")

# Step 3: Plot the mean depth per 1kb bin and save it to a file
def plot_coverage(input_file, plot_file=None):
    # Load depth file into a DataFrame
    depth_data = pd.read_csv(input_file, sep='\t', header=None, names=['chr', 'pos', 'depth'])

    # Create 1kb bins
    depth_data['bin'] = (depth_data['pos'] - 1) // 1000

    # Calculate the mean depth for each bin
    bin_means = depth_data.groupby(['chr', 'bin'])['depth'].mean().reset_index()

# Plot the mean depth per 1kb bin
    plt.figure(figsize=(10, 6))
    plt.plot(bin_means['bin'] * 1000 + 1, bin_means['depth'], marker='o', linestyle='-', color='b')
    plt.xlabel('Position (start of 1kb bin)')
    plt.ylabel('Mean Depth')
    plt.title('Mean Depth of Coverage in 1kb Bins')
    plt.tight_layout()

    # Save the plot to the specified file
    if plot_file:
        plt.savefig(plot_file)
        print(f"Plot saved to '{plot_file}'.")
    else:
        # If no file is specified, show the plot interactively
        plt.show()


# Step 4: Main function to handle command-line arguments
def main():
    args = parse_args()  # Parse command-line arguments
    compute_coverage(args.input, args.output)  # Compute and save coverage
    plot_coverage(args.input, args.plot)  # Plot and save to file if specified

if __name__ == '__main__':
    main()
    
