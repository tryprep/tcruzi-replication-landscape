import pandas as pd
import argparse

def calculate_rfd(watson_file, crick_file, output_file, window_size, step_size):
    # Load bedgraph files
    watson_df = pd.read_csv(watson_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'value'])
    crick_df = pd.read_csv(crick_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'value'])

    # Remove rows with NaN values
    watson_df.dropna(subset=['start', 'end', 'value'], inplace=True)
    crick_df.dropna(subset=['start', 'end', 'value'], inplace=True)

    results = []

    # Iterate through each unique chromosome
    for chrom in watson_df['chrom'].unique():
        watson_chrom_df = watson_df[watson_df['chrom'] == chrom]
        crick_chrom_df = crick_df[crick_df['chrom'] == chrom]

        if watson_chrom_df.empty or crick_chrom_df.empty:
            continue  # Skip empty chromosomes

        # Get the min and max values for start positions
        start_min = watson_chrom_df['start'].min()
        start_max = watson_chrom_df['start'].max()

        if pd.isna(start_min) or pd.isna(start_max):
            continue  # Skip if start values are not valid

        # Iterate through each position in the chromosome using the step size
        for position in range(int(start_min), int(start_max) - window_size, step_size):
            left_bound = position
            right_bound = position + window_size

            # Extract data for the window
            watson_window = watson_chrom_df[(watson_chrom_df['start'] >= left_bound) & (watson_chrom_df['start'] < right_bound)]
            crick_window = crick_chrom_df[(crick_chrom_df['start'] >= left_bound) & (crick_chrom_df['start'] < right_bound)]

            # Sum values for Watson and Crick within the window
            W = watson_window['value'].sum()
            C = crick_window['value'].sum()

            # Calculate RFD
            RFD = (C - W) / (C + W) if (C + W) != 0 else 0

            results.append([chrom, position, right_bound, RFD])

    # Convert results to DataFrame and save
    result_df = pd.DataFrame(results, columns=['chrom', 'start', 'end', 'RFD'])
    result_df.to_csv(output_file, sep='\t', index=False, header=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate RFD from Watson and Crick bedgraph files.")
    parser.add_argument('watson_file', help='Bedgraph file for the Watson strand.')
    parser.add_argument('crick_file', help='Bedgraph file for the Crick strand.')
    parser.add_argument('output_file', help='Output file for the RFD values.')
    parser.add_argument('-w', '--window_size', type=int, default=20000, help='Size of the sliding window in base pairs.')
    parser.add_argument('-s', '--step_size', type=int, default=1, help='Step size for the sliding window in base pairs.')
    args = parser.parse_args()

    calculate_rfd(args.watson_file, args.crick_file, args.output_file, args.window_size, args.step_size)
