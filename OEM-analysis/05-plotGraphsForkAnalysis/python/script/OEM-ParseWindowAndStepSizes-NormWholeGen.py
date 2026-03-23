import pandas as pd
import argparse

def calculate_oem(watson_file, crick_file, output_file, window_size, step_size):
    # Load bedgraph files
    watson_df = pd.read_csv(watson_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'value'])
    crick_df = pd.read_csv(crick_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'value'])

    # Remove NaN values in 'start' column
    watson_df.dropna(subset=['start'], inplace=True)
    crick_df.dropna(subset=['start'], inplace=True)

    # Ensure 'start' column contains integers
    watson_df['start'] = watson_df['start'].astype(int)
    crick_df['start'] = crick_df['start'].astype(int)

    # Calculate global sums for normalization
    global_watson_sum = watson_df['value'].sum()
    global_crick_sum = crick_df['value'].sum()

    results = []

    # Iterate through each unique chromosome
    for chrom in watson_df['chrom'].unique():
        watson_chrom_df = watson_df[watson_df['chrom'] == chrom]
        crick_chrom_df = crick_df[crick_df['chrom'] == chrom]

        # Define half of the window size for easy calculation of bounds
        half_window = window_size // 2

        # Iterate through each position in the chromosome using the step size
        for position in range(watson_chrom_df['start'].min(), watson_chrom_df['start'].max(), step_size):
            left_bound = position - half_window
            right_bound = position + half_window

            # Extract data for the window (left and right bounds)
            watson_left = watson_chrom_df[(watson_chrom_df['start'] >= left_bound) & (watson_chrom_df['start'] < position)]
            watson_right = watson_chrom_df[(watson_chrom_df['start'] >= position) & (watson_chrom_df['start'] < right_bound)]
            crick_left = crick_chrom_df[(crick_chrom_df['start'] >= left_bound) & (crick_chrom_df['start'] < position)]
            crick_right = crick_chrom_df[(crick_chrom_df['start'] >= position) & (crick_chrom_df['start'] < right_bound)]

            # Calculate values for each quadrant
            WL = watson_left['value'].sum() / global_watson_sum
            WR = watson_right['value'].sum() / global_watson_sum
            CL = crick_left['value'].sum() / global_crick_sum
            CR = crick_right['value'].sum() / global_crick_sum

            # Normalize Watson values
            WLn = WL / (WL + CL) if (WL + CL) != 0 else 0
            WRn = WR / (WR + CR) if (WR + CR) != 0 else 0

            # Calculate OEM
            OEM = WLn - WRn

            # Append result
            results.append([chrom, position, position + 1, OEM])

    # Convert results to DataFrame and save to file
    result_df = pd.DataFrame(results, columns=['chrom', 'start', 'end', 'OEM'])
    result_df.to_csv(output_file, sep='\t', index=False, header=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate OEM from Watson and Crick bedgraph files.")
    parser.add_argument('watson_file', help='Bedgraph file for the Watson strand.')
    parser.add_argument('crick_file', help='Bedgraph file for the Crick strand.')
    parser.add_argument('output_file', help='Output file for the OEM values.')
    parser.add_argument('-w', '--window_size', type=int, default=20000, help='Size of the sliding window in base pairs.')
    parser.add_argument('-s', '--step_size', type=int, default=1, help='Step size for the sliding window in base pairs.')
    args = parser.parse_args()

    calculate_oem(args.watson_file, args.crick_file, args.output_file, args.window_size, args.step_size)
