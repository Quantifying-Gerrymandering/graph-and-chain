import pickle
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import os

# --- Configuration ---
DEM_VOTES_FILE = 'dem_votes.pkl'
REP_VOTES_FILE = 'rep_votes.pkl'
OUTPUT_FIGURE_MARGIN_FILE = 'district_competitiveness_histogram.png'
OUTPUT_FIGURE_VOTE_SHARE_FILE = 'district_dem_vote_share_histogram.png'
NUM_BINS = 50 # Number of bins for the histogram

# --- Helper function to load data (with error handling) ---
def load_pickle_data(filename):
    """Loads data from a pickle file."""
    try:
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        if not isinstance(data, list):
            raise TypeError(f"Expected a list in {filename}, but got {type(data)}")
        for i, item in enumerate(data):
            if not isinstance(item, defaultdict):
                # Attempt to convert if it's a regular dict
                if isinstance(item, dict):
                    print(f"Warning: Item {i} in {filename} is a dict, converting to defaultdict(int).")
                    data[i] = defaultdict(int, item)
                else:
                    raise TypeError(f"Expected items in the list to be defaultdict (or dict), but got {type(item)} in {filename}")
        return data
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found. Please ensure the file is in the correct directory or provide the full path.")
        return None
    except (pickle.UnpicklingError, EOFError) as e:
        print(f"Error: Could not unpickle data from '{filename}'. File might be corrupted or not a pickle file. Details: {e}")
        return None
    except TypeError as e:
        print(f"Error: Data format incorrect in '{filename}'. {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while loading {filename}: {e}")
        return None

# --- Main script ---
def main():
    print(f"Attempting to load Democratic votes from: {DEM_VOTES_FILE}")
    dem_votes_data = load_pickle_data(DEM_VOTES_FILE)

    print(f"Attempting to load Republican votes from: {REP_VOTES_FILE}")
    rep_votes_data = load_pickle_data(REP_VOTES_FILE)

    if dem_votes_data is None or rep_votes_data is None:
        print("Exiting due to data loading errors. Please check the file paths and file integrity.")
        return

    if len(dem_votes_data) != len(rep_votes_data):
        print("Error: The number of partitionings (entries in the lists) in dem_votes.pkl and rep_votes.pkl do not match.")
        print(f"  Number of Democratic partitionings: {len(dem_votes_data)}")
        print(f"  Number of Republican partitionings: {len(rep_votes_data)}")
        return

    if not dem_votes_data: # Check if the lists are empty
        print("Error: No data found in the pickle files (lists are empty). Cannot generate visualizations.")
        return

    all_margins = []
    all_dem_vote_shares = []

    num_partitionings = len(dem_votes_data)
    print(f"Processing {num_partitionings} partitionings...")

    for i in range(num_partitionings):
        dem_partition_votes = dem_votes_data[i]
        rep_partition_votes = rep_votes_data[i]

        # Ensure they are defaultdicts (should be handled by loader, but as a safeguard)
        if not isinstance(dem_partition_votes, defaultdict):
            dem_partition_votes = defaultdict(int, dem_partition_votes)
        if not isinstance(rep_partition_votes, defaultdict):
            rep_partition_votes = defaultdict(int, rep_partition_votes)

        all_district_ids_in_partition = set(dem_partition_votes.keys()) | set(rep_partition_votes.keys())

        if not all_district_ids_in_partition:
            print(f"Warning: Partitioning {i+1} (index {i}) contains no districts. Skipping.")
            continue

        for district_id in all_district_ids_in_partition:
            dem_count = dem_partition_votes[district_id] # defaultdict provides 0 if key is missing
            rep_count = rep_partition_votes[district_id] # defaultdict provides 0 if key is missing
            total_votes = dem_count + rep_count

            if total_votes > 0:
                margin = (dem_count - rep_count) / total_votes
                dem_vote_share = dem_count / total_votes
                all_margins.append(margin)
                all_dem_vote_shares.append(dem_vote_share)
            # else:
            #    print(f"Info: District '{district_id}' in partitioning {i+1} (index {i}) has zero total votes. Skipping margin calculation for this district.")


    if not all_margins:
        print("Error: No valid district margins or vote shares could be calculated. This might happen if all districts across all partitionings had zero total votes or no districts were found.")
        return

    print(f"Calculated margins for {len(all_margins)} districts across all partitionings.")

    # 1. Plot the Histogram of Margins
    plt.figure(figsize=(12, 7)) # Increased height slightly for better label spacing
    counts, bins, patches = plt.hist(all_margins, bins=NUM_BINS, color='skyblue', edgecolor='black', alpha=0.7)
    plt.title(f'Distribution of Partisan Margins Across All Sampled Districts ({len(all_margins)} total districts from {num_partitionings} partitionings)', fontsize=14)
    plt.xlabel('Partisan Margin ((Dem Votes - Rep Votes) / Total Votes)\n(Negative = Republican Leaning, Positive = Democratic Leaning)', fontsize=12)
    plt.ylabel(f'Number of Districts (Frequency)', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xlim([-1.05, 1.05]) # Give a little space at the ends

    # Add a vertical line at 0 for reference
    plt.axvline(0, color='red', linestyle='dashed', linewidth=1.5, label='Perfectly Competitive (Margin = 0)')
    plt.legend()
    plt.tight_layout() # Adjust layout to prevent labels from overlapping

    try:
        plt.savefig(OUTPUT_FIGURE_MARGIN_FILE)
        print(f"Histogram of partisan margins saved to {OUTPUT_FIGURE_MARGIN_FILE}")
    except Exception as e:
        print(f"Error saving partisan margin histogram: {e}")
    # plt.show() # Comment out if running in a non-interactive environment or saving multiple plots

    # 2. Plot Histogram of Democratic Vote Shares
    plt.figure(figsize=(12, 7))
    plt.hist(all_dem_vote_shares, bins=NUM_BINS, color='lightcoral', edgecolor='black', alpha=0.7)
    plt.title(f'Distribution of Democratic Vote Share Across All Sampled Districts ({len(all_dem_vote_shares)} total districts from {num_partitionings} partitionings)', fontsize=14)
    plt.xlabel('Democratic Vote Share (Dem Votes / Total Votes)', fontsize=12)
    plt.ylabel(f'Number of Districts (Frequency)', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xlim([-0.05, 1.05]) # Give a little space at the ends

    # Add a vertical line at 0.5 for reference
    plt.axvline(0.5, color='blue', linestyle='dashed', linewidth=1.5, label='Split (Vote Share = 0.5)')

    # Add lines for common competitiveness thresholds (e.g., 45% and 55%)
    plt.axvline(0.45, color='darkgrey', linestyle=':', linewidth=1, label='Competitive Range (45%/55%)')
    plt.axvline(0.55, color='darkgrey', linestyle=':', linewidth=1)

    plt.legend()
    plt.tight_layout()

    try:
        plt.savefig(OUTPUT_FIGURE_VOTE_SHARE_FILE)
        print(f"Histogram of Democratic vote shares saved to {OUTPUT_FIGURE_VOTE_SHARE_FILE}")
    except Exception as e:
        print(f"Error saving Democratic vote share histogram: {e}")
    # plt.show() # Display the plot

    print("\n--- Summary ---")
    if all_margins:
        print(f"Average Partisan Margin: {np.mean(all_margins):.3f} (Std: {np.std(all_margins):.3f})")
        median_margin = np.median(all_margins)
        print(f"Median Partisan Margin: {median_margin:.3f}")
        # Count competitive districts, e.g., margin between -0.1 and 0.1 (i.e. 10% margin, or 45-55 split)
        competitive_districts_margin = sum(1 for m in all_margins if -0.1 <= m <= 0.1)
        print(f"Number of districts with margin between -0.1 and 0.1 (approx. 45-55%): {competitive_districts_margin} ({competitive_districts_margin/len(all_margins)*100:.2f}%)")

    if all_dem_vote_shares:
        print(f"Average Democratic Vote Share: {np.mean(all_dem_vote_shares):.3f} (Std: {np.std(all_dem_vote_shares):.3f})")
        median_dem_vote_share = np.median(all_dem_vote_shares)
        print(f"Median Democratic Vote Share: {median_dem_vote_share:.3f}")
        # Count competitive districts, e.g., Dem vote share between 0.45 and 0.55
        competitive_districts_share = sum(1 for s in all_dem_vote_shares if 0.45 <= s <= 0.55)
        print(f"Number of districts with Dem vote share between 0.45 and 0.55: {competitive_districts_share} ({competitive_districts_share/len(all_dem_vote_shares)*100:.2f}%)")

    print(f"\nTo view the plots, open '{OUTPUT_FIGURE_MARGIN_FILE}' and '{OUTPUT_FIGURE_VOTE_SHARE_FILE}'.")
    # **How to Interpret the Visualizations:**

    # * **Concentration around the center (margin approximately 0 or vote share approximately 0.5):** Indicates that many of the randomly sampled districts are competitive.
    # * **Peaks at the extremes (margin approximately -1 or +1; vote share approximately 0 or 1):** Suggests that many districts are heavily skewed towards one party (i.e., "safe" seats).
    # * **Bimodal distribution (peaks at both extremes, dip in the middle):** Can indicate a highly polarized set of districts, where few are competitive.
    # * **Shape of the distribution:** Will give you an overall sense of how competitiveness is distributed across your sampled partitionings.


if __name__ == '__main__':
    main()