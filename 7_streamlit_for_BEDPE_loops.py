# streamlit run file.py --server.maxUploadSize=204800

import streamlit as st
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.stats.multitest as multitest
import subprocess
import os

# Title of the Web App
st.title("HiC Data Analysis - Flemington's Lab")

# File uploader for the BEDPE file
uploaded_file = st.file_uploader("Upload BEDPE File", type="bedpe")

if uploaded_file is not None:
    # Display file path
    st.write(f"Processing file: {uploaded_file.name}")
    
    # Extract the file name without extension
    input_filename = os.path.splitext(uploaded_file.name)[0]

    # Reading BEDPE file
    data_mutu = pd.read_csv(uploaded_file, delimiter='\t', header=None, names=["chrom_H", "start1", "end1", "EBV", "start2", "end2", "strand1", "strand2", "mapped"])
    
    # Show first few rows of the dataset
    st.write("First few rows of the uploaded file:")
    st.dataframe(data_mutu.head())

    # Parameters
    bin_size = st.number_input("Bin Size", value=1500)
    genome_size = st.number_input("Genome Size", value=171323)

    ### RIGHT PART ###
    st.subheader("Right Part Analysis")

    # Calculate the bin for start2
    data_mutu['bin_right'] = data_mutu['start2'].apply(lambda x: (x - 1) // bin_size + 1)

    # Count duplicates and add the count to a new column "dups_right"
    data_mutu['dups_right'] = data_mutu.groupby(['EBV', 'bin_right'])['bin_right'].transform('count')

    # Normalize each data point by dividing by 10,000
    data_mutu['dups_right'] /= 10000

    # Show the "dups_right" column
    st.write("Right Part: Normalized Duplicates:")
    st.dataframe(data_mutu[['start2', 'dups_right']].head())

    # Perform one-sample Wilcoxon signed-rank test for Right part
    w_statistic_right, p_value_right = stats.wilcoxon(data_mutu['dups_right'])
    st.write(f"Wilcoxon Test (Right Part) - Statistic: {w_statistic_right}, p-value: {p_value_right}")

    # BH Correction for Right part
    _, adjusted_p_values_right, _, _ = multitest.multipletests([p_value_right], alpha=0.05, method='fdr_bh')

    # Filter the DataFrame based on the statistical test for Right part
    mean_right = np.mean(data_mutu['dups_right'])
    std_dev_right = np.std(data_mutu['dups_right'])
    z_scores_right = (data_mutu['dups_right'] - mean_right) / std_dev_right
    mask_right = (adjusted_p_values_right < 0.05) & ((z_scores_right > 2) | (z_scores_right < -2))
    filtered_df_bh_right = data_mutu[mask_right]

    # Show filtered results for Right part
    st.write("Filtered Data (Right Part):")
    st.dataframe(filtered_df_bh_right)

    ### LEFT PART ###
    st.subheader("Left Part Analysis")

    # Calculate the bin for start1
    data_mutu['bin_left'] = data_mutu['start1'].apply(lambda x: (x - 1) // bin_size + 1)

    # Count duplicates and add the count to a new column "dups_left"
    data_mutu['dups_left'] = data_mutu.groupby(['EBV', 'bin_left'])['bin_left'].transform('count')

    # Normalize each data point by dividing by 10,000
    data_mutu['dups_left'] /= 10000

    # Show the "dups_left" column
    st.write("Left Part: Normalized Duplicates:")
    st.dataframe(data_mutu[['start1', 'dups_left']].head())

    # Perform one-sample Wilcoxon signed-rank test for Left part
    w_statistic_left, p_value_left = stats.wilcoxon(data_mutu['dups_left'])
    st.write(f"Wilcoxon Test (Left Part) - Statistic: {w_statistic_left}, p-value: {p_value_left}")

    # BH Correction for Left part
    _, adjusted_p_values_left, _, _ = multitest.multipletests([p_value_left], alpha=0.05, method='fdr_bh')

    # Filter the DataFrame based on the statistical test for Left part
    mean_left = np.mean(data_mutu['dups_left'])
    std_dev_left = np.std(data_mutu['dups_left'])
    z_scores_left = (data_mutu['dups_left'] - mean_left) / std_dev_left
    mask_left = (adjusted_p_values_left < 0.05) & ((z_scores_left > 2) | (z_scores_left < -2))
    filtered_df_bh_left = data_mutu[mask_left]

    # Show filtered results for Left part
    st.write("Filtered Data (Left Part):")
    st.dataframe(filtered_df_bh_left)

    ### COMBINE LEFT AND RIGHT PARTS ###
    st.subheader("Combining Left and Right Parts")
    
    # Combine both left and right parts
    combined_df = pd.concat([filtered_df_bh_left, filtered_df_bh_right], ignore_index=True)
    
    # Show combined data
    st.write("Combined Data (Left + Right):")
    st.dataframe(combined_df)

    # Save the combined data as BEDPE file with input file name appended
    combined_file_path = f"Cat_L_R_filtered_{input_filename}_df_bh_final_EBV.bedpe"
    combined_df.to_csv(combined_file_path, sep='\t', index=False, header=False)
    st.write(f"Combined file saved as: '{combined_file_path}'")

    # Option to download the combined file in BEDPE format
    with open(combined_file_path, "rb") as file:
        st.download_button(f"Download Combined Data as BEDPE ({combined_file_path})", file, file_name=combined_file_path)

    ### EXECUTE EXTERNAL SCRIPT ###
    st.subheader("Execute External Script for Clustering")

    # User-defined parameters for the clustering script
    min_dist = st.number_input("Minimum Distance", value=1500)
    max_window = st.number_input("Maximum Window", value=146000)
    min_number = st.number_input("Minimum Number", value=25)
    cluster_dist = st.number_input("Cluster Distance", value=2000)

    # Running the external script (as a subprocess)
    output_file = f"{combined_file_path}.min_dist.{min_dist}.max_window.{max_window}.min_number.{min_number}.cluster_dist.{cluster_dist}_pyp_para.bedpe"
    command = f"python /Users/truong/Documents/MicroC_2024/1_Script/5_Parallel_Clustering_loops_from_bedpe.py {combined_file_path} {min_dist} {max_window} {min_number} {cluster_dist}"
    
    try:
        subprocess.run(command, shell=True, check=True)
        st.success(f"External clustering script executed successfully! Output saved to: {output_file}")
        
        # Option to download the output file generated by the external script
        if os.path.exists(output_file):
            with open(output_file, "rb") as clustered_file:
                st.download_button("Download Clustered Data", clustered_file, file_name=os.path.basename(output_file))
        else:
            st.error(f"Error: Output file '{output_file}' not found!")
    except subprocess.CalledProcessError as e:
        st.error(f"Error while executing the external script: {str(e)}")

    st.success("Analysis complete for both Wilcoxon Signed-Rank Test and the Benjamini-Hochberg Procedure!\n Truong Nguyen")
