import streamlit as st
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.stats.multitest as multitest
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

st.title("HiC Data Analysis - Flemington's Lab")

def median(lst):
    lst = sorted(lst)
    n = len(lst)
    if n % 2 == 1:
        return lst[n // 2]
    else:
        return (lst[n // 2 - 1] + lst[n // 2]) / 2

def collapse_cluster(cluster):
    start1_vals = []
    end1_vals = []
    start2_vals = []
    end2_vals = []

    for loop in cluster:
        split_line = loop.split("\t")
        start1_vals.append(int(split_line[1]))
        end1_vals.append(int(split_line[2]))
        start2_vals.append(int(split_line[4]))
        end2_vals.append(int(split_line[5]))

    start1 = median(start1_vals)
    end1 = median(end1_vals)
    start2 = median(start2_vals)
    end2 = median(end2_vals)

    representative_loop = cluster[0].split("\t")
    representative_loop[1] = str(int(start1))
    representative_loop[2] = str(int(end1))
    representative_loop[4] = str(int(start2))
    representative_loop[5] = str(int(end2))

    return "\t".join(representative_loop)

def process_chunk(chunk, min_distance_away, max_window, min_number_of_proximal_connections, cluster_distance):
    clustered_loops = []
    bedpe_length = len(chunk)

    i = 0
    while i < bedpe_length:
        split_bedpei = chunk[i].split("\t")
        temp_array = [chunk[i]]

        for j in range(i + 1, bedpe_length):
            split_bedpej = chunk[j].split("\t")

            if (split_bedpei[0] == split_bedpej[0] and abs(int(split_bedpei[1]) - int(split_bedpej[1])) > cluster_distance) or \
                    split_bedpei[0] != split_bedpej[0]:
                break

            if (split_bedpei[0] == split_bedpej[0] and split_bedpei[3] == split_bedpej[3] and
                abs(int(split_bedpei[1]) - int(split_bedpej[1])) <= cluster_distance and
                abs(int(split_bedpei[4]) - int(split_bedpej[4])) <= cluster_distance):
                temp_array.append(chunk[j])

        if len(temp_array) > min_number_of_proximal_connections:
            collapsed_cluster = collapse_cluster(temp_array)
            clustered_loops.append(collapsed_cluster)

            i += len(temp_array) - 1

        i += 1

    return clustered_loops

def process_bedpe(file, min_distance_away, max_window, min_number_of_proximal_connections, cluster_distance):
    print(f"Opening file: {file}")
    try:
        with open(file, "r") as inf:
            bedpe = inf.readlines()
    except Exception as e:
        print(f"Error opening file: {e}")
        return []

    print(f"Number of loops read: {len(bedpe)}")

    filtered_bedpe = []
    for line in bedpe:
        split_line = line.strip().split("\t")

        if len(split_line) < 6:  
            continue

        start1 = int(split_line[1])
        start2 = int(split_line[4])
        distance = abs(start2 - start1)

        if distance > max_window:
            continue

        if split_line[0] == "chrEBV_Akata_inverted" and split_line[3] == "chrEBV_Akata_inverted" and \
                (int(split_line[1]) + (171323 - int(split_line[4]))) < min_distance_away:
            continue

        if (split_line[0] == split_line[3] and int(split_line[4]) >= int(split_line[1]) + min_distance_away) or \
                (split_line[0] != split_line[3]):
            filtered_bedpe.append(line.strip())

    print(f"Number of loops after filtering: {len(filtered_bedpe)}")

    if not filtered_bedpe:
        return []

    chunk_size = 1000
    chunks = [filtered_bedpe[i:i + chunk_size] for i in range(0, len(filtered_bedpe), chunk_size)]

    clustered_loops = []
    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_chunk, chunk, min_distance_away, max_window, min_number_of_proximal_connections, cluster_distance): chunk for chunk in chunks}

        for future in as_completed(futures):
            result = future.result()
            clustered_loops.extend(result)

    return clustered_loops

uploaded_file = st.file_uploader("Upload BEDPE File", type="bedpe")

if uploaded_file is not None:
    st.write(f"Processing file: {uploaded_file.name}")

    
    input_filename = os.path.splitext(uploaded_file.name)[0]

    
    temp_file_path = f"./temp_{input_filename}.bedpe"
    with open(temp_file_path, "wb") as temp_file:
        temp_file.write(uploaded_file.getvalue())

    
    data_mutu = pd.read_csv(temp_file_path, delimiter='\t', header=None, names=["chrom_H", "start1", "end1", "EBV", "start2", "end2", "strand1", "strand2", "mapped"])
    
   
    st.write("First few rows of the uploaded file:")
    st.dataframe(data_mutu.head())

   
    bin_size = st.number_input("Bin Size", value=1500)
    genome_size = st.number_input("Genome Size", value=171323)

    
    st.subheader("Right Part Analysis")

    
    data_mutu['bin_right'] = data_mutu['start2'].apply(lambda x: (x - 1) // bin_size + 1)

    
    data_mutu['dups_right'] = data_mutu.groupby(['EBV', 'bin_right'])['bin_right'].transform('count')

   
    data_mutu['dups_right'] /= 10000

   
    st.write("Right Part: Normalized Duplicates:")
    st.dataframe(data_mutu[['start2', 'dups_right']].head())

   
    w_statistic_right, p_value_right = stats.wilcoxon(data_mutu['dups_right'])
    st.write(f"Wilcoxon Test (Right Part) - Statistic: {w_statistic_right}, p-value: {p_value_right}")

    
    _, adjusted_p_values_right, _, _ = multitest.multipletests([p_value_right], alpha=0.05, method='fdr_bh')

    
    mean_right = np.mean(data_mutu['dups_right'])
    std_dev_right = np.std(data_mutu['dups_right'])
    z_scores_right = (data_mutu['dups_right'] - mean_right) / std_dev_right
    mask_right = (adjusted_p_values_right < 0.05) & ((z_scores_right > 2) | (z_scores_right < -2))
    filtered_df_bh_right = data_mutu[mask_right]

   
    st.write("Filtered Data (Right Part):")
    st.dataframe(filtered_df_bh_right)

    
    st.subheader("Left Part Analysis")

   
    data_mutu['bin_left'] = data_mutu['start1'].apply(lambda x: (x - 1) // bin_size + 1)

   
    data_mutu['dups_left'] = data_mutu.groupby(['EBV', 'bin_left'])['bin_left'].transform('count')

    
    data_mutu['dups_left'] /= 10000

    
    st.write("Left Part: Normalized Duplicates:")
    st.dataframe(data_mutu[['start1', 'dups_left']].head())

    
    w_statistic_left, p_value_left = stats.wilcoxon(data_mutu['dups_left'])
    st.write(f"Wilcoxon Test (Left Part) - Statistic: {w_statistic_left}, p-value: {p_value_left}")

    
    _, adjusted_p_values_left, _, _ = multitest.multipletests([p_value_left], alpha=0.05, method='fdr_bh')

    
    mean_left = np.mean(data_mutu['dups_left'])
    std_dev_left = np.std(data_mutu['dups_left'])
    z_scores_left = (data_mutu['dups_left'] - mean_left) / std_dev_left
    mask_left = (adjusted_p_values_left < 0.05) & ((z_scores_left > 2) | (z_scores_left < -2))
    filtered_df_bh_left = data_mutu[mask_left]

    
    st.write("Filtered Data (Left Part):")
    st.dataframe(filtered_df_bh_left)

    
    st.subheader("Combining Left and Right Parts")
    
    
    combined_df = pd.concat([filtered_df_bh_left, filtered_df_bh_right], ignore_index=True)
    
    
    st.write("Combined Data (Left + Right):")
    st.dataframe(combined_df)

    
    combined_file_path = f"Cat_L_R_filtered_{input_filename}_df_bh_final_EBV.bedpe"
    combined_df.to_csv(combined_file_path, sep='\t', index=False, header=False)
    st.write(f"Combined file saved as: '{combined_file_path}'")

    
    with open(combined_file_path, "rb") as file:
        st.download_button(f"Download Combined Data as BEDPE ({combined_file_path})", file, file_name=combined_file_path)

    
    st.subheader("Execute Clustering Process")

    
    min_dist = st.number_input("Minimum Distance", value=1500)
    max_window = st.number_input("Maximum Window", value=146000)
    min_number = st.number_input("Minimum Number", value=25)
    cluster_dist = st.number_input("Cluster Distance", value=2000)

    
    clustered_loops = process_bedpe(combined_file_path, min_dist, max_window, min_number, cluster_dist)

    if not clustered_loops:
        st.error("No clustered loops were generated.")
    else:
        
        output_file = f"{combined_file_path}.min_dist.{min_dist}.max_window.{max_window}.min_number.{min_number}.cluster_dist.{cluster_dist}_pyp_para.bedpe"
        st.write(f"Writing clustered output to: {output_file}")

        
        with open(output_file, "w") as out:
            for loop in clustered_loops:
                out.write(loop + "\n")

        st.success(f"Clustering process completed. Output saved to: {output_file}")

       
        if os.path.exists(output_file):
            with open(output_file, "rb") as clustered_file:
                st.download_button("Download Clustered Data", clustered_file, file_name=os.path.basename(output_file))
