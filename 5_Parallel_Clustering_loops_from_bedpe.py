import sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

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
    with open(file, "r") as inf:
        bedpe = inf.readlines()

    print(f"Number of loops read: {len(bedpe)}")

    filtered_bedpe = []
    for line in bedpe:
        split_line = line.strip().split("\t")
        
        # Calculate the distance between start1 and start2
        start1 = int(split_line[1])
        start2 = int(split_line[4])
        distance = abs(start2 - start1)

        # Skip loops where the distance between start1 and start2 exceeds max_window
        if distance > max_window:
            continue

        # Existing filtering logic based on min_distance_away
        if split_line[0] == "chrEBV_Akata_inverted" and split_line[3] == "chrEBV_Akata_inverted" and \
                (int(split_line[1]) + (171323 - int(split_line[4]))) < min_distance_away:
            continue

        if (split_line[0] == split_line[3] and int(split_line[4]) >= int(split_line[1]) + min_distance_away) or \
                (split_line[0] != split_line[3]):
            filtered_bedpe.append(line.strip())

    print(f"Number of loops after filtering: {len(filtered_bedpe)}")

    # Split the data into chunks for parallel processing
    chunk_size = 1000  # Adjust based on system performance
    chunks = [filtered_bedpe[i:i + chunk_size] for i in range(0, len(filtered_bedpe), chunk_size)]

    clustered_loops = []
    
    # Use ProcessPoolExecutor to parallelize the loop processing
    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_chunk, chunk, min_distance_away, max_window, min_number_of_proximal_connections, cluster_distance): chunk for chunk in chunks}

        for future in as_completed(futures):
            result = future.result()
            clustered_loops.extend(result)

    output_file = f"{file}.min_dist.{min_distance_away}.max_window.{max_window}.min_number.{min_number_of_proximal_connections}.cluster_dist.{cluster_distance}_pyp_para.bedpe"
    print(f"Writing output to: {output_file}")
    
    with open(output_file, "w") as out:
        for loop in clustered_loops:
            out.write(loop + "\n")

    print(f"Number of clustered loops written to output: {len(clustered_loops)}")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script.py <file> <min_distance_away> <max_window> <min_number_of_proximal_connections> <cluster_distance>")
        sys.exit(1)

    file = sys.argv[1]
    min_distance_away = int(sys.argv[2])
    max_window = int(sys.argv[3])
    min_number_of_proximal_connections = int(sys.argv[4])
    cluster_distance = int(sys.argv[5])

    process_bedpe(file, min_distance_away, max_window, min_number_of_proximal_connections, cluster_distance)