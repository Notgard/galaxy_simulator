import os
import csv
import numpy as np

# Define paths to timing directories
SEQ_DIR = "out/out_seq"
OMP_DIR = "out/out_omp"
MPI_DIR = "out/out_mpi/threads_1"
MPI_OMP_DIR = "out/out_mpi_omp/threads_4"
MPI_OMP_128_DIR = "out/out_mpi_omp_128/threads_2"

# Output CSV filenames
SEQ_CSV = "sequential_results.csv"
OMP_CSV = "openmp_results.csv"
MPI_CSV = "mpi_results.csv"
MPI_OMP_CSV = "mpi_omp_results.csv"
MPI_OMP_128_CSV = "mpi_omp_128_results.csv"


def parse_time_file(file_path):
    """Extract execution times from a file (seconds, milliseconds, microseconds)."""
    with open(file_path, "r") as f:
        content = f.read().strip().split(", ")
    
    if len(content) == 3:
        return float(content[0].replace("s", "")), int(content[1].replace("ms", "")), int(content[2].replace("us", ""))
    return None, None, None


def average_times(data):
    """Compute the average of execution times for each unique (nb_particles, *other params*)."""
    averaged_data = {}

    for key, values in data.items():
        print(f"Processing values {key}: {values}")
        avg_seconds = np.mean([v[0] for v in values])
        avg_milliseconds = np.mean([v[1] for v in values])
        avg_microseconds = np.mean([v[2] for v in values])
        averaged_data[key] = [avg_seconds, avg_milliseconds, avg_microseconds]

    return averaged_data


def process_sequential():
    """Process the sequential version and save results to CSV."""
    data = {}

    for file in sorted(os.listdir(SEQ_DIR)):
        if file.startswith("times_"):
            parts = file.replace("times_", "").replace(".txt", "").split("_")
            nb_particles = int(parts[0])
            
            seconds, milliseconds, microseconds = parse_time_file(os.path.join(SEQ_DIR, file))
            key = (nb_particles,)

            if key not in data:
                data[key] = []
            data[key].append([seconds, milliseconds, microseconds])

    averaged_data = average_times(data)

    with open(SEQ_CSV, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["nb_particles", "seconds", "milliseconds", "microseconds"])
        for key, values in averaged_data.items():
            writer.writerow([*key, *values])


def process_openmp():
    """Process the OpenMP version and save results to CSV."""
    data = {}

    for file in sorted(os.listdir(OMP_DIR)):
        if file.startswith("times_"):
            parts = file.replace("times_", "").replace(".txt", "").split("_")
            nb_particles, nb_threads = int(parts[0]), int(parts[1])
            
            seconds, milliseconds, microseconds = parse_time_file(os.path.join(OMP_DIR, file))
            key = (nb_particles, nb_threads)

            if key not in data:
                data[key] = []
            data[key].append([seconds, milliseconds, microseconds])

    averaged_data = average_times(data)

    with open(OMP_CSV, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["nb_particles", "nb_threads", "seconds", "milliseconds", "microseconds"])
        for key, values in averaged_data.items():
            writer.writerow([*key, *values])


def process_mpi(directory, output_file):
    """Process the MPI versions (subtree and chunk) and save results to CSV."""
    data = {}

    for mpi_type in ["sub_tree", "chunk"]:
        mpi_path = os.path.join(directory, mpi_type)
        if not os.path.exists(mpi_path):
            continue

        for file in sorted(os.listdir(mpi_path)):
            if file.startswith("times_"):
                parts = file.replace("times_", "").replace(".txt", "").split("_")
                nb_particles, nb_process = int(parts[0]), int(parts[1])
                
                seconds, milliseconds, microseconds = parse_time_file(os.path.join(mpi_path, file))
                key = (nb_particles, nb_process, mpi_type)

                if key not in data:
                    data[key] = []
                data[key].append([seconds, milliseconds, microseconds])

    averaged_data = average_times(data)

    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["nb_particles", "nb_process", "version", "seconds", "milliseconds", "microseconds"])
        for key, values in averaged_data.items():
            writer.writerow([*key, *values])

if __name__ == "__main__":
    process_sequential()
    process_openmp()
    process_mpi(MPI_DIR, MPI_CSV)
    process_mpi(MPI_OMP_DIR, MPI_OMP_CSV)
    process_mpi(MPI_OMP_128_DIR, MPI_OMP_128_CSV)
    print("CSV files generated successfully!")
