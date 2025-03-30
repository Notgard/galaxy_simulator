import os
import csv

# Define paths to timing directories
SEQ_DIR = "out/out_seq"
OMP_DIR = "out/out_omp"
MPI_DIR = "out/out_mpi"

# Output CSV filenames
SEQ_CSV = "sequential_results.csv"
OMP_CSV = "openmp_results.csv"
MPI_CSV = "mpi_results.csv"


def parse_time_file(file_path):
    """Extract execution times from a file (seconds, milliseconds, microseconds)."""
    with open(file_path, "r") as f:
        content = f.read().strip().split(", ")
    
    if len(content) == 3:
        return float(content[0].replace("s", "")), int(content[1].replace("ms", "")), int(content[2].replace("us", ""))
    return None, None, None


def process_sequential():
    """Process the sequential version and save results to CSV."""
    results = []
    
    for file in sorted(os.listdir(SEQ_DIR)):
        if file.startswith("times_"):
            parts = file.replace("times_", "").replace(".txt", "").split("_")
            nb_particles, n_try = int(parts[0]), int(parts[1])
            seconds, milliseconds, microseconds = parse_time_file(os.path.join(SEQ_DIR, file))
            
            results.append([nb_particles, n_try, seconds, milliseconds, microseconds])

    with open(SEQ_CSV, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["nb_particles", "n_try", "seconds", "milliseconds", "microseconds"])
        writer.writerows(results)


def process_openmp():
    """Process the OpenMP version and save results to CSV."""
    results = []
    
    for file in sorted(os.listdir(OMP_DIR)):
        if file.startswith("times_"):
            parts = file.replace("times_", "").replace(".txt", "").split("_")
            nb_particles, nb_threads, n_try = int(parts[0]), int(parts[1]), int(parts[2])
            seconds, milliseconds, microseconds = parse_time_file(os.path.join(OMP_DIR, file))

            results.append([nb_particles, nb_threads, n_try, seconds, milliseconds, microseconds])

    with open(OMP_CSV, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["nb_particles", "nb_threads", "n_try", "seconds", "milliseconds", "microseconds"])
        writer.writerows(results)


def process_mpi():
    """Process the MPI version (both chunked and subtree) and save results to CSV."""
    results = []
    
    for mpi_type in ["sub_tree", "chunk"]:
        mpi_path = os.path.join(MPI_DIR, "threads_1", mpi_type)
        if not os.path.exists(mpi_path):
            continue
        
        for file in sorted(os.listdir(mpi_path)):
            if file.startswith("times_"):
                parts = file.replace("times_", "").replace(".txt", "").split("_")
                nb_particles, nb_process, n_try = int(parts[0]), int(parts[1]), int(parts[2])
                seconds, milliseconds, microseconds = parse_time_file(os.path.join(mpi_path, file))

                results.append([nb_particles, nb_process, n_try, seconds, milliseconds, microseconds, mpi_type])

    with open(MPI_CSV, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["nb_particles", "nb_process", "n_try", "seconds", "milliseconds", "microseconds", "version"])
        writer.writerows(results)


if __name__ == "__main__":
    process_sequential()
    process_openmp()
    process_mpi()
    print("CSV files generated successfully!")
