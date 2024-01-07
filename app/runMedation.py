import csv
import subprocess
import os

# Function to run an external R script
script = 'RScript/Mediation.R'

def run_script(params):
    command = ['Rscript', script] + params

    # Get the current working directory
    current_directory = os.getcwd()

    # Run the script in the current directory
    subprocess.run(command, check=True, cwd=current_directory)


def process_file(input_file, completed_file):
    with open(input_file, 'r') as file:
        reader = csv.reader(file)
        rows = list(reader)

    for row in rows[1:]:  # Skipping the header
        try:
            # Run the script with parameters from the row
            run_script(row)

            # Log the completed run
            with open(completed_file, 'a', newline='') as comp_file:
                writer = csv.writer(comp_file)
                writer.writerow(row)

            # Re-read the file to get the most updated list
            with open(input_file, 'r') as file:
                reader = csv.reader(file)
                updated_rows = list(reader)

            # Remove the processed line from the updated list
            updated_rows.remove(row)

            # Rewrite the input file without the processed line
            with open(input_file, 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerows(updated_rows)

        except subprocess.CalledProcessError as e:
            print(f"An error occurred while processing row {row}: {e}")


# File paths
input_csv = 'run.csv'
completed_csv = 'completed.csv'

# Check if completed.csv exists, if not create it with the header of run.csv
if not os.path.isfile(completed_csv):
    with open(input_csv, 'r') as infile, open(completed_csv, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        writer.writerow(next(reader))  # Copy header

process_file(input_csv, completed_csv)