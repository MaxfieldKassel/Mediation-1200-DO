import csv
import subprocess
import os
from datetime import datetime

# Function to run an external R script
script = 'Mediation.R'

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
        status = "Success"
        time_started = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        for i in range(len(row)):
            row[i] = row[i].replace("\"", "")
        folder_name = row[0].split('.')[1] + '_' + row[1].split('.')[1] + '_' + row[2].split('.')[1]
        try:
            # Run the script with parameters from the row
            run_script(row + [folder_name])
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while processing row {row}: {e}")
            status = "Error"

        time_ended = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        time_ran = datetime.strptime(time_ended, "%Y-%m-%d %H:%M:%S") - datetime.strptime(time_started, "%Y-%m-%d %H:%M:%S")

        # Log the completed run
        with open(completed_file, 'a', newline='') as comp_file:
            writer = csv.writer(comp_file)
            writer.writerow([folder_name] + row + [status, time_ran])

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

        


# File paths
input_csv = 'Memory/run.csv'
completed_csv = 'Memory/completed.csv'

# Check if completed.csv exists, if not create it with the header of run.csv
if not os.path.isfile(completed_csv):
    with open(input_csv, 'r') as infile, open(completed_csv, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        writer.writerow(next(reader))  # Copy header

process_file(input_csv, completed_csv)