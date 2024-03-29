import shutil
from quart import Quart, request, jsonify, send_file
import asyncio
import subprocess
import os
import json
from aiofiles import open as aio_open
import logging
import sys

# Logging configuration
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s %(name)s %(threadName)s : %(message)s')

app = Quart(__name__)
csv_filename = 'Memory/run.csv'
completed_filename = 'Memory/completed.csv'
script_process = None


async def run_script():
    global script_process
    logging.info("Starting script...")
    try:
        script_process = await asyncio.create_subprocess_exec(
            'python3', 'runMediation.py',
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            cwd=os.getcwd()
        )

        stdout, stderr = await script_process.communicate()
        logging.debug(stdout.decode())
        logging.error(stderr.decode())

        script_process = None
        logging.info("Script completed")
    except Exception as e:
        logging.error(f"Error running script: {e}")


@app.route('/')
async def index():
    return await send_file('index.html')


@app.route('/upload_csv', methods=['POST'])
async def upload_csv():
    files = await request.files
    uploaded_file = files.get('csvfile')
    if uploaded_file:
        filename = uploaded_file.filename.replace(",", "")
        filepath = os.path.join('Input', filename)
        await uploaded_file.save(filepath)
        return jsonify({"message": "File uploaded successfully"}), 200
    return jsonify({"message": "No file uploaded"}), 400


@app.route('/get_csv_file_names', methods=['GET'])
async def get_csv_file_names():
    files = [f for f in os.listdir('Input') if f.endswith('.csv') or f.endswith('.xlsx')]
    return jsonify(files)


@app.route('/add', methods=['POST'])
async def add_row():
    global script_process
    try:
        data = await request.form
        job = [data[field] for field in ['csv_filename', 'trait_name', 'chromosome', 'position', 'flanking', 'minexp', 'percent']]
        logging.debug(f"Job to add: {job}")

        async with aio_open(csv_filename, 'r') as file:
            lines = await file.readlines()
            for line in lines[1:]:
                if line.rstrip().split(",") == job:
                    return jsonify({"message": "Job already exists"}), 400

        async with aio_open(csv_filename, 'a') as file:
            if lines[-1][-1] != '\n':
                await file.write('\n')
            await file.write(','.join(job) + '\n')

        if script_process is None:
            asyncio.create_task(run_script())

        return jsonify({"message": "Job added and row written successfully"}), 200
    except Exception as e:
        logging.error(f"Error in add_row: {e}")
        return jsonify({"message": "Error adding job"}), 500

@app.route('/download_folder', methods=['GET'])
async def download_folder():
    try:
        folder_name = request.args.get('folder_name')
        if folder_name is None:
            return jsonify({"message": "Folder name not found"}), 400

        # Path to the folder
        folder_path = os.path.join('Output', folder_name)

        # Path for the zip file
        zip_path = f"{folder_path}.zip"

        # Create a zip file if it doesn't exist
        if not os.path.exists(zip_path):
            shutil.make_archive(folder_path, 'zip', folder_path)

        # Send the zip file
        return await send_file(zip_path, as_attachment=True)
    except Exception as e:
        logging.error(f"Error in download_folder: {e}")
        return jsonify({"message": "Error downloading folder"}), 500

@app.route('/jobs', methods=['GET'])
async def get_jobs():
    async with aio_open(csv_filename, 'r') as file:
        lines = await file.readlines()
    jobs = [line.rstrip().split(",") for line in lines[1:]]
    return jsonify(jobs)


@app.route('/completed_jobs', methods=['GET'])
async def get_completed_jobs():
    async with aio_open(completed_filename, 'r') as file:
        lines = await file.readlines()
    jobs = [line.rstrip().split(",") for line in lines[1:]]
    return jsonify(jobs)


@app.route('/remove_job', methods=['POST'])
async def remove_job():
    global script_process
    form_data = await request.form
    json_str = list(form_data.keys())[0]
    try:
        data = json.loads(json_str)
        job_to_remove = data.get('job')
        if job_to_remove is None:
            return jsonify({"message": "Job data not found"}), 400
    except json.JSONDecodeError:
        return jsonify({"message": "Invalid JSON format"}), 400

    async with aio_open(csv_filename, 'r') as file:
        lines = await file.readlines()

    if job_to_remove == lines[1].rstrip().split(",") and script_process is not None:
        script_process.terminate()
        script_process = None
        lines = lines[1:]
        asyncio.create_task(run_script())
    else:
        lines = [line for line in lines if line.rstrip().split(",") != job_to_remove]

    async with aio_open(csv_filename, 'w') as file:
        await file.writelines(lines)

    return jsonify({"message": "Job removed successfully"}), 200


@app.route('/download_all_folders', methods=['GET'])
async def download_all_folders():
    try:
        output_dir = 'Output'

        print(os.listdir(output_dir))
        sys.stdout.flush()
        # Remove existing .zip files in the Output directory
        for item in os.listdir(output_dir):
            if item.endswith('.zip'):
                os.remove(os.path.join(output_dir, item))

        # Path for the final zip file
        final_zip_path = 'all_jobs.zip'

        # Create a zip file containing all the contents of the Output folder
        shutil.make_archive(os.path.splitext(final_zip_path)[0], 'zip', output_dir)

        return await send_file(final_zip_path, as_attachment=True)
    except Exception as e:
        logging.error(f"Error in download_all_folders: {e}")
        return jsonify({"message": "Error downloading all folders"}), 500
    

@app.route('/clear_all_jobs', methods=['GET'])
async def clear_all_jobs():
    try:
        # Remove all .zip files and folders in the Output directory
        for item in os.listdir('Output'):
            item_path = os.path.join('Output', item)
            if item.endswith('.zip') or os.path.isdir(item_path):
                if os.path.isdir(item_path):
                    shutil.rmtree(item_path)
                else:
                    os.remove(item_path)

        # Preserve the first line of the completed jobs file
        async with aio_open(completed_filename, 'r') as file:
            first_line = await file.readline()

        async with aio_open(completed_filename, 'w') as file:
            await file.write(first_line)

        return jsonify({"message": "All completed jobs cleared"}), 200
    except Exception as e:
        logging.error(f"Error in clear_all_jobs: {e}")
        return jsonify({"message": "Error clearing jobs"}), 500

@app.route('/is_script_running', methods=['GET'])
async def is_script_running():
    return jsonify({"script_running": script_process is not None})

@app.route('/start_script', methods=['GET'])
async def start_script():
    global script_process
    if script_process is None:
        asyncio.create_task(run_script())
        return jsonify({"message": "Script started successfully"}), 200
    return jsonify({"message": "Script already running"}), 400

@app.route('/stop_script', methods=['GET'])
async def stop_script():
    global script_process
    if script_process is not None:
        script_process.terminate()
        script_process = None
        return jsonify({"message": "Script stopped successfully"}), 200
    return jsonify({"message": "Script not running"}), 400

if __name__ == '__main__':
    app.run()
