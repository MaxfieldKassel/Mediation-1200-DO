FROM python:3.8

# Set the working directory
WORKDIR /app

# Copy the requirements file and install dependencies
COPY requirements.txt /app/
RUN pip install -r requirements.txt

# Copy the rest of the application's code
COPY . /app

# Command to run the application
CMD ["hypercorn", "app:app", "--bind", "0.0.0.0:5000"]
