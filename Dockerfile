# Start from the Python 3.8 image
FROM python:3.8

# Install R
RUN apt-get update && apt-get install -y r-base

# Set the working directory
WORKDIR /app

# Copy the Python requirements file and install Python dependencies
COPY requirements.txt /app/
RUN pip install -r requirements.txt

# Install R dependencies
RUN R -e "install.packages(c('gplots', 'ggplot2', 'RColorBrewer', 'data.table'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('qtl2', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('qtl2convert', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('parallel', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('readxl', repos='http://cran.rstudio.com/')"

# Copy the rest of the application's code
COPY ./app /app

# Command to run the application
CMD ["hypercorn", "app:app", "--bind", "0.0.0.0:5000"]
