# Use an official Python runtime as a parent image
FROM python:3

# Don't change this section. It ensures docker can find requirements.txt
# Set the working directory to /app
WORKDIR /app
# Copy the current directory contents into the container at /app
COPY . /app

RUN apt-get update -q && apt-get install -qy \
    texlive-full \
    python-pygments gnuplot \
    make git \
    && rm -rf /var/lib/apt/lists/*

# Install any needed packages specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt

# Define environment variable -- example only
ENV NAME World

# Run bash when the container launches
CMD ["bash"]