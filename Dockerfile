# Use an official Python runtime as a parent image
FROM python:3.8-slim-buster

# Set the working directory in the container
WORKDIR /usr/src/app

# Copy the current directory contents into the container at /usr/src/app
COPY . .

RUN pip install --upgrade pip

# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt


# Install dependencies required for building vsearch
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    autoconf \
    automake \
    libtool \
    zlib1g-dev \
    make \
    gcc \
    zlib1g-dev \
    nano \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Clone the minimap2 repository
RUN git clone https://github.com/lh3/minimap2.git

# Change to the minimap2 directory
WORKDIR /usr/src/app/minimap2

# Build minimap2
RUN make

# Install minimap2
RUN cp minimap2 /usr/local/bin/

# Reset the working directory
WORKDIR /usr/src/app

# Clone the vsearch repository
RUN git clone https://github.com/torognes/vsearch.git

# Change to the vsearch directory
WORKDIR /usr/src/app/vsearch

# Build vsearch
RUN ./autogen.sh \
    && ./configure \
    && make \
    && make install

# Reset the working directory
WORKDIR /usr/src/app

# Download the 32-bit version of USEARCH
RUN curl -L -o usearch32.gz http://drive5.com/downloads/usearch8.1.1861_i86linux32.gz

RUN gunzip usearch32

RUN cp usearch32 /usr/local/bin/

RUN chmod +x /usr/local/bin/usearch32

# Run the Python script when the container launches
# CMD ["python3", "./src/run_steps.py"]
ENTRYPOINT [ "/bin/bash" ] 