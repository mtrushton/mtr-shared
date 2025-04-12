#!/bin/bash
#
# Creates a Docker container for the use of the Starlink Astronomical Software
# The details are as follows
# - basic ubuntu installation:
# - user USERNAME home directory
# - Starlink astronomical software and dependencies
# - Volume mounting for host directory access
# - hostname for prompt
#
# Notes: 
# To restart container
# sudo docker start -ai $CONTAINER_NAME
# To remove container:
# sudo docker rm $CONTAINER_NAME
# To save the container, create a new image from the running
# modified container
# sudo docker commit $CONTAINER_NAME docker-ubuntu-image
# This can then be used to create a new container 
# sudo docker run -it -v "$PWD/shared_folder:/home/myuser/shared" docker-ubuntu-image
# To remove the image:
# sudo docker rmi $CONTAINER_NAME-image
#
# The X11 installation and forwarding proved problematic. However, the setup outlined
# here works!

# Configs variables (edit as required)
USERNAME="ubuntu_user"
HOST_DIR="$PWD/ubuntu_docker_share"  # shared directory with host system
CONTAINER_NAME="ubuntu-docker-container"
CONTAINER_HOSTNAME="ubuntu-docker"  # This will be used in the prompt

# Create the shared folder if it doesn't exist and assign permissions
echo "Creating shared directory at $HOST_DIR..."
mkdir -p "$HOST_DIR"
chmod 777 "$HOST_DIR" 

# Create Dockerfile
echo "Creating Dockerfile..."
cat > Dockerfile << 'EOF'
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Update and install a basic setup
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    sudo \
    ssh \
    xorg \
    openbox \
    x11-apps \
    vim \
    curl \
    wget \
    git \
    build-essential \
    gfortran \
    libgl1-mesa-glx \
    libx11-dev \
    libxext-dev \
    libpng-dev \
    libz-dev \
    mesa-utils \
    perl \
    tcsh \
    subversion \
    python3 \
    python3-dev \
    python3-pip \
    python3-numpy \
    python3-matplotlib \
    python3-scipy \
    xauth && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
    
# Create a new user with sudo privileges
ARG USERNAME=docker-user
ARG USER_UID=1000
ARG USER_GID=1000

RUN groupadd --gid $USER_GID $USERNAME && \
    useradd --uid $USER_UID --gid $USER_GID -m $USERNAME && \
    echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/$USERNAME && \
    chmod 0440 /etc/sudoers.d/$USERNAME

# Set up a bash prompt
RUN echo 'PS1="\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ "' >> /etc/bash.bashrc

# Set the working directory to the user's home
WORKDIR /home/$USERNAME

# Switch to the new user
USER $USERNAME

# Install Starlink software
RUN mkdir -p /home/$USERNAME/starlink && \
    cd /home/$USERNAME/starlink && \
 #   wget http://star-www.rl.ac.uk/star/starlink-latest.tar.gz && \
     wget https://ftp.eao.hawaii.edu/starlink/2023A/starlink-2023A-Linux-Ubuntu22.tar.gz  && \
     tar xzf starlink-2023A-Linux-Ubuntu22.tar.gz #&& \

    echo 'export STARLINK_DIR=/home/$USERNAME/starlink/star-2023A' >> /etc/bash.bashrc && \
    echo 'source $STARLINK_DIR/etc/profile' >> /etc/bash.bashrc && \
    echo 'echo "Starlink software is ready to use"' >> /etc/bash.bashrc

# Add Starlink variables to bashrc
RUN echo "export STARLINK_DIR=/home/$USERNAME/starlink/star-2023A" >> /home/$USERNAME/.bashrc && \
   echo "source \$STARLINK_DIR/etc/profile" >> /home/$USERNAME/.bashrc

# Set the default command to bash
CMD ["/bin/bash"]
EOF

# Build Docker image
echo "Building Docker image..."
sudo docker build -t "$CONTAINER_NAME-image" \
  --build-arg USERNAME="$USERNAME" \
  --build-arg USER_UID="$(id -u)" \
  --build-arg USER_GID="$(id -g)" .

# Remove container if it already exists
echo "Checking for existing container..."
if sudo docker ps -a --format '{{.Names}}' | grep -q "^$CONTAINER_NAME$"; then
    echo "Container $CONTAINER_NAME already exists. Removing it..."
    sudo docker rm -f "$CONTAINER_NAME"
fi

# Allow X11 forwarding to the host
xhost +local:docker

# Run the container with the shared volume, custom hostname, and pass the DISPLAY 
# environment variable to the container and mount X11 socket from host to container
echo "Starting Docker container with shared volume, custom hostname..."
sudo docker run -it \
  --name "$CONTAINER_NAME" \
  --hostname "$CONTAINER_HOSTNAME" \
  --network=host \
  -e DISPLAY=$DISPLAY \
  -e XAUTHORITY=$HOME/.Xauthority \
  -v $HOME/.Xauthority:/home/$USERNAME/.Xauthority:ro \
  -v /tmp/.X11-unix:/tmp/.X11-unix \
  -v "$HOST_DIR:/home/$USERNAME/shared" \
  "$CONTAINER_NAME-image"
  

  
