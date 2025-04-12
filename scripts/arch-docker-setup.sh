#!/bin/bash
#
# Creates a Docker container for a basic arch installation with:
# - Custom user with home directory
# - Required packages installed
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
# sudo docker commit $CONTAINER_NAME docker-arch-image
# This can then be used to create a new container 
# sudo docker run -it -v "$PWD/shared_folder:/home/myuser/shared" docker-arch-image
# To remove the image:
# sudo docker rmi $CONTAINER_NAME-image


# Configs variables (edit as required)
USERNAME="arch_user"
HOST_DIR="$PWD/arch_docker_share"  # shared directory with host system
CONTAINER_NAME="arch-docker-container"
CONTAINER_HOSTNAME="arch-docker"  # This will be used in the prompt

# Create the shared folder if it doesn't exist and assign permissions
echo "Creating shared directory at $HOST_DIR..."
mkdir -p "$HOST_DIR"
chmod 777 "$HOST_DIR" 

# Create Dockerfile
echo "Creating Dockerfile..."
cat > Dockerfile << 'EOF'
FROM archlinux:latest

# Update and install a basic setup
RUN pacman -Syu --noconfirm && \
    pacman -S --noconfirm \
    sudo \
    openssh \
    xorg \
    openbox \
    vim \
    curl \
    wget \
    git && \
    pacman -Scc --noconfirm

# Create a new user with sudo privileges
ARG USERNAME=docker-user
ARG USER_UID=1000
ARG USER_GID=1000

RUN groupadd --gid $USER_GID $USERNAME && \
    useradd --uid $USER_UID --gid $USER_GID -m $USERNAME && \
    echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/$USERNAME && \
    chmod 0440 /etc/sudoers.d/$USERNAME

# Set up a nicer bash prompt
RUN echo 'PS1="\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ "' >> /etc/bash.bashrc

# Set the working directory to the user's home
WORKDIR /home/$USERNAME

# Switch to the new user
USER $USERNAME

# Set the default command to bash
CMD ["/bin/bash"]
EOF

# Build Docker image
echo "Building Docker image..."
sudo docker build -t "$CONTAINER_NAME-image" \
  --build-arg USERNAME="$USERNAME" \
  --build-arg USER_UID="$(id -u)" \
  --build-arg USER_GID="$(id -g)" .

# Run the container with the shared volume and custom hostname
echo "Starting Docker container with shared volume and custom hostname..."
sudo docker run -it \
  --name "$CONTAINER_NAME" \
  --hostname "$CONTAINER_HOSTNAME" \
  -v "$HOST_DIR:/home/$USERNAME/shared" \
  "$CONTAINER_NAME-image"
