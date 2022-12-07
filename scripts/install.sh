sudo apt update
# Install basic packages for building c++ applications
sudo apt install -y build-essential cmake
# Install vtk headers
sudo apt install libvtk9.1 libvtk9-dev
# Install qt dependencies required by vtk
sudo apt install qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools
# Install Google Test suite (optional)
sudo apt install libgtest-dev
