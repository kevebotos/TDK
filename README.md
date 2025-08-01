# TDK - 2D Transient Heat Transfer Simulation

A complete C++ finite element method (FEM) implementation for simulating 2D transient heat transfer on a stove surface using Gmsh and Eigen libraries.

## 🔥 Features

- **Real-time heat transfer simulation** with 4 active burners
- **Different material properties** for stove body vs. burner elements
- **Realistic temperature ranges** (20°C to ~114°C)
- **Animated visualization** of temperature distribution over time
- **Gmsh mesh integration** with physical groups
- **Eigen sparse matrix** computations for efficiency

## 📋 Prerequisites

### Required Software

- **macOS** with Homebrew
- **Xcode Command Line Tools**
- **CMake** (3.15 or higher)
- **Git**

### Required Libraries

- **Gmsh** (4.14.0+) - Mesh processing and visualization
- **Eigen** (3.4.0+) - Linear algebra operations
- **Python** (3.8+) with matplotlib, numpy, pillow

## 🚀 Quick Start

### 1. Clone and Setup

```bash
# Clone the repository
git clone https://github.com/kevebotos/TDK.git
cd TDK
```

### 2. Install Dependencies

```bash
# Install required libraries via Homebrew
brew install eigen gmsh

# Set up Python environment and install packages
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install matplotlib numpy pillow
```

### 3. Build the Project

```bash
# Navigate to cmake-build-debug directory (or create build directory)
cd cmake-build-debug

# If using custom build directory instead:
# mkdir -p build && cd build && cmake ..

# Build the executable
make

# Verify build
ls -la tdk  # Should show the executable
```

### 4. Run the Simulation

```bash
# Copy mesh file to build directory
cp ../src/stove.msh .

# Run the heat transfer simulation
./tdk
```

**Expected output:**

```
Successfully opened stove.msh
Number of nodes: 9840
Number of triangular elements: 19318
Stove Body nodes: 8404
Burner 1 nodes: 359
...
Simulation completed successfully!
Generated 31 output files
```

### 5. Create Animation

```bash
# Activate Python environment (if not already active)
source ../.venv/bin/activate

# Create animated GIF from simulation results
python ../create_matplotlib_animation.py
```

## 📊 Understanding the Output

### Generated Files

- **`temperature_X.pos`** - Gmsh post-processing files (31 files, every 10 seconds)
- **`heat_transfer_animation.gif`** - Animated visualization of temperature evolution

### File Structure After Simulation

```
cmake-build-debug/  (or build/)
├── tdk                           # Main executable
├── stove.msh                     # Input mesh file
├── temperature_0.pos             # Initial temperature (t=0s)
├── temperature_1.pos             # Temperature at t=10s
├── ...                           # ... continuing every 10s
├── temperature_30.pos            # Final temperature (t=300s)
└── heat_transfer_animation.gif   # Complete animation
```

## 🎯 Detailed Command Reference

### Building Commands

```bash
# Clean previous build
make clean

# Full rebuild
make clean && make

# Parallel build (faster)
make -j$(nproc)

# Configure CMake with debug info (if using custom build directory)
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Configure CMake with release optimization
cmake -DCMAKE_BUILD_TYPE=Release ..
```

### Running Commands

```bash
# Basic simulation run
./tdk

# Run with output redirection
./tdk > simulation_log.txt 2>&1

# Run and time execution
time ./tdk

# Run in background
./tdk &
```

### Analysis Commands

```bash
# Count generated files
ls *.pos | wc -l

# Show file sizes
ls -lh *.pos *.gif

# Check latest temperatures
tail -10 simulation_log.txt

# View a specific .pos file in Gmsh GUI
gmsh temperature_10.pos

# Open all .pos files for animation in Gmsh
gmsh temperature_*.pos
```

### Animation Commands

```bash
# Create animation (basic)
python ../create_matplotlib_animation.py

# Create animation with custom Python interpreter
/path/to/python ../create_matplotlib_animation.py

# Check animation file size and properties
ls -lh *.gif
file heat_transfer_animation.gif
```

### File Management Commands

```bash
# Remove all generated files
rm -f *.pos *.gif

# Remove just animation files
rm -f *.gif

# Archive simulation results
tar -czf simulation_results.tar.gz *.pos *.gif

# View file contents (first .pos file)
head -20 temperature_0.pos
```

## 🛠️ Development Commands

### VS Code Integration

```bash
# Open project in VS Code
code ..

# Build using VS Code task (Cmd+Shift+P → "Tasks: Run Task")
# Select "C/C++: clang build active file"

# Debug in VS Code (F5 or Run → Start Debugging)
```

### Debugging

```bash
# Build with debug symbols
make clean && make

# Run with lldb debugger
lldb ./tdk
(lldb) run
(lldb) bt  # Show backtrace if crash occurs
(lldb) quit
```

### Testing and Validation

```bash
# Check if libraries are properly linked
otool -L tdk

# Verify Gmsh can read mesh
gmsh ../src/stove.msh -info

# Test Python environment
python -c "import matplotlib, numpy; print('✅ Libraries OK')"

# Check mesh file details
head -20 ../src/stove.msh
```

## 📈 Simulation Parameters

### Current Configuration

- **Total time**: 300 seconds (5 minutes)
- **Time step**: 1 second
- **Output interval**: Every 10 seconds (31 frames)
- **Thermal diffusivity**:
  - Stove body: 4.2×10⁻⁶ m²/s (stainless steel)
  - Burners: 1.2×10⁻⁵ m²/s (cast iron)
- **Heat power**: 2kW per burner (realistic)

### Customizing Parameters

Edit `src/main.cpp` and modify:

```cpp
double total_time = 300.0;      // Simulation duration
double dt = 1.0;                // Time step size
int output_interval = 10;       // Output frequency
double heat_power_per_burner = 2000.0;  // Heat source power
```

## 🔍 Troubleshooting

### Common Issues and Solutions

#### Build Errors

```bash
# If Eigen not found
brew reinstall eigen

# If Gmsh not found
brew reinstall gmsh

# If CMake configuration fails in custom build
rm -rf build && mkdir build && cd build && cmake ..

# If make fails in cmake-build-debug
cd cmake-build-debug && make clean && make
```

#### Runtime Errors

```bash
# If mesh file not found
cp ../src/stove.msh .

# If permission denied
chmod +x tdk

# If segmentation fault
lldb ./tdk
(lldb) run
```

#### Animation Issues

```bash
# If Python packages missing
pip install --upgrade matplotlib numpy pillow

# If animation script fails
python -c "import sys; print(sys.path)"  # Check Python path

# If virtual environment issues
deactivate
python3 -m venv .venv
source .venv/bin/activate
pip install matplotlib numpy pillow
```

### Performance Tips

- Use `make -j$(nproc)` for parallel compilation
- For large meshes, increase available memory
- Use release build for faster execution: `cmake -DCMAKE_BUILD_TYPE=Release ..`

## 📁 Project Structure

```
TDK/
├── src/
│   ├── main.cpp              # Main simulation code
│   └── stove.msh             # Gmsh mesh file
├── cmake-build-debug/        # Build directory (CLion default)
├── .vscode/                  # VS Code configuration
├── CMakeLists.txt            # CMake configuration
├── README.md                 # This file
├── create_matplotlib_animation.py  # Animation script
└── .venv/                    # Python virtual environment
```

## 🎓 Learning Resources

### Understanding the Physics

- **Heat equation**: ∂T/∂t = α∇²T + Q
- **Finite Element Method**: Discretization using triangular elements
- **Backward Euler**: Implicit time integration for stability

### Key Concepts

- **Thermal diffusivity (α)**: How quickly heat spreads through material
- **Heat sources**: Applied as Gaussian distributions at burner centers
- **Material properties**: Different α values for stove body vs. burners

### Complete Workflow Summary

```bash
# 1. Navigate to build directory
cd cmake-build-debug

# 2. Copy mesh file
cp ../src/stove.msh .

# 3. Build project
make clean && make

# 4. Run simulation
./tdk

# 5. Activate Python environment
source ../.venv/bin/activate

# 6. Create animation
python ../create_matplotlib_animation.py

# 7. View results
ls -lh *.gif *.pos
```

## 📄 License

This project is for educational purposes. Feel free to modify and experiment!

---

**Happy simulating! 🔥**

_Temperature ranges from 20°C (room temperature) to ~114°C (boiling point vicinity)_
