#!/usr/bin/env python3
"""
Create animation from .pos files using matplotlib
This approach reads the .pos files directly and creates visualizations
"""

import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.colors as colors

def parse_pos_file(filename):
    """Parse a Gmsh .pos file to extract coordinates and temperature data"""
    coordinates = []
    temperatures = []
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # Find all SP (scalar point) entries
    pattern = r'SP\(([^,]+),([^,]+),([^)]+)\)\{([^}]+)\}'
    matches = re.findall(pattern, content)
    
    for match in matches:
        x, y, z, temp = match
        coordinates.append([float(x), float(y), float(z)])
        temperatures.append(float(temp))
    
    return np.array(coordinates), np.array(temperatures)

def create_animation_matplotlib():
    """Create animation using matplotlib"""
    
    # Find all .pos files
    pos_files = sorted(glob.glob("temperature_*.pos"))
    
    if not pos_files:
        print("No .pos files found!")
        return
    
    print(f"Found {len(pos_files)} .pos files")
    
    # Read the first file to get the mesh structure
    coords, temps = parse_pos_file(pos_files[0])
    
    if len(coords) == 0:
        print("No data found in .pos files!")
        return
    
    print(f"Found {len(coords)} data points")
    
    # Extract x, y coordinates
    x = coords[:, 0]
    y = coords[:, 1]
    
    # Create triangulation for smooth visualization
    triangulation = tri.Triangulation(x, y)
    
    # Set up the figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Find global temperature range for consistent colorbar
    all_temps = []
    all_coordinates = []
    
    for pos_file in pos_files:
        coords_temp, temps_temp = parse_pos_file(pos_file)
        all_temps.extend(temps_temp)
        all_coordinates.append(coords_temp)
    
    temp_min = min(all_temps)
    temp_max = max(all_temps)
    
    print(f"Temperature range: {temp_min:.1f} to {temp_max:.1f}")
    
    # Create colormap
    norm = colors.Normalize(vmin=temp_min, vmax=temp_max)
    cmap = plt.cm.hot
    
    def animate(frame):
        ax.clear()
        
        # Read current frame data
        coords_frame, temps_frame = parse_pos_file(pos_files[frame])
        
        if len(coords_frame) > 0:
            x_frame = coords_frame[:, 0]
            y_frame = coords_frame[:, 1]
            
            # Create triangulation
            tri_frame = tri.Triangulation(x_frame, y_frame)
            
            # Create the plot
            tcf = ax.tripcolor(tri_frame, temps_frame, cmap=cmap, norm=norm, shading='gouraud')
            
            # Add contour lines
            ax.tricontour(tri_frame, temps_frame, levels=10, colors='black', alpha=0.3, linewidths=0.5)
            
            # Set equal aspect ratio and labels
            ax.set_aspect('equal')
            ax.set_title(f'Heat Transfer Simulation - Time Step {frame}', fontsize=14, fontweight='bold')
            ax.set_xlabel('X Position', fontsize=12)
            ax.set_ylabel('Y Position', fontsize=12)
            
            # Add grid
            ax.grid(True, alpha=0.3)
    
    # Create animation
    anim = FuncAnimation(fig, animate, frames=len(pos_files), interval=500, repeat=True)
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label('Temperature', fontsize=12)
    
    # Save as GIF
    print("Creating animated GIF...")
    writer = PillowWriter(fps=2)
    anim.save('heat_transfer_animation.gif', writer=writer)
    print("Animation saved as heat_transfer_animation.gif")
    
    # Show the animation
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    try:
        create_animation_matplotlib()
    except ImportError as e:
        print(f"Missing required library: {e}")
        print("Install with: pip install matplotlib numpy")
    except Exception as e:
        print(f"Error: {e}")
