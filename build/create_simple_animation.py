#!/usr/bin/env python3
"""
Simple animation script using Gmsh Python API
Creates PNG frames and combines them into a GIF
"""

import gmsh
import os
import glob
import subprocess

def create_simple_animation():
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    
    # Find all .pos files
    pos_files = sorted(glob.glob("temperature_*.pos"))
    print(f"Found {len(pos_files)} .pos files")
    
    # Create frames directory
    os.makedirs("frames", exist_ok=True)
    
    # Process each file
    for i, pos_file in enumerate(pos_files):
        print(f"Processing {pos_file} -> frame_{i:03d}.png")
        
        # Clear and load file
        gmsh.clear()
        gmsh.open(pos_file)
        
        # Set graphics options
        gmsh.option.setNumber("General.GraphicsWidth", 1024)
        gmsh.option.setNumber("General.GraphicsHeight", 768)
        
        # Configure view if available
        views = gmsh.view.getTags()
        if views:
            # Set consistent color scale
            gmsh.view.option.setNumber(views[0], "RangeType", 2)  # Custom range
            gmsh.view.option.setNumber(views[0], "CustomMin", -3000)
            gmsh.view.option.setNumber(views[0], "CustomMax", 60000)
            gmsh.view.option.setNumber(views[0], "ShowScale", 1)
        
        # Export frame
        frame_path = f"frames/frame_{i:03d}.png"
        gmsh.write(frame_path)
    
    gmsh.finalize()
    
    # Create GIF using ImageMagick (if available)
    try:
        print("Creating GIF animation...")
        cmd = ["convert", "-delay", "50", "frames/frame_*.png", "heat_transfer_animation.gif"]
        subprocess.run(cmd, check=True)
        print("Animation saved as heat_transfer_animation.gif")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("ImageMagick not found. Install it with: brew install imagemagick")
        print("Or manually combine the PNG files in the 'frames' directory")

if __name__ == "__main__":
    create_simple_animation()
