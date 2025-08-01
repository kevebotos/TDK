#!/usr/bin/env python3
"""
Script to create an animation from Gmsh .pos files
This script loads all temperature .pos files and exports them as images,
then combines them into an animated GIF or MP4
"""

import gmsh
import os
import glob
import sys
from PIL import Image
import numpy as np

def create_animation():
    # Initialize Gmsh
    gmsh.initialize()
    
    # Find all .pos files
    pos_files = sorted(glob.glob("temperature_*.pos"))
    
    if not pos_files:
        print("No .pos files found!")
        return
    
    print(f"Found {len(pos_files)} .pos files")
    
    # Create output directory for images
    if not os.path.exists("animation_frames"):
        os.makedirs("animation_frames")
    
    images = []
    
    for i, pos_file in enumerate(pos_files):
        print(f"Processing {pos_file}...")
        
        # Clear previous data
        gmsh.clear()
        
        # Open the .pos file
        gmsh.open(pos_file)
        
        # Set up the view
        views = gmsh.view.getTags()
        if views:
            view_tag = views[0]
            
            # Configure the view
            gmsh.option.setNumber("View.ShowScale", 1)
            gmsh.option.setNumber("View.ScaleType", 1)  # Linear scale
            gmsh.option.setNumber("View.RangeType", 2)  # Custom range
            
            # Set consistent color range for all frames
            gmsh.option.setNumber("View.CustomMin", -3000)
            gmsh.option.setNumber("View.CustomMax", 60000)
            
            # Set up graphics
            gmsh.option.setNumber("General.GraphicsWidth", 800)
            gmsh.option.setNumber("General.GraphicsHeight", 600)
            gmsh.option.setNumber("View.Light", 1)
            
            # Export as image
            image_filename = f"animation_frames/frame_{i:03d}.png"
            gmsh.write(image_filename)
            
            # Load image for GIF creation
            if os.path.exists(image_filename):
                img = Image.open(image_filename)
                images.append(img)
    
    # Create animated GIF
    if images:
        print("Creating animated GIF...")
        images[0].save(
            "temperature_animation.gif",
            save_all=True,
            append_images=images[1:],
            duration=500,  # 500ms per frame
            loop=0
        )
        print("Animation saved as temperature_animation.gif")
    
    # Finalize Gmsh
    gmsh.finalize()

if __name__ == "__main__":
    try:
        create_animation()
    except Exception as e:
        print(f"Error: {e}")
        gmsh.finalize()
