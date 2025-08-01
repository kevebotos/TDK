#!/bin/bash
# Script to create animation frames using Gmsh command line

echo "Creating animation frames from .pos files..."

# Create frames directory
mkdir -p frames

# Counter for frame numbering
counter=0

# Process each .pos file
for pos_file in temperature_*.pos; do
    if [ -f "$pos_file" ]; then
        frame_name=$(printf "frames/frame_%03d.png" $counter)
        echo "Processing $pos_file -> $frame_name"
        
        # Use Gmsh to export PNG (without -batch option)
        gmsh "$pos_file" -format png -o "$frame_name" -
        
        ((counter++))
    fi
done

echo "Generated $counter frames"

# Check if frames were actually created
if [ -f "frames/frame_000.png" ]; then
    echo "Frames created successfully!"
    
    # Try to create GIF using ImageMagick
    if command -v convert &> /dev/null; then
        echo "Creating GIF animation..."
        convert -delay 50 frames/frame_*.png heat_animation.gif
        echo "Animation saved as heat_animation.gif"
    elif command -v magick &> /dev/null; then
        echo "Creating GIF animation with magick..."
        magick -delay 50 frames/frame_*.png heat_animation.gif
        echo "Animation saved as heat_animation.gif"
    else
        echo "ImageMagick not found. Install with: brew install imagemagick"
    fi
    
    # Try to create MP4 using ffmpeg
    if command -v ffmpeg &> /dev/null; then
        echo "Creating MP4 animation..."
        ffmpeg -y -framerate 2 -i 'frames/frame_%03d.png' -c:v libx264 -pix_fmt yuv420p heat_animation.mp4
        echo "Animation saved as heat_animation.mp4"
    else
        echo "ffmpeg not found. Install with: brew install ffmpeg"
    fi
else
    echo "No frames were created. Gmsh might not support command-line image export."
    echo "Try using the Python script or opening .pos files manually in Gmsh GUI."
fi
