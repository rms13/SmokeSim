# Smoke Simulator

ffmpeg usage:
`ffmpeg -framerate 24 -i smoke_%04d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4`

