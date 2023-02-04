from PIL import Image
import numpy as np
import math

# Open the image

# free surface flow
free_surface = True

# resize option
scaling_factor = 1.0

#boundary conditions. i = inflow, f = outflow, n = no-slip
top = "i"
bottom = "o"
left = "n"
right = "n"

# set domain name (input and output file name)
domainname = "trinity"

image = Image.open(f"{domainname}.png")

image_size = image.size

new_size = (int(math.ceil(image_size[0]/(scaling_factor * 10.0))*10), int(math.ceil(image_size[1]/(scaling_factor * 10.0))*10))

resized_image = image.resize(new_size, Image.LANCZOS)

# Convert image to grayscale
gray_image = resized_image.convert("L")

if free_surface:
    # Convert PIL image to numpy array with three values: " ", "x" and "a"
    # Define the threshold for the conversion
    lower_threshold = 85
    upper_threshold = 170

    # values assigned to the numpy array
    obstacle = "x"
    water = " "
    air = "a"

    # Convert grayscale image to three values using numpy where function
    gray_image_np = np.array(gray_image)
    bw_image_np = np.where(gray_image_np < lower_threshold, obstacle, np.where(gray_image_np > upper_threshold, water, air))
else:
    # Convert PIL image to numpy array
    gray_image_np = np.array(gray_image)

    # Convert grayscale image to binary using numpy where function
    bw_image_np = np.where(gray_image_np < 128, "x", " ")

# replace all "x" with with " " where left and right neighours are " "
iteration_change = True
while iteration_change:
    iteration_change = False
    for i in range(0, bw_image_np.shape[0]):
        for j in range(1, bw_image_np.shape[1]):
            if bw_image_np[i][j] == "x":
                if (j != 0 and j != bw_image_np.shape[1] - 1) and bw_image_np[i][j-1] == " " and bw_image_np[i][j+1] == " ":
                    bw_image_np[i][j] = " "
                    iteration_change = True
                elif (i != 0 and i != bw_image_np.shape[0] -1) and bw_image_np[i-1][j] == " " and bw_image_np[i+1][j] == " ":
                    bw_image_np[i][j] = " "
                    iteration_change = True

# Write the 2D binary array to a text file
with open(domainname + ".txt", "w") as text_file:
    text_file.write("n" + "".join((top if not bw_image_np[0][i] == "x" else "n") for i in range(bw_image_np.shape[1])) + "n" + "\n")
    for row in bw_image_np:
        text_file.write((left if not row[0] == "x" else "n") + "".join(str(i) for i in row) + (right if not row[row.shape[0] -1] == "x" else "n") + "\n")
    text_file.write("n" + "".join((bottom if not bw_image_np[bw_image_np.shape[0] -1][i] == "x" else "n") for i in range(bw_image_np.shape[1] )) + "n" + "\n")

print("File written to: \n" + domainname + ".txt")
print("physical dimensions: \n" + str(new_size[0] / 10 ) + " x " + str(new_size[1] / 10) + " m")
print("nCellsX = ", new_size[0], "\nnCellsY = ", new_size[1])