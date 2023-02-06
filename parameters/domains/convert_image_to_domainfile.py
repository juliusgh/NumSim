from PIL import Image
import numpy as np
import math

# Open the image

# free surface flow
free_surface = False

# resize option
scaling_factor = 40.0

#boundary conditions. i = inflow, o = outflow, n = no-slip
top = "i"
bottom = "i"
left = "i"
right = "i"

# set domain name (input and output file name)
domainname = "catamaran_color"

image = Image.open(f"{domainname}.png")

image_size = image.size

if scaling_factor != 1.0:
    new_size = (int(math.ceil(image_size[0]/(scaling_factor * 10.0))*10), int(math.ceil(image_size[1]/(scaling_factor * 10.0))*10))
else:
    new_size = image_size
resized_image = image.resize(new_size, Image.LANCZOS)

# Convert image to grayscale
gray_image = resized_image.convert("L")

obstacle = "x"
fluid = "-"

if free_surface:
    # Convert PIL image to numpy array with three values: " ", "x" and "a"
    # Define the threshold for the conversion
    lower_threshold = 85
    upper_threshold = 170

    # values assigned to the numpy array
    empty = " "

    # Convert grayscale image to three values using numpy where function
    gray_image_np = np.array(gray_image)
    bw_image_np = np.where(gray_image_np < lower_threshold, obstacle, np.where(gray_image_np > upper_threshold, fluid, empty))
else:
    # Convert PIL image to numpy array
    gray_image_np = np.array(gray_image)

    # Convert grayscale image to binary using numpy where function
    bw_image_np = np.where(gray_image_np < 240, obstacle, fluid)

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

print("physical dimensions: \n" + str(new_size[0] / 10 ) + " x " + str(new_size[1] / 10) + " m")
print("nCellsX = ", new_size[0], "\nnCellsY = ", new_size[1])

# open parameter file template and exchange parameters
lines = []
with open("parameter_template.txt", "r") as template_file:
    lines = template_file.readlines()
for i, line in enumerate(lines):
    if line.startswith("domainFile ="):
        lines[i] = "domainFile = domains/" + domainname + ".txt\n"
    if line.startswith("# ./src/numsim"):
        lines[i] = "# ./src/numsim  ../parameters/" + domainname + "_parameterfile.txt\n"
    if line.startswith("dirichletBottomX") and bottom == "i":
        lines[i] = "dirichletBottomX = 0.0 # here we have inflow\n"
        lines[i+1] = "dirichletBottomY = 0.0 # here we have inflow\n"
    if line.startswith("dirichletTopX") and top == "i":
        lines[i] = "dirichletTopX = 0.0 # here we have inflow\n"
        lines[i+1] = "dirichletTopY = 0.0 # here we have inflow\n"
    if line.startswith("dirichletLeftX") and left == "i":
        lines[i] = "dirichletLeftX = 0.0 # here we have inflow\n"
        lines[i+1] = "dirichletLeftY = 0.0 # here we have inflow\n"
    if line.startswith("dirichletRightX") and right == "i":
        lines[i] = "dirichletRightX = 0.0 # here we have inflow\n"
        lines[i+1] = "dirichletRightY = 0.0 # here we have inflow\n"
with open(f"{domainname}_parameterfile.txt", "w") as parameter_file:
    parameter_file.writelines(lines)

print("Domain file written to: \n" + domainname + ".txt")
print("Parameter file written to: \n" + domainname + "_parameterfile.txt")