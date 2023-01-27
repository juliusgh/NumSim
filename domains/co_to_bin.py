from PIL import Image
import numpy as np
import math

# Open the image

#boundary conditions. i = inflow, f = outflow, n = no-slip
top = "i"
bottom = "n"
left = "n"
right = "n"

domainname = "vase"

image = Image.open("vase_sw.png")

image_size = image.size

# Resize based on image size
if image_size[0] > image_size[1]:
    if 1000 > image_size[0] >= 100:
        new_size = (int(math.ceil(image_size[0]/(10.0 * 10.0))*10), int(math.ceil(image_size[1]/(10.0 * 10.0))*10))
    elif image_size[0] >= 1000:
        new_size = (int(math.ceil(image_size[0]/(100.0 * 10.0))*10), int(math.ceil(image_size[1]/(100.0 * 10.0))*10))
    else: new_size = image_size
else:
    if 1000 > image_size[1] >= 100:
        new_size = (int(math.ceil(image_size[0]/(10.0 * 10.0))*10), int(math.ceil(image_size[1]/(10.0 * 10.0))*10))
    elif image_size[0] >= 1000:
        new_size = (int(math.ceil(image_size[0]/(100.0 * 10.0))*10), int(math.ceil(image_size[1]/(100.0 * 10.0))*10))
    else: new_size = image_size

resized_image = image.resize(new_size, Image.LANCZOS)

# Convert image to grayscale
gray_image = resized_image.convert("L")




# Convert PIL image to numpy array
gray_image_np = np.array(gray_image)

# Convert grayscale image to binary using numpy where function
bw_image_np = np.where(gray_image_np < 128, "x", " ")

print(bw_image_np.shape)
print(bw_image_np.shape[0])
print(bw_image_np.shape[1])

# replace all "x" with with " " where left and right neighours are " "
for i in range(1, bw_image_np.shape[0] - 2):
    for j in range(1, bw_image_np.shape[1] - 2):
        if bw_image_np[i][j] == "x":
            if bw_image_np[i][j-1] == " " and bw_image_np[i][j+1] == " ":
                bw_image_np[i][j] = " "
            elif bw_image_np[i-1][j] == " " and bw_image_np[i+1][j] == " ":
                bw_image_np[i][j] = " "

# Write the 2D binary array to a text file
with open(domainname + ".txt", "w") as text_file:
    text_file.write("".join(top for i in range(new_size[0] + 2)) + "\n")
    for row in bw_image_np:
        text_file.write(left + "".join(str(i) for i in row) + right + "\n")
    text_file.write("".join(bottom for i in range(new_size[0] + 2)) + "\n")

print("nCellsX = ", new_size[0], "\nnCellsX = ", new_size[1])
print("File written to: \n" + domainname + ".txt")
print("physical dimensions: \n" + str(new_size[0] / 10 ) + " x " + str(new_size[1] / 10) + " m")
print("nCellsX = ", new_size[0], "\nnCellsY = ", new_size[1])
print("Note: points surrounded by only one neighbour of the same type have to be manually fixed")