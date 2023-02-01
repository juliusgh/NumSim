from PIL import Image
import numpy as np
import math

# Open the image

# resize option
resize = False

#boundary conditions. i = inflow, f = outflow, n = no-slip
top = "n"
bottom = "n"
left = "i"
right = "o"

# set domain name (input and output file name)
domainname = "newyorkbig"

image = Image.open(f"{domainname}.jpg")

image_size = image.size
if resize:
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
else:
    new_size = image_size


resized_image = image.resize(new_size, Image.LANCZOS)

# Convert image to grayscale
gray_image = resized_image.convert("L")




# Convert PIL image to numpy array
gray_image_np = np.array(gray_image)

# Convert grayscale image to binary using numpy where function
bw_image_np = np.where(gray_image_np < 128, "x", " ")


# replace all "x" with with " " where left and right neighours are " "
for i in range(0, bw_image_np.shape[0]):
    for j in range(1, bw_image_np.shape[1]):
        if bw_image_np[i][j] == "x":
            if (j != 0 and j != bw_image_np.shape[1] - 1) and bw_image_np[i][j-1] == " " and bw_image_np[i][j+1] == " ":
                bw_image_np[i][j] = " "
            elif (i != 0 and i != bw_image_np.shape[0] -1) and bw_image_np[i-1][j] == " " and bw_image_np[i+1][j] == " ":
                bw_image_np[i][j] = " "

# Write the 2D binary array to a text file
with open(domainname + ".txt", "w") as text_file:
    text_file.write("n" + "".join((top if bw_image_np[0][i] == " " else "n") for i in range(bw_image_np.shape[1])) + "n" + "\n")
    for row in bw_image_np:
        text_file.write((left if row[0] == " " else "n") + "".join(str(i) for i in row) + (right if row[row.shape[0] -1] == " " else "n") + "\n")
    text_file.write("n" + "".join((bottom if bw_image_np[bw_image_np.shape[0] -1][i] == " " else "n") for i in range(bw_image_np.shape[1] )) + "n" + "\n")

print("File written to: \n" + domainname + ".txt")
print("physical dimensions: \n" + str(new_size[0] / 10 ) + " x " + str(new_size[1] / 10) + " m")
print("nCellsX = ", new_size[0], "\nnCellsY = ", new_size[1])