from PIL import Image
import numpy as np
import math

# Open the image

#boundary conditions. i = inflow, f = outflow, n = no-slip
top = "i"
bottom = "n"
left = "n"
right = "n"

domainname = "ferarritest"

image = Image.open("ferrari2.png")

image_size = image.size

# Resize based on image size
if image_size[0] > image_size[1]:
    if 1000 > image_size[0] >= 100:
        new_size = (int(math.ceil(image_size[0]/(10.0 * 10.0))*10), int(math.ceil(image_size[1]/(10.0 * 10.0))*10))
    elif image_size[0] >= 1000:
        new_size = (int(math.ceil(image_size[0]/(100.0 * 10.0))*10), int(math.ceil(image_size[1]/(100.0 * 10.0))*10))
else:
    if 1000 > image_size[1] >= 100:
        new_size = (int(math.ceil(image_size[0]/(10.0 * 10.0))*10), int(math.ceil(image_size[1]/(10.0 * 10.0))*10))
    elif image_size[0] >= 1000:
        new_size = (int(math.ceil(image_size[0]/(100.0 * 10.0))*10), int(math.ceil(image_size[1]/(100.0 * 10.0))*10))


resized_image = image.resize(new_size, Image.LANCZOS)

# Convert image to grayscale
gray_image = resized_image.convert("L")
print("nCellsX = ", new_size[0], "\nnCellsX = ", new_size[1])



# Convert PIL image to numpy array
gray_image_np = np.array(gray_image)

# Convert grayscale image to binary using numpy where function
bw_image_np = np.where(gray_image_np < 128, "o", "x")


# Write the 2D binary array to a text file
with open(domainname + ".txt", "w") as text_file:
    text_file.write("".join(top for i in range(new_size[0] + 2)) + "\n")
    for row in bw_image_np:
        text_file.write(left + "".join(str(i) for i in row) + right + "\n")
    text_file.write("".join(bottom for i in range(new_size[0] + 2)))

#replace symbols with only one neighbour of the same type with the majority neighbour symbol
with open(domainname + ".txt", "r+") as text_file:
    lines = text_file.readlines()
    for i in range(1, len(lines) - 1):
        print(lines[i])
        for j in range(1, len(lines[i]) - 1):
            if lines[i][j] == "x":
                if lines[i][j - 1] != "x" and lines[i][j + 1] != "x" and lines[i+1][j] != "x":
                    lines[i] = lines[:i] + "o" + lines[i + 1:]
                elif lines[i - 1][j] != "x" and lines[i + 1][j] != "x" and lines[i][j+1] != "x":
                    lines[i] = lines[:i] + "o" + lines[i + 1:]
                elif lines[i][j - 1] != "x" and lines[i + 1][j] != "x" and lines[i][j+1] != "x":
                    lines[i] = lines[:i] + "o" + lines[i + 1:]
                elif lines[i][j + 1] != "x" and lines[i + 1][j] != "x" and lines[i][j-1] != "x":
                    lines[i] = lines[:i] + "o" + lines[i + 1:]
with open(domainname + ".txt", "w") as text_file:
    for line in lines:
        text_file.write(line)
# with open(domainname + ".txt", "r") as text_file:
#     lines = text_file.readlines()
#     for i in range(1, len(lines) - 1):
#         for j in range(1, len(lines[i]) - 1):
#             if lines[i][j] == "x":
#                 if lines[i][j - 1] == "x" and lines[i][j + 1] == "x":
#                     lines[i] = lines[i][:j] + "x" + lines[i][j + 1:]
#                 elif lines[i - 1][j] == "x" and lines[i + 1][j] == "x":
#                     lines[i] = lines[i][:j] + "x" + lines[i][j + 1:]
#                 elif lines[i][j - 1] == "x" and lines[i + 1][j] == "x":
#                     lines[i] = lines[i][:j] + "x" + lines[i][j + 1:]
#                 elif lines[i][j + 1] == "x" and lines[i + 1][j] == "x":
#                     lines[i] = lines[i][:j] + "x" + lines[i][j + 1:]
#                 elif lines[i][j - 1] == "x" and lines[i - 1][j] == "x":
#                     lines[i] = lines[i][:j] + "x" + lines[i][j + 1:]
#                 elif lines[i][j + 1] == "x" and lines[i - 1][j] == "x":
#                     lines[i] = lines[i][:j] + "x" + lines[i][j + 1:]    

print("File written to: \n" + domainname + ".txt")
print("physical dimensions: \n" + str(new_size[0] / 10 ) + " x " + str(new_size[1] / 10) + " m")
print("nCellsX = ", new_size[0], "\nnCellsX = ", new_size[1])
print("Note: points surrounded by only one neighbour of the same type have to be manually fixed")