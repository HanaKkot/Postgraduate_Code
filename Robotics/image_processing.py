import os
import sys
import numpy as np
import cv2




image_size = 64 #64




def resize_image(image, height = image_size, width = image_size):
    return cv2.resize(image,(height, width))

def read_path(path_name, label):
    images = []
    labels = []
    for files in os.listdir(path_name):
        full_path = os.path.join(path_name, files)
        image = resize_image(cv2.imread(full_path, cv2.IMREAD_GRAYSCALE))
        images.append(image)
        labels.append(label)

    return images, labels
