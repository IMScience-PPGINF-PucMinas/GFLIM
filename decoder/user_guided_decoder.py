import pyift.pyift as ift
import skimage
import os
import cv2
import numpy as np
import csv
import matplotlib.pyplot as plt
from skimage.filters import threshold_otsu
import math
import argparse
import sys
from tqdm import tqdm


def normalize_by_band(feature):
    max_ = feature.max()
    min_ = feature.min()
    norm = 255 * ((feature - min_) / (max_ - min_))
    
    return norm

def show_two_fig(image1, image2):
    fig = plt.figure(figsize=(20,20))
    fig.add_subplot(1, 3, 1)
    plt.imshow(image1, cmap="gray")
    plt.axis('off')
    plt.title("First")
    
    fig.add_subplot(1, 3, 2)
    plt.imshow(image2, cmap="gray")
    plt.axis('off')
    plt.title("Second")
    
def decoder_fb(feature, kernel_importances = None, threshold=0.5, filter_by_size=False, size_range=[2000,15000], foreground_kernels=None):
    kernel  = np.ones((5,5),np.uint8)

    if(kernel_importances is None):
        kernel_importances = np.ones(feature.shape[3])
        
    img_size = [feature.shape[1], feature.shape[2]]

    positives = np.zeros(img_size)
    negatives = np.zeros(img_size)

    counter = 0
    n_counter = 0
    
    for b in range(feature.shape[3]):
        if(feature[0,:,:,b].max() != 0.0):
            band_image = feature[0,:,:,b]
            band_mean = (band_image/band_image.max()).mean()
            norm_band_image = normalize_by_band(band_image)

            if(foreground_kernels is None):
                if(band_mean < threshold):
                    positives += (band_image*kernel_importances[b])
                    counter+=kernel_importances[b]
                else:
                    negatives += (band_image*kernel_importances[b])
                    n_counter+=kernel_importances[b]
            else:
                if(b in foreground_kernels):
                    positives += (band_image*kernel_importances[b])
                    counter+=kernel_importances[b]
                else:
                    negatives += (band_image*kernel_importances[b])
                    n_counter+=kernel_importances[b]

    if(counter != 0):
        positives/=counter
        positives = normalize_by_band(positives)
    if(n_counter != 0):
        negatives/=n_counter
        negatives = normalize_by_band(negatives)
       
    if(not math.isnan(negatives.max())):
        decoded = positives - negatives
    else:
        decoded = positives
        
    
    decoded[decoded < 0] = 0
    decoded = decoded.astype(np.uint8)
    
    if(filter_by_size):
        nb_components, components_image, stats, centroids = cv2.connectedComponentsWithStats((decoded>0).astype(np.uint8), connectivity=8)
        for c in range(1,nb_components):
            indices = np.argwhere(components_image==c)
            if(indices.shape[0] < size_range[0] or indices.shape[0] > size_range[1]):
                for i in indices:
                    decoded[i[0], i[1]] = 0
                
    if(decoded.max() == 0):
        return decoded
    decoded = normalize_by_band(decoded)
    
    return decoded
    
if __name__ == "__main__":
    try:
        ap = argparse.ArgumentParser()
    except:
        ap.print_help()
        sys.exit(0)
    ap.add_argument("-i", "--input_dataset", required=True,	help="path to the folder with the feature images")
    ap.add_argument("-o", "--output_dataset", required=True,	help="path to the folder for the decoded images")
    ap.add_argument("-t", "--threshold", required=False, default=0.5, help="Threshold to determine based on the mean saliency if the image is to be considered foreground")
    ap.add_argument("-fg", "--foreground", required=False, default=None, help="A list of foreground kernel numbers")
    
    args = vars(ap.parse_args())
    
    input_dataset = args["input_dataset"]
    threshold = float(args["threshold"])
    output_dataset = args["output_dataset"]

    foreground_kernels = [int(x) for x in args["foreground"].split(",")]
    
    feature_ext = "mimg"
    out_ext = "png"
    
    if not os.path.exists(output_dataset):
        os.makedirs(output_dataset)
    
    dataset_size = len(os.listdir(input_dataset))
    toolbar_width = dataset_size
    print("Decoding images")
    dataset_files = os.listdir(input_dataset)

    for i in tqdm(range(len(dataset_files))):

        filename = dataset_files[i]
        filename_,ext = filename.split(".")
        featuref = os.path.join(input_dataset, filename)
        mfeature = ift.ReadMImage(featuref)
        
        feature = mfeature.AsNumPy()
        decoded = decoder_fb(feature, threshold=threshold, foreground_kernels=foreground_kernels)
        
        cv2.imwrite(output_dataset+"/"+filename_+"."+out_ext, decoded)
    sys.stdout.write("]\n") # this ends the progress bar
    print("Done!")
