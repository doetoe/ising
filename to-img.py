#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import PIL.Image

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("infile", help="State file dumped from Ising execution.")
parser.add_argument("outfile", help="Image file visualizing state.", nargs="?")
parser.add_argument("-s", "--show", help="Show image")
args = parser.parse_args()

try: # text format
    b = np.loadtxt(args.infile)
except UnicodeDecodeError: # binary format
    a = np.fromfile(args.infile, dtype=np.int8)
    rows = np.uint8(a[0]) * 256 + np.uint8(a[1])
    cols = np.uint8(a[2]) * 256 + np.uint8(a[3])
    b = a[4:].reshape((rows, cols))
    # plt.imshow(b)
    # plt.show()

rgb = np.moveaxis(
    np.array([b == -1, b == 1, b == 0], dtype=np.uint8) * 255, 0, -1)

def plt_img(rgb):
    fig = plt.figure()
    plt.axis("off")
    plt.imshow(rgb)
    fig.gca().set_aspect(1)
    plt.tight_layout()

if args.show or not args.outfile:
    plt_img(rgb)
    plt.show()

if args.outfile:
    plt_img(rgb)
    plt.savefig(args.outfile)
    PIL.Image.fromarray(rgb).save("pil_" + args.outfile)
