import time
import os
from PIL import Image
from progressBar import *


def compose_frames(id, seq, movie_width, movie_height, matrix_frame_width, mut_matrices_dir, png_dir, composite_dir):
    start_time = time.time()
    total = len(seq) * 19
    printProgressBar(0, total, prefix='Composing final frames:', suffix='Complete', length=50)
    counter = 0
    for i in range(len(seq)):
        for aa in list("ACDEFGHIKLMNPQRSTVWY"):
            if aa == seq[i]:
                continue
            temp_seq = seq.copy()
            temp_seq[i] = aa

            composite_file = os.path.join(composite_dir, "{}.png".format(counter))

            png_file = os.path.join(png_dir, "{}_{}{}{}.png".format(id, seq[i], i + 1, aa))
            img_3d = Image.open(png_file)

            mut_matrix_file = os.path.join(mut_matrices_dir, "{}_matrix_{}{}{}.png".format(id, seq[i], i + 1, aa))
            img_mut_matrix = Image.open(mut_matrix_file)

            tar_img = Image.new('RGB', (movie_width, movie_height), (255, 255, 255))
            tar_img.paste(img_3d, (0, 0))
            tar_img.paste(img_mut_matrix, (movie_width - matrix_frame_width, 0))
            tar_img.save(composite_file)

            counter += 1
            printProgressBar(counter, total, prefix='Composing final frames:', suffix='Complete', length=50)

    end_time = time.time()
    print("Elapsed time: {}".format(end_time - start_time))