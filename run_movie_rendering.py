import sys
import os
import time
import numpy as np
import shutil
from Bio import SeqIO
from PIL import Image, ImageFont, ImageDraw
from progressBar import *
from render_3d_frames import render_3d_frames
from render_mutation_matrices import render_mutation_matrices


def main():
    input_fasta = sys.argv[1]
    pdb_dir = sys.argv[2]
    
    movie_width = 1280
    movie_height = 720
    matrix_frame_width = 250
    
    tmp_dir = "./tmp/"
    
    if len(sys.argv) > 3:
        movie_width = int(sys.argv[3])
    if len(sys.argv) > 4:
        movie_height = int(sys.argv[4])
    if len(sys.argv) > 5:
        tmp_dir = int(sys.argv[4])

    aa_list = list("ACDEFGHIKLMNPQRSTVWY")
    for record in SeqIO.parse(input_fasta, "fasta"):
        id = record.id
        seq = list(record.seq)

        print("Processing {}".format(id))

        png_dir = os.path.join(tmp_dir, "png")
        mut_matrices_dir = os.path.join(tmp_dir, "mut_matrices_png")
        composite_dir = os.path.join(tmp_dir, "composite_png")
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)
        if not os.path.isdir(png_dir):
            os.makedirs(png_dir)
        if not os.path.isdir(mut_matrices_dir):
            os.makedirs(mut_matrices_dir)
        if not os.path.isdir(composite_dir):
            os.makedirs(composite_dir)
            
        render_3d_frames(id, seq, pdb_dir, png_dir, movie_width - matrix_frame_width, movie_height)
        render_mutation_matrices(id, seq, movie_height, pdb_dir, mut_matrices_dir)

        start_time = time.time()
        total = len(seq) * 19
        printProgressBar(0, total, prefix = 'Composing final frames:', suffix = 'Complete', length = 50)
        counter = 0
        for i in range(len(seq)):
            for aa in aa_list:
                if aa == seq[i]:
                    continue
                temp_seq = seq.copy()
                temp_seq[i] = aa

                composite_file = os.path.join(composite_dir, "{}.png".format(counter))
               
                png_file = os.path.join(png_dir, "{}_{}{}{}.png".format(id, seq[i], i+1, aa))
                img_3d = Image.open(png_file)
                
                mut_matrix_file = os.path.join(mut_matrices_dir, "{}_matrix_{}{}{}.png".format(id, seq[i], i+1, aa))
                img_mut_matrix = Image.open(mut_matrix_file)
                
                tar_img = Image.new('RGB', (movie_width, movie_height), (255,255,255))
                tar_img.paste(img_3d, (0,0))
                tar_img.paste(img_mut_matrix, (movie_width - matrix_frame_width, 0))
                tar_img.save(composite_file)

                counter += 1
                printProgressBar(counter, total, prefix = 'Composing final frames:', suffix = 'Complete', length = 50)
                
        end_time = time.time()
        print("Elapsed time: {}".format(end_time - start_time))

        print("Rendering movie for {}".format(id))
        start_time = time.time()
        os.system("ffmpeg -y -f image2 -framerate 19 -i {}/%d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p {}.mp4 > /dev/null 2>&1".format(composite_dir, id))
        end_time = time.time()
        print("Elapsed time: {}".format(end_time - start_time))

        # Cleanup
        print("Done, cleaning up..")
        shutil.rmtree(tmp_dir)
        print("")

if __name__ == "__main__":
    main()
