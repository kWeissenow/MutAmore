from Bio import SeqIO
import os
import sys


def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


def check_movie_resolution(input_fasta, movie_height, matrix_margin_vertical):
    print("Checking protein lengths in input file..")
    total_min_height = 0
    for record in SeqIO.parse(input_fasta, "fasta"):
        id = record.id
        seq = list(record.seq)
        min_height = len(seq) + 2 * matrix_margin_vertical
        if min_height > movie_height:
            print("You will need a vertical resolution of at least {} to render the protein with identifier {}".format(min_height, id))
            if min_height > total_min_height:
                total_min_height = min_height
    if total_min_height != 0:
        print("")
        print("To render all sequences in your input file, you would need a vertical resolution of at least {}.".format(total_min_height))
        if total_min_height <= 1080:
            print("Consider rendering in Full-HD by passing the command line arguments: --width 1920 --height 1080")
        elif total_min_height <= 2160:
            print("Consider rendering in 4K resolution by passing the command line arguments: --width 3840 --height 2160")
            print("This will take a while and use a lot of disk space. Instead, consider removing the longest sequences from your input file.")
        else:
            print("You seem to have particularly long proteins in your input file! You can try setting a movie resolution as indicated above, but rendering times and disk space requirements will be massive.")
            print("Consider removing particulary long sequences from your input file.")
        quit()


def prepare_predictor_script(template_script, output_script, input_file, output_dir):
    with open(template_script, "r") as f_in:
        with open(output_script, "w") as f_out:
            for line in f_in:
                f_out.write(line.replace("MUTAMORE_INPUT", input_file).replace("MUTAMORE_OUTPUT", output_dir))
