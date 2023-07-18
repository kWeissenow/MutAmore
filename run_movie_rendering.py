import argparse
import os
import time
import shutil
import stat
from utils import *


def main():
    parser = argparse.ArgumentParser(
                    prog='MutaMoRe',
                    description='Rendering protein mutation movies from predicted 3D structures')

    parser.add_argument('-i', '--input_fasta', type=str, help="Input FASTA file", required=True)
    parser.add_argument('-o', '--output_dir', type=str, help="Output directory for movies (default: current working directory)", default="./")
    parser.add_argument('--width', type=int, help="Movie horizontal resolution (default: 1280)", default=1280)
    parser.add_argument('--height', type=int, help="Movie vertical resolution (default: 720)", default=720)
    parser.add_argument('-t', '--temp_dir', type=str, help="Directory for temporary files (default: ./tmp)", default="./tmp")
    parser.add_argument('-p', '--prediction_dir', type=str, help="Directory for predicted structures (default: ./tmp/predictions)", default="./tmp/predictions")
    parser.add_argument('-s', '--predictor_script', type=str, dest="template_script", help="Structure prediction script")
    parser.add_argument('--only-predict', dest="only_predict", action="store_true", help="Only run structure prediction")
    parser.add_argument('--only-render', dest="only_render", action="store_true", help="Only run movie rendering")
    args = parser.parse_args()

    input_fasta = args.input_fasta
    output_dir = args.output_dir
    movie_width = args.width
    movie_height = args.height
    tmp_dir = args.temp_dir
    prediction_dir = args.prediction_dir
    matrix_frame_width = 250
    matrix_margin_horizontal = 25
    matrix_margin_vertical = 20
    
    if args.only_predict and args.only_render:
        print("Please only set either --only-predict or --only-render")
    doPredict = True
    doRender = True
    if args.only_predict:
        doRender = False
    elif args.only_render:
        doPredict = False

    if doPredict and args.template_script is None:
        print("Please select a structure prediction script from the directory 'predictor_scripts/', e.g. using")
        print("-s predictor_scripts/esmfold.sh")
        quit()

    # Check if mutation matrices can be rendered in the selected movie resolution
    check_movie_resolution(input_fasta, movie_height, matrix_margin_vertical)
    
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)


    # Structure prediction
    
    if doPredict:
        if not os.path.isdir(prediction_dir):
            os.makedirs(prediction_dir)

        # Create mutated sequences
        mutated_fasta = os.path.join(tmp_dir, "mutated_sequences.fasta")
        with open(mutated_fasta, "w") as f:
            for record in SeqIO.parse(input_fasta, "fasta"):
                id = record.id
                seq = list(record.seq)

                f.write(">{}\n{}\n".format(id, "".join(seq)))
                for i in range(len(seq)):
                    for aa in list("ACDEFGHIKLMNPQRSTVWY"):
                        if aa == seq[i]:
                            continue
                        temp_seq = seq.copy()
                        temp_seq[i] = aa

                        f.write(">{}_{}{}{}\n{}\n".format(id, seq[i], i+1, aa, "".join(temp_seq)))

        # Prepare predictor script
        predictor_script = os.path.join(tmp_dir, "predictor.sh")
        prepare_predictor_script(args.template_script, predictor_script, mutated_fasta, prediction_dir)

        # Make script executable
        st = os.stat(predictor_script)
        os.chmod(predictor_script, st.st_mode | stat.S_IEXEC)

        # Call structure predictor
        print("Starting structure predictions")
        os.system("bash -i " + predictor_script)


    # Movie rendering

    if doRender:
        from render_3d_frames import render_3d_frames
        from render_mutation_matrices import render_mutation_matrices
        from compose_frames import compose_frames

        for record in SeqIO.parse(input_fasta, "fasta"):
            id = record.id
            seq = list(record.seq)

            print("Processing {}".format(id))

            # Prepare directory structure

            current_tmp_dir = os.path.join(tmp_dir, id)
            if not os.path.isdir(current_tmp_dir):
                os.makedirs(current_tmp_dir)

            png_dir = os.path.join(current_tmp_dir, "png")
            mut_matrices_dir = os.path.join(current_tmp_dir, "mut_matrices_png")
            composite_dir = os.path.join(current_tmp_dir, "composite_png")
            if not os.path.isdir(png_dir):
                os.makedirs(png_dir)
            if not os.path.isdir(mut_matrices_dir):
                os.makedirs(mut_matrices_dir)
            if not os.path.isdir(composite_dir):
                os.makedirs(composite_dir)

            # Render 3D frames and mutation matrices

            render_3d_frames(id, seq, prediction_dir, png_dir, movie_width - matrix_frame_width, movie_height)
            render_mutation_matrices(id, seq, movie_height, prediction_dir, mut_matrices_dir, width=matrix_frame_width, margin_horiz=matrix_margin_horizontal, margin_vert=matrix_margin_vertical)

            # Compose final frames

            compose_frames(id, seq, movie_width, movie_height, matrix_frame_width, mut_matrices_dir, png_dir, composite_dir)

            # Render output file

            print("Rendering movie for {}".format(id))
            start_time = time.time()
            out_file = os.path.join(output_dir, "{}.mp4".format(id))
            os.system("ffmpeg -y -f image2 -framerate 19 -i {}/%d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p {} > /dev/null 2>&1".format(composite_dir, out_file))
            end_time = time.time()
            print("Elapsed time: {}".format(end_time - start_time))

    # Cleanup
    print("Done, cleaning up..")
    shutil.rmtree(tmp_dir)
    print("")


if __name__ == "__main__":
    main()
