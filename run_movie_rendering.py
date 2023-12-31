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
    parser.add_argument('--experimental_dir', type=str, help="Directory for experimental structures (optional)")
    parser.add_argument('-s', '--predictor_script', type=str, dest="template_script", help="Structure prediction script")
    parser.add_argument('-z', '--zoom_factor', type=float, help="Specific zoom level to be used for 3D rendering (optional)")
    parser.add_argument('--top', type=int, help="Only show top-N mutants with structural difference")
    parser.add_argument('--only-predict', dest="only_predict", action="store_true", help="Only run structure prediction")
    parser.add_argument('--only-render', dest="only_render", action="store_true", help="Only run movie rendering")
    args = parser.parse_args()

    input_fasta = args.input_fasta
    output_dir = args.output_dir
    movie_width = args.width
    movie_height = args.height
    tmp_dir = args.temp_dir
    prediction_dir = args.prediction_dir
    
    scale_factor = movie_height / 360

    matrix_frame_width = int(200 * scale_factor)
    matrix_margin_horizontal = int(5 * scale_factor)
    matrix_margin_vertical = int(10 * scale_factor)
    
    if args.only_predict and args.only_render:
        print("Please only set either --only-predict or --only-render")
    doPredict = True
    doRender = True
    if args.only_predict:
        doRender = False
    elif args.only_render:
        doPredict = False
    zoom_factor = None
    if args.zoom_factor:
        zoom_factor = args.zoom_factor
        
    topN = None
    framerate = 19
    if args.top:
        topN = args.top
        framerate = topN / 5
        if framerate < 2:
            framerate = 2
        if framerate > 19:
            framerate = 19

    if doPredict and args.template_script is None:
        print("Please select a structure prediction script from the directory 'predictor_scripts/', e.g. using")
        print("-s predictor_scripts/esmfold.sh")
        quit()

    # Check if mutation matrices can be rendered in the selected movie resolution
    check_movie_resolution(input_fasta, movie_height, matrix_margin_vertical)
    
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)

    # parse optional experimental structures
    experimental_mutations = None
    if args.experimental_dir:
        print("Parsing experimental structures..")
        for file_name in os.listdir(args.experimental_dir):
            if file_name.endswith(".pdb"):
                mut_name = file_name[:-4]
                if experimental_mutations is None:
                    experimental_mutations = []
                experimental_mutations.append(mut_name)

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
                        
                        # skip experimental structures
                        mut_name = seq[i] + str(i+1) + aa
                        if experimental_mutations is not None and mut_name in experimental_mutations:
                            continue

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

            topN_indices = render_mutation_matrices(id, seq, movie_height, prediction_dir, mut_matrices_dir, width=matrix_frame_width, margin_horiz=matrix_margin_horizontal, margin_vert=matrix_margin_vertical, scale_factor=scale_factor, experimental_mutations=experimental_mutations, experimental_dir=args.experimental_dir, topN=topN)
            render_3d_frames(id, seq, prediction_dir, png_dir, movie_width - matrix_frame_width, movie_height, scale_factor, zoom_factor, experimental_mutations, args.experimental_dir, topN_indices)

            # Compose final frames

            compose_frames(id, seq, movie_width, movie_height, matrix_frame_width, mut_matrices_dir, png_dir, composite_dir, topN_indices)

            # Render output file

            print("Rendering movie for {}".format(id))
            start_time = time.time()
            out_file = os.path.join(output_dir, "{}.mp4".format(id))
            os.system("ffmpeg -y -f image2 -framerate {} -i {}/%d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p {} > /dev/null 2>&1".format(framerate, composite_dir, out_file))
            end_time = time.time()
            print("Elapsed time: {}".format(end_time - start_time))

    # Cleanup
    print("Done, cleaning up..")
    shutil.rmtree(tmp_dir)
    print("")


if __name__ == "__main__":
    main()
