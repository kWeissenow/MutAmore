import os
import sys
import time
from PIL import Image, ImageFont, ImageDraw, ImageEnhance
from progressBar import *
from utils import get_script_path

if "pymol" not in  "\t".join(sys.path):
    sys.path.append("/usr/bin/pymol")
    sys.path.append("/usr/lib/python")
    import __main__
    __main__.pymol_argv = ['pymol','-Qqc']
    import pymol
    pymol.finish_launching()
else: # This IF avoids re-appending to PATH in case you execute this snippet multiple times in the same session
    import pymol

pymol.cmd.feedback("disable", "all", "actions")
pymol.cmd.feedback("disable", "all", "results")


def write_png(id, pdb_file, png_file, width=720, height=720, scale_factor=1.0, zoom_factor=None, transparency=False):
    if os.path.isfile(png_file):
        return
    pdb_name = os.path.basename(pdb_file).split('.')[0]
    parts = pdb_name.split("_")
    suffix = parts[-1]
    from_res = suffix[0]
    to_res = suffix[-1]
    mid = suffix[1:-1]
    try:
        residx = int(mid)
    except ValueError:
        residx = 0

    pymol.cmd.load(pdb_file, pdb_name)
    pymol.cmd.align(pdb_name, "wt")
    pymol.cmd.reset()
    pymol.cmd.disable("all")
    pymol.cmd.enable(pdb_name)
    pymol.cmd.hide('all')
    if zoom_factor is not None:
        pymol.cmd.zoom("center", zoom_factor)
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 1)
    pymol.cmd.set('depth_cue', 0)
    pymol.cmd.bg_color('white')
    pymol.cmd.spectrum('b', palette='red red red orange yellow cyan blue', minimum=0, maximum=100)
    pymol.cmd.color('black', 'resi {}'.format(residx))
    pymol.cmd.png(str(png_file), width=width, height=height, quiet=1)
    pymol.cmd.delete(pdb_name)

    img = Image.open(png_file)
    if transparency:
        contrast_enhancer = ImageEnhance.Contrast(img)
        img = contrast_enhancer.enhance(0.6)
    draw = ImageDraw.Draw(img)
    font_file = os.path.join(get_script_path(), "font.ttf")
    font = ImageFont.truetype(font_file, int(16*scale_factor))
    draw.text((1, 1), suffix, (0,0,0), font=font)
    img.save(png_file)
    
    
def render_3d_frames(id, seq, pdb_dir, png_dir, width, height, scale_factor, zoom_factor=None, experimental_mutations=None, experimental_dir=None, topN_indices=None):
    start_time = time.time()

    wt_file = os.path.join(pdb_dir, "{}.pdb".format(id))
    pymol.cmd.load(wt_file, "wt")

    total = len(seq) * 19
    printProgressBar(0, total, prefix = 'Rendering 3D structures:', suffix = 'Complete', length = 50)
    counter = 0
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")
    for i in range(len(seq)):
        for aa in list("ACDEFGHIKLMNPQRSTVWY"):
            if aa == seq[i]:
                continue
            if topN_indices is not None and topN_indices[aa_list.index(aa), i] == False:
                counter += 1
                printProgressBar(counter, total, prefix = 'Rendering 3D structures:', suffix = 'Complete', length = 50)
                continue
            temp_seq = seq.copy()
            temp_seq[i] = aa
            
            mut_name = seq[i] + str(i+1) + aa
            transparency = False
            if experimental_mutations is not None and not mut_name in experimental_mutations:
                transparency = True
            
            if experimental_mutations is not None and mut_name in experimental_mutations:
                pdb_file = os.path.join(experimental_dir, mut_name + ".pdb")
            else:
                pdb_file = os.path.join(pdb_dir, "{}_{}{}{}.pdb".format(id, seq[i], i+1, aa))
            png_file = os.path.join(png_dir, "{}_{}{}{}.png".format(id, seq[i], i+1, aa))
            write_png(id, pdb_file, png_file, width=width, height=height, scale_factor=scale_factor, zoom_factor=zoom_factor, transparency=transparency)
            
            counter += 1
            printProgressBar(counter, total, prefix = 'Rendering 3D structures:', suffix = 'Complete', length = 50)
            
    pymol.cmd.reinitialize()
            
    end_time = time.time()
    print("Elapsed time: {}".format(end_time - start_time))
