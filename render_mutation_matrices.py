from Bio.PDB import *
import numpy as np
import os
import time
from scipy.spatial.distance import cdist
from PIL import Image, ImageFont, ImageDraw
from progressBar import *
from utils import get_script_path


def distance_map_from_pdb(parser, filename, seqlen):
    structure = parser.get_structure(filename, filename)

    c_alpha = np.zeros((seqlen, 3))
    for chain in structure[0]:
        if len(chain) != seqlen:
            return None
        for i,residue in enumerate(chain):
            c_alpha[i] = residue["CA"].get_coord()
        break

    return cdist(c_alpha, c_alpha)


def distance_diff(d1, d2):
    length = d1.shape[-1]
    cutoff = 15.0
    mask_2d = (d1 < cutoff) * (1.0 - np.eye(length))

    dist_l1 = np.absolute(d1 - d2)

    lddt_score = (dist_l1 < 0.5).astype(dist_l1.dtype) + \
                 (dist_l1 < 1.0).astype(dist_l1.dtype) + \
                 (dist_l1 < 2.0).astype(dist_l1.dtype) + \
                 (dist_l1 < 4.0).astype(dist_l1.dtype)
    lddt_score *= 0.25
    lddt_score *= mask_2d

    lddt_score = (1e-10 + np.sum(lddt_score, axis=(-1,))) / \
                 (1e-10 + np.sum(mask_2d, axis=(-1,)))
    return np.average(lddt_score)


def get_mutation_matrix(id, seq, pdb_dir, experimental_mutations, experimental_dir):
    parser = PDBParser(QUIET=True)
    wt_file = os.path.join(pdb_dir, "{}.pdb".format(id))
    ref_distance_map = distance_map_from_pdb(parser, wt_file, len(seq))
    
    start_time = time.time()
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")
    total = len(seq) * 19
    counter = 0
    printProgressBar(0, total, prefix = 'Calculating structural similarity:', suffix = 'Complete', length = 50)
    mut_matrix = np.ones((20, len(seq)))
    for i in range(len(seq)):
        for aa in aa_list:
            if aa == seq[i]:
                continue
            temp_seq = seq.copy()
            temp_seq[i] = aa    
            mut_name = seq[i] + str(i+1) + aa

            if experimental_mutations is not None and mut_name in experimental_mutations:
                filename = os.path.join(experimental_dir, mut_name + ".pdb")
            else:
                filename = os.path.join(pdb_dir, "{}_{}{}{}.pdb".format(id, seq[i], i+1, aa))
            dist_map = distance_map_from_pdb(parser, filename, len(seq))
            
            if dist_map is None:
                print("WARNING: ignoring experimental structure {} since it does not have the same length as the wild-type structure".format(filename))
                mut_matrix[aa_list.index(aa), i] = 0.0
                continue

            delta = distance_diff(ref_distance_map, dist_map)
            mut_matrix[aa_list.index(aa), i] = delta
            
            counter += 1
            printProgressBar(counter, total, prefix = 'Calculating structural similarity:', suffix = 'Complete', length = 50)
            
    end_time = time.time()
    print("Elapsed time: {}".format(end_time - start_time))
            
    return mut_matrix


def gradient_color(minval, maxval, val, color_palette=((0,0,0), (255,0,0), (255, 165, 0), (255,255,255))):
    """ Computes intermediate RGB color of a value in the range of minval
        to maxval (inclusive) based on a color_palette representing the range.
    """
    max_index = len(color_palette)-1
    delta = maxval - minval
    if delta == 0:
        delta = 1
    v = float(val-minval) / delta * max_index
    i1, i2 = int(v), min(int(v)+1, max_index)
    (r1, g1, b1), (r2, g2, b2) = color_palette[i1], color_palette[i2]
    f = v - i1
    return int(r1 + f*(r2-r1)), int(g1 + f*(g2-g1)), int(b1 + f*(b2-b1))
    
    
def render_matrix_frames(id, seq, mut_matrix, legend, aa_labels, height, out_dir, width, margin_horiz, margin_vert, scale_factor, topN_indices):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    start_time = time.time()

    length = len(seq)

    ticks = 10
    if length > 400:
        ticks = 50
    
    font_size = int(12 * scale_factor)
    font_file = os.path.join(get_script_path(), "font.ttf")
    font = ImageFont.truetype(font_file, font_size)
    
    cell_width = int(5 * scale_factor)
    cell_height = int((height - margin_vert*2) / length)
    
    # calculate offset to center the mutation matrix vertically
    matrix_height = length*cell_height + 1
    offset_y = int(height / 2 - matrix_height / 2)
    
    if cell_height == 0:
        print("ERROR: Can not render mutation matrix. Your protein is too long for the requested movie resolution.")
        print("Try setting a higher (vertical) resolution.")
        quit()

    aa_list = list("ACDEFGHIKLMNPQRSTVWY")
    total = len(seq) * 19
    counter = 0
    printProgressBar(0, total, prefix = 'Rendering mutation matrices:', suffix = 'Complete', length = 50)
    for i in range(len(seq)):
        for aa in aa_list:
            if aa == seq[i]:
                continue
            if topN_indices is not None and topN_indices[aa_list.index(aa), i] == False:
                counter += 1
                printProgressBar(counter, total, prefix = 'Rendering mutation matrices:', suffix = 'Complete', length = 50)
                continue
            temp_seq = seq.copy()
            temp_seq[i] = aa    
            
            im = Image.new('RGBA', (width, height), (255, 255, 255, 255))
            draw = ImageDraw.Draw(im)
            
            for idx in range(length):
                if idx+1 == 1 or (idx+1) % ticks == 0:
                    x0 = margin_horiz + 20*cell_width + 1
                    y0 = offset_y + idx * cell_height - font_size/2 - 1
                    draw.text((x0, y0), "-{}".format(idx+1), (0,0,0), font=font)
        
                for m in range(20):
                    x0 = margin_horiz + m * cell_width
                    y0 = offset_y + idx * cell_height
                    x1 = x0 + cell_width
                    y1 = y0 + cell_height
                    color = gradient_color(0.3, 1.0, mut_matrix[m][idx])
                    outline = None
                    if idx == i and m == aa_list.index(aa):
                        outline = (0,0,0)
                        x1 -= 1
                        y1 -= 1
                    draw.rectangle((x0,y0,x1,y1), fill=color, outline=outline)
            
            # outline for the whole mutation matrix
            x0 = margin_horiz - 1
            y0 = offset_y - 1
            x1 = x0 + 20*cell_width + 1
            y1 = y0 + length*cell_height + 1
            draw.rectangle((x0,y0,x1,y1), outline=(0,0,0))
            
            # amino acid labels
            im.paste(aa_labels, (margin_horiz, offset_y - int(11*scale_factor)))
    
            # left boundary line
            draw.line([(0,0), (0,height)], fill=(100,100,100))
    
            # draw legend
            legend_x_offset = int(60 * scale_factor)
            legend_y_offset = int(30 * scale_factor)
            im.paste(legend, (width - margin_horiz - legend_x_offset, int(height / 2) - legend_y_offset))
            
            filename = os.path.join(out_dir, "{}_matrix_{}{}{}.png".format(id, seq[i], i+1, aa))
            im.save(filename)
            counter += 1
            printProgressBar(counter, total, prefix = 'Rendering mutation matrices:', suffix = 'Complete', length = 50)
            
    end_time = time.time()
    print("Elapsed time: {}".format(end_time - start_time))
    
    
def draw_legend(scale_factor):
    font_size = int(12 * scale_factor)

    font_file = os.path.join(get_script_path(), "font.ttf")
    font = ImageFont.truetype(font_file, font_size)

    size_x = int(60 * scale_factor)
    size_y = int(110 * scale_factor)
    im = Image.new('RGBA', (size_x, size_y), (255, 255, 255, 255))
    draw = ImageDraw.Draw(im)
    
    legend_width = int(10 * scale_factor)
    legend_height = int(100 * scale_factor)
    legend_y_offset = int(5 * scale_factor)
    for i in range(legend_height):
        x0 = 1
        y0 = 1 + i + legend_y_offset
        x1 = x0 + legend_width
        y1 = y0
        color = gradient_color(0.0, 1.0, (legend_height-i) / legend_height)
        draw.rectangle((x0,y0,x1,y1), fill=color)
    
    # outline    
    x0 = 0
    y0 = legend_y_offset
    x1 = x0 + legend_width + 1
    y1 = y0 + legend_height
    draw.rectangle((x0,y0,x1,y1), outline=(0,0,0))
    
    # draw axis ticks    
    x0 = legend_width + 1
    y0 = legend_y_offset - font_size/2 - 1
    draw.text((x0, y0), "100%", (0,0,0), font=font)
    draw.text((x0, y0+legend_height), "30%", (0,0,0), font=font)
    
    # draw axis annotation
    size_x = int(100 * scale_factor)
    size_y = int(15 * scale_factor)
    im_annotation = Image.new('RGBA', (size_x, size_y), (255, 255, 255, 255))
    draw_annotation = ImageDraw.Draw(im_annotation)
    draw_annotation.text((0,0), "structure similarity", (0,0,0), font=font)
    im_annotation = im_annotation.transpose(method=Image.Transpose.ROTATE_90)
    pos_x = int(45 * scale_factor)
    pos_y = int(5 * scale_factor)
    im.paste(im_annotation, (pos_x, pos_y))
    
    return im
    

def draw_amino_acid_labels(scale_factor):
    font_size = int(7 * scale_factor)
    font_file = os.path.join(get_script_path(), "font.ttf")
    font = ImageFont.truetype(font_file, font_size)
    
    size_x = int(100 * scale_factor)
    size_y = int(10 * scale_factor)
    im = Image.new('RGBA', (size_x, size_y), (255, 255, 255, 255))
    draw = ImageDraw.Draw(im)
    for i,c in enumerate("ACDEFGHIKLMNPQRSTVWY"):
        pos = int(i * 5 * scale_factor)
        draw.text((pos,0), c, (0,0,0), font=font)
    
    return im
    

def render_mutation_matrices(id, seq, height, pdb_dir, out_dir, width=250, margin_horiz=25, margin_vert=20, scale_factor=1.0, experimental_mutations=None, experimental_dir=None, topN=None):
    # draw legend
    legend = draw_legend(scale_factor)
    
    # draw amino acid labels
    aa_labels = draw_amino_acid_labels(scale_factor)

    mut_matrix = get_mutation_matrix(id, seq, pdb_dir, experimental_mutations, experimental_dir)
    
    # if topN is set, get positions with highest structural difference
    topN_indices = None
    if topN is not None:
        topN_similarity = bottom_k = np.partition(mut_matrix.flatten(), topN)[:topN]
        threshold = np.sort(topN_similarity)[-1]
        topN_indices = mut_matrix <= threshold
    
    render_matrix_frames(id, seq, mut_matrix, legend, aa_labels, height, out_dir, width=width, margin_horiz=margin_horiz, margin_vert=margin_vert, scale_factor=scale_factor, topN_indices=topN_indices)
    return topN_indices
