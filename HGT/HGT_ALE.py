import os
from Bio import SeqIO
from ete3 import Tree
from PIL import Image
from itolapi import Itol
import multiprocessing as mp
from PyPDF3.pdf import PageObject
from PyPDF3 import PdfFileWriter, PdfFileReader
from ete3 import TextFace, TreeStyle, NodeStyle


def select_seq(arg_list):

    seq_file    = arg_list[0]
    id_file     = arg_list[1]
    output_file = arg_list[2]

    seq_id_set = {i.strip() for i in open(id_file)}
    output_file_handle = open(output_file, 'w')
    for seq_record in SeqIO.parse(seq_file, 'fasta'):
        if seq_record.id in seq_id_set:
            SeqIO.write(seq_record, output_file_handle, 'fasta-2line')
    output_file_handle.close()


def root_with_out_group(tree_file, out_group_txt, tree_file_rooted):

    out_group_set = set()
    for each_og in open(out_group_txt):
        out_group_set.add(each_og.strip())

    tre = Tree(tree_file, format=1)
    out_group_lca = tre.get_common_ancestor(out_group_set)
    tre.set_outgroup(out_group_lca)
    tre.write(outfile=tree_file_rooted)


def subset_tree(tree_file_in, leaves_to_keep_list, tree_file_out):

    input_tree = Tree(tree_file_in)
    subset_tree = input_tree.copy()
    subset_tree.prune(leaves_to_keep_list, preserve_branch_length=True)
    if tree_file_out is None:
        return subset_tree.write()
    else:
        subset_tree.write(outfile=tree_file_out)


def plot_tree(input_tree, tree_title, node_label_dict, node_label_color_dict, align_leaf_label, show_scale, output_plot):

    if os.path.isfile(input_tree) is False:
        print('Tree file not found, program exited!')
        print(input_tree)
        exit()

    t = Tree(input_tree)
    ts = TreeStyle()
    ts.mode = "r"                       # tree model: 'r' for rectangular, 'c' for circular
    ts.show_border = False              # set tree image border
    ts.show_leaf_name = False           # show/hide leaf name, hide here, so you can customise it below with node.add_face()
    ts.title.add_face(TextFace(tree_title, fsize=9, fgcolor='black', ftype='Arial', tight_text=False), column=0)  # add tree title

    # set node style
    for each_node in t.traverse():
        ns = NodeStyle()
        ns["shape"]         = "circle"  # dot shape: circle, square or sphere
        ns["fgcolor"]       = "black"   # color of shape(not label)
        ns['size']          = 0         # node shape size
        ns['hz_line_type']  = 0         # horizontal branch line type: 0 for solid, 1 for dashed, 2 for dotted
        ns['vt_line_type']  = 0         # vertical branch line type:   0 for solid, 1 for dashed, 2 for dotted
        ns['hz_line_width'] = 0.5       # horizontal branch line width
        ns['vt_line_width'] = 0.5       # vertical branch line width

        leaf_label_position = 'branch-right'
        if align_leaf_label is True:
            leaf_label_position = 'aligned'

        if each_node.is_leaf():
            node_id = each_node.name
            node_label_color = node_label_color_dict.get(node_id, 'black')
            node_label_text  = node_label_dict.get(node_id, node_id)
            each_node.add_face(TextFace(node_label_text, fsize=8, fgcolor=node_label_color, tight_text=False, bold=False),
                               column=0, position=leaf_label_position)  # aligned, branch-right
        else:
            pass
        each_node.set_style(ns)

    # set layout
    ts.rotation                 = 0             # from 0 to 360
    ts.margin_top               = 10            # top tree image margin
    ts.margin_bottom            = 10            # bottom tree image margin
    ts.margin_left              = 10            # left tree image margin
    ts.margin_right             = 10            # right tree image margin
    ts.branch_vertical_margin   = 3             # 3 pixels between adjancent branches
    ts.show_scale               = show_scale    # show_scale
    ts.show_border              = False         # set tree image border

    # write out tree
    t.render(output_plot, w=1200, units="px", tree_style=ts)


def merge_image(image_file_list, output_image):

    images = [Image.open(x) for x in image_file_list]
    widths, heights = zip(*(i.size for i in images))
    total_width = sum(widths)
    max_height = max(heights)
    new_im = Image.new('RGB', (total_width, max_height), color='white')

    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset, 0))
        x_offset += im.size[0]

    new_im.save(output_image)


def merge_pdf(pdf_1, pdf_2, op_pdf):

    page1 = PdfFileReader(open(pdf_1, "rb"), strict=False).getPage(0)
    page2 = PdfFileReader(open(pdf_2, "rb"), strict=False).getPage(0)

    total_width = page1.mediaBox.upperRight[0] + page2.mediaBox.upperRight[0]
    total_height = max([page1.mediaBox.upperRight[1], page2.mediaBox.upperRight[1]])

    new_page = PageObject.createBlankPage(None, total_width, total_height)

    # Add first page at the 0,0 position
    new_page.mergePage(page1)
    # Add second page with moving along the axis x
    new_page.mergeTranslatedPage(page2, page1.mediaBox.upperRight[0], 0)

    output = PdfFileWriter()
    output.addPage(new_page)
    output.write(open(op_pdf, "wb"))


def uts_to_itol_connections(genome_tree_file, ale_formatted_gnm_tree, interal_node_prefix, uts_file, freq_cutoff, ignore_leaf_hgt, ignore_vertical_hgt, min_donor_node_leaf_num, min_recipient_node_leaf_num, itol_connection_txt):

    # get internal_node_to_leaf_dict
    internal_node_to_leaf_dict = get_node_to_leaf_dict(ale_formatted_gnm_tree)

    paired_donor_to_recipient_leaf_dict = dict()
    qualified_hgt_num = 0

    leaf_id_set = []
    if os.path.isfile(genome_tree_file):
        leaf_id_set = [i.name for i in Tree(genome_tree_file, format=3).get_leaves()]
    else:
        print('%s not found!' % genome_tree_file)

    with open(itol_connection_txt, 'w') as itol_connection_txt_handle:
        itol_connection_txt_handle.write('DATASET_CONNECTION\nSEPARATOR TAB\nDATASET_LABEL\tdemo_connections\n')
        itol_connection_txt_handle.write('COLOR\t#ff0ff0\nDRAW_ARROWS\t1\nARROW_SIZE\t20\nLOOP_SIZE\t100\n')
        itol_connection_txt_handle.write('MAXIMUM_LINE_WIDTH\t10\nCURVE_ANGLE\t45\nCENTER_CURVES\t1\nALIGN_TO_LABELS\t0\nDATA\n')
        for each_line in open(uts_file):
            if not each_line.startswith('#'):
                each_line_split = each_line.strip().split('\t')
                donor = each_line_split[0]
                recipient = each_line_split[1]

                # add prefix to internal donor node
                if donor in leaf_id_set:
                    donor_with_prefix = donor
                else:
                    donor_with_prefix = interal_node_prefix + donor

                # add prefix to internal recipient node
                if recipient in leaf_id_set:
                    recipient_with_prefix = recipient
                else:
                    recipient_with_prefix = interal_node_prefix + recipient

                freq = float(each_line_split[2])
                if freq >= freq_cutoff:
                    if ignore_leaf_hgt is False:
                        if ignore_vertical_hgt is False:
                            itol_connection_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s->%s(%s)\n' % (donor_with_prefix, recipient_with_prefix, freq, '#EB984E', 'normal', donor_with_prefix, recipient_with_prefix, freq))
                            qualified_hgt_num += 1
                        else:
                            donor_is_ancestor_of_recipient = check_a_is_ancestor_of_b(ale_formatted_gnm_tree, donor, recipient)
                            donor_is_child_of_recipient    = check_a_is_child_of_b(ale_formatted_gnm_tree, donor, recipient)
                            if (donor_is_ancestor_of_recipient is False) and (donor_is_child_of_recipient is False):
                                itol_connection_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s->%s(%s)\n' % (donor_with_prefix, recipient_with_prefix, freq, '#EB984E', 'normal', donor_with_prefix, recipient_with_prefix, freq))
                                qualified_hgt_num += 1
                    else:
                        if (each_line_split[0] not in leaf_id_set) and (each_line_split[1] not in leaf_id_set):
                            donor_node_leaf_num = len(internal_node_to_leaf_dict.get(donor, []))
                            recipient_node_leaf_num = len(internal_node_to_leaf_dict.get(recipient, []))
                            if (donor_node_leaf_num >= donor_node_min_leaf_num) and (recipient_node_leaf_num >= recipient_node_min_leaf_num):
                                if ignore_vertical_hgt is False:
                                    itol_connection_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s->%s(%s)\n' % (donor_with_prefix, recipient_with_prefix, freq, '#EB984E', 'normal', donor_with_prefix, recipient_with_prefix, freq))
                                    qualified_hgt_num += 1
                                else:
                                    donor_is_ancestor_of_recipient = check_a_is_ancestor_of_b(ale_formatted_gnm_tree, donor, recipient)
                                    donor_is_child_of_recipient    = check_a_is_child_of_b(ale_formatted_gnm_tree, donor, recipient)
                                    if (donor_is_ancestor_of_recipient is False) and (donor_is_child_of_recipient is False):
                                        itol_connection_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s->%s(%s)\n' % (donor_with_prefix, recipient_with_prefix, freq, '#EB984E', 'normal', donor_with_prefix, recipient_with_prefix, freq))
                                        qualified_hgt_num += 1
                                        key_str = '%s___%s' % (donor_with_prefix, recipient_with_prefix)
                                        paired_donor_to_recipient_leaf_dict[key_str] = [internal_node_to_leaf_dict.get(donor, []), internal_node_to_leaf_dict.get(recipient, [])]

    return qualified_hgt_num, internal_node_to_leaf_dict, paired_donor_to_recipient_leaf_dict


def itol_tree(tree_file, annotation_file_list, project_name, APIkey, display_mode, op_plot):

    # https://github.com/albertyw/itolapi
    # http://itol.embl.de/help.cgi#batch

    op_plot_ext = op_plot.split('.')[-1]

    # upload tree to iTOL
    itol_uploader = Itol()
    itol_uploader.params['projectName'] = project_name  # better to create a project with a unique name.
    itol_uploader.params['APIkey'] = APIkey  # sine we are the same account, we can use the same APIkey
    itol_uploader.params['treeName'] = tree_file
    itol_uploader.add_file(tree_file)

    # upload annotation files to iTOL
    for annotation_file in annotation_file_list:
        itol_uploader.add_file(annotation_file)

    status = itol_uploader.upload()
    # import pdb;pdb.set_trace()
    assert status != False

    # the following parameters are optional, refer to https://itol.embl.de/help.cgi#batchExp
    if len(annotation_file_list) == 1:
        datasets_visible_str = '0'
    elif len(annotation_file_list) == 2:
        datasets_visible_str = '0,1'
    elif len(annotation_file_list) == 3:
        datasets_visible_str = '0,1,2'
    else:
        datasets_visible_str = ','.join([str(i) for i in list(range(0, len(annotation_file_list)))])
    itol_exporter = itol_uploader.get_itol_export()
    itol_exporter.set_export_param_value('datasets_visible', datasets_visible_str)
    itol_exporter.set_export_param_value('display_mode', display_mode)
    itol_exporter.set_export_param_value('range_mode', '2')
    itol_exporter.set_export_param_value('dashed_lines', '0')
    # itol_exporter.set_export_param_value('current_font_size', '96')
    itol_exporter.set_export_param_value('line_width', '3')
    # itol_exporter.set_export_param_value('vertical_shift_factor', '3')
    # itol_exporter.set_export_param_value('horizontal_scale_factor', '3')
    itol_exporter.set_export_param_value('format', op_plot_ext)
    itol_exporter.export(op_plot)


def get_node_to_leaf_dict(tree_file):
    internal_node_to_leaf_dict = dict()
    for node in Tree(tree_file, format=1).traverse():
        if not node.is_leaf():
            node_name = node.name
            node_leaf_list = node.get_leaf_names()
            internal_node_to_leaf_dict[node_name] = node_leaf_list
    return internal_node_to_leaf_dict


def combine_trees(t1_with_len, t2_with_name, op_tree_with_both):

    # assume t1 has brancn length
    # assume t2 has internal node name

    t1 = Tree(t1_with_len, format=0)
    t2 = Tree(t2_with_name, format=1)

    t1_leaves_to_node_dict = dict()
    for t1_node in t1.traverse():
        leaf_str = '__'.join(sorted(list(t1_node.get_leaf_names())))
        t1_leaves_to_node_dict[leaf_str] = t1_node

    t2_leaves_to_node_dict = dict()
    for t2_node in t2.traverse():
        leaf_str = '__'.join(sorted(list(t2_node.get_leaf_names())))
        t2_leaves_to_node_dict[leaf_str] = t2_node

    t1_node_to_t2_node_dict = dict()
    for index, t1_node in t1_leaves_to_node_dict.items():
        t2_node = t2_leaves_to_node_dict[index]
        t1_node_to_t2_node_dict[t1_node] = t2_node

    merged_tree = t1.copy()
    for node, t1_node in zip(merged_tree.traverse(), t1.traverse()):
        node.name = t1_node_to_t2_node_dict[t1_node].name
    merged_tree.write(outfile=op_tree_with_both, format=3)


def prefix_internal_nodes(tree_in, prefix_str, tree_out):
    t = Tree(tree_in, format=3)
    t_renamed = t.copy()
    for node in t_renamed.traverse():
        if not node.is_leaf():
            node_name_prefixed = '%s%s' % (prefix_str, node.name)
            node.name = node_name_prefixed
        t_renamed.write(outfile=tree_out, format=3)


def check_a_is_ancestor_of_b(tree_file, node_a, node_b):

    a_is_ancestor_of_b = False
    for node in Tree(tree_file, format=1).traverse():
        node_name = node.name
        if node_name == node_b:
            node_ancestor_list = [i.name for i in node.get_ancestors()]
            if node_a in node_ancestor_list:
                a_is_ancestor_of_b = True

    return a_is_ancestor_of_b


def check_a_is_child_of_b(tree_file, node_a, node_b):

    a_is_child_of_b = False
    for node in Tree(tree_file, format=1).traverse():
        node_name = node.name
        if node_name == node_b:
            node_children_list = [i.name for i in node.get_descendants()]
            if node_a in node_children_list:
                a_is_child_of_b = True

    return a_is_child_of_b


########################################################################################################################

# file in
wd                              = '/Users/songweizhi/Desktop/DateArTree/0_HGT_MetaCHIP'
orthogroups_op_txt              = '%s/Orthogroups.txt'                                                  % wd
taxon_for_MetaCHIP_txt          = '%s/taxon_for_MetaCHIP.txt'                                           % wd
combined_faa                    = '%s/Archaea_133_HGT_pc_pc_combined_faa.fasta'                         % wd
genome_tree_file                = '%s/concatenated.treefile'                                            % wd
outgroup                        = '%s/out_group.txt'                                                    % wd
min_og_genome_num               = 50
min_og_phylum_num               = 2
align_leaf_name                 = True
show_scale                      = False
ignore_leaf_hgt                 = True
project_name                    = 'batch_access_tmp'
API_key                         = 'S1kZZuDHc0d5M7J5vLnUNQ'
display_mode                    = '1'  # 1=rectangular, 2=circular, 3=unrooted
hgt_freq_cutoff                 = 0.3
num_threads                     = 10
js_num_threads                  = 2
interal_node_prefix             = 'IN'
ale_splitter_py                 = '/home-user/wzsong/Tests/ALE/ALEtutorial/ale_splitter_modified.py'
ignore_vertical_hgt             = True
donor_node_min_leaf_num         = 5
recipient_node_min_leaf_num     = 5
d_color                         = '#FF0000'
r_color                         = '#0000FF'

# file out
op_dir                          = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/op_qualified_OGs'
gene_tree_dir                   = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/op_qualified_OGs_gene_tree_dir'
genome_tree_file_rooted         = '%s/concatenated_rooted.treefile'                                                 % op_dir
ale_wd                          = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_wd'
ale_op_dir                      = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_op_dir'
ale_hgt_plot_dir                = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_hgt_plot_dir'

########################################################################################################################

force_create_op_dir             = True
force_create_ale_wd             = True
extract_sequence                = False
prepare_ale_input_files         = False

# specify OGs to process
# designate_ogs                   = ['OG0000006', 'OG0000007', 'OG0000011', 'OG0000012', 'OG0000014', 'OG0000015', 'OG0000017']
designate_ogs                   = next(os.walk(gene_tree_dir))[1]
#designate_ogs                   = ['OG0000408']

'''
The number of orthogroups spanning >= 50 genomes and >= 2 phyla is 763.

OG0000294 is very interesting!!!
'''

########################################################################################################################

if force_create_op_dir is True:
    if os.path.isdir(op_dir) is True:
        os.system('rm -r %s' % op_dir)
os.system('mkdir %s' % op_dir)

# root genome tree with outgroup
root_with_out_group(genome_tree_file, outgroup, genome_tree_file_rooted)

# read in genome taxonomy
gnm_p_dict = dict()
gnm_c_dict = dict()
gnm_o_dict = dict()
gnm_pco_dict = dict()
for each_gnm in open(taxon_for_MetaCHIP_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id     = each_gnm_split[0]
    taxon_str  = each_gnm_split[1]
    gnm_phylum = taxon_str.split(';')[1]
    gnm_class  = taxon_str.split(';')[2]
    gnm_order  = taxon_str.split(';')[3]
    gnm_id_new = gnm_id
    if '.gtdb' in gnm_id_new:
        gnm_id_new = gnm_id_new.replace('.gtdb', '')
    gnm_p_dict[gnm_id_new] = gnm_phylum
    gnm_c_dict[gnm_id_new] = gnm_class
    gnm_o_dict[gnm_id_new] = gnm_order
    gnm_pco_dict[gnm_id_new] = '%s__%s__%s__%s' % (gnm_phylum[3:], gnm_class[3:], gnm_order[3:], gnm_id_new)

# read in OrthoFinder output
ortho_to_gene_dict = dict()
gene_to_ortho_dict = dict()
for each_og in open(orthogroups_op_txt):
    each_og_split = each_og.strip().split(' ')
    og_id = each_og_split[0][:-1]
    gene_list = each_og_split[1:]
    ortho_to_gene_dict[og_id] = gene_list
    for each_gene in gene_list:
        gene_to_ortho_dict[each_gene] = og_id

# get qualified orthogroups
qualified_og_set = set()
for each_ortho in ortho_to_gene_dict:
    ortho_gene_set = ortho_to_gene_dict[each_ortho]
    ortho_p_set = set()
    ortho_gnm_set = set()
    for each_gene in ortho_gene_set:
        gene_gnm = '_'.join(each_gene.split('_')[:-1]).replace('.gtdb', '')
        gnm_taxon = gnm_p_dict[gene_gnm]
        ortho_gnm_set.add(gene_gnm)
        ortho_p_set.add(gnm_taxon)
    if (len(ortho_gnm_set) >= min_og_genome_num) and (len(ortho_p_set) >= min_og_phylum_num):
        qualified_og_set.add(each_ortho)
print('The total number of identified orthogroups is %s.' % len(ortho_to_gene_dict))
print('The number of orthogroups spanning >= %s genomes and >= %s phyla is %s.' % (min_og_genome_num, min_og_phylum_num, len(qualified_og_set)))

# process qualified OG
og_to_process = sorted([i for i in qualified_og_set])
if len(designate_ogs) > 0:
    print('The number of designated OGs to process: %s' % len(designate_ogs))
    og_to_process = designate_ogs

# extract gene sequences and prepare commands for building gene tree
print('Extracting sequences and preparing commands for building gene trees')
extract_seq_arg_lol = []
for qualified_og in og_to_process:
    qualified_og_gene_set         = ortho_to_gene_dict[qualified_og]
    qualified_og_gene_txt         = '%s/%s.txt'           % (op_dir, qualified_og)
    qualified_og_gene_faa         = '%s/%s.faa'           % (op_dir, qualified_og)
    qualified_og_gene_aln         = '%s/%s.aln'           % (op_dir, qualified_og)
    qualified_og_gene_aln_trimmed = '%s/%s_trimmed.aln'   % (op_dir, qualified_og)
    js_file                       = '%s/zjs_%s.sh'        % (op_dir, qualified_og)

    # write out the id of genes
    with open(qualified_og_gene_txt, 'w') as qualified_og_gene_txt_handle:
        qualified_og_gene_txt_handle.write('\n'.join(qualified_og_gene_set))

    # add to mp lol
    extract_seq_arg_lol.append([combined_faa, qualified_og_gene_txt, qualified_og_gene_faa])

    # write out js for mafft, trimal and iqtree
    mafft_cmd          = 'mafft-einsi --thread %s --quiet %s > %s'                          % (js_num_threads, qualified_og_gene_faa, qualified_og_gene_aln)
    trimal_cmd         = 'trimal -in %s -out %s -automated1'                                % (qualified_og_gene_aln, qualified_og_gene_aln_trimmed)
    iqtree_cmd         = 'iqtree -m LG+G+I -bb 1000 --wbtl -nt %s -s %s -pre %s'            % (js_num_threads, qualified_og_gene_aln, qualified_og)
    iqtree_cmd_trimmed = 'iqtree -m LG+G+I -bb 1000 --wbtl -nt %s -s %s -pre %s_trimmed'    % (js_num_threads, qualified_og_gene_aln_trimmed, qualified_og)
    mv_file_cmd_1      = 'mkdir %s'                                                         % (qualified_og)
    mv_file_cmd_2      = 'mv %s.* %s/'                                                      % (qualified_og, qualified_og)
    js_file_handle = open(js_file, 'w')
    js_file_handle.write('#!/bin/bash\n#SBATCH --ntasks 1\n#SBATCH --cpus-per-task %s\n' % js_num_threads)
    js_file_handle.write(mafft_cmd.replace((op_dir + '/'), '') + '\n')
    # js_file_handle.write(trimal_cmd.replace((op_dir + '/'), '') + '\n')
    # js_file_handle.write(iqtree_cmd_trimmed.replace((op_dir + '/'), '') + '\n')
    js_file_handle.write(iqtree_cmd.replace((op_dir + '/'), '') + '\n')
    js_file_handle.write(mv_file_cmd_1 + '\n')
    js_file_handle.write(mv_file_cmd_2 + '\n')
    js_file_handle.close()

# extract gene sequences with multiprocessing
if extract_sequence is True:
    print('Extracting gene sequences with %s cores' % num_threads)
    pool = mp.Pool(processes=num_threads)
    pool.map(select_seq, extract_seq_arg_lol)
    pool.close()
    pool.join()

########################################################################################################################

# prepare input files and job script for running ALE
if prepare_ale_input_files is True:

    if force_create_ale_wd is True:
        if os.path.isdir(ale_wd) is True:
            os.system('rm -r %s' % ale_wd)
    os.system('mkdir %s' % ale_wd)
    n = 1
    for qualified_og in og_to_process:

        genome_tree_file_subset             = '%s_genome_tree.treefile'             % qualified_og
        genome_tree_file_subset_for_ale     = '%s_genome_tree_for_ALE.treefile'     % qualified_og
        gene_tree_ufboot                    = '%s.ufboot'                           % qualified_og
        gene_tree_ufboot_for_ale            = '%s_for_ALE.ufboot'                   % qualified_og
        gene_tree_treefile                  = '%s.treefile'                         % qualified_og
        gene_tree_treefile_subset           = '%s_subset.treefile'                  % qualified_og
        js_ale                              = '/zjs_%s_ALE.sh'                      % qualified_og
        #genome_tree_file_subset_png         = '%s_genome_tree.treefile.png'         % qualified_og
        #gene_tree_treefile_subset_png       = '%s_subset.treefile.png'              % qualified_og
        #combined_image                      = '%s_combined_trees.png'               % qualified_og
        #genome_tree_file_subset_pdf         = '%s_genome_tree.treefile.pdf'         % qualified_og
        #gene_tree_treefile_subset_pdf       = '%s_subset.treefile.pdf'              % qualified_og
        #combined_pdf                        = '%s_combined_trees.pdf'               % qualified_og

        pwd_genome_tree_file_subset         = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, genome_tree_file_subset)
        pwd_genome_tree_file_subset_for_ale = '%s/%s'                               % (ale_wd, genome_tree_file_subset_for_ale)
        pwd_gene_tree_ufboot                = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, gene_tree_ufboot)
        pwd_gene_tree_ufboot_for_ale        = '%s/%s'                               % (ale_wd, gene_tree_ufboot_for_ale)
        pwd_gene_tree_treefile              = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, gene_tree_treefile)
        pwd_gene_tree_treefile_subset       = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, gene_tree_treefile_subset)
        #pwd_genome_tree_file_subset_png     = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, genome_tree_file_subset_png)
        #pwd_gene_tree_treefile_subset_png   = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, gene_tree_treefile_subset_png)
        #pwd_combined_image                  = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, combined_image)
        #pwd_genome_tree_file_subset_pdf     = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, genome_tree_file_subset_pdf)
        #pwd_gene_tree_treefile_subset_pdf   = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, gene_tree_treefile_subset_pdf)
        #pwd_combined_pdf                    = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, combined_pdf)
        pwd_js_ale                          = '%s/%s'                               % (ale_wd, js_ale)

        if os.path.isfile(pwd_gene_tree_ufboot) is False:
            print('%s not found, please build gene tree first!' % gene_tree_ufboot)
        else:
            print('%s (%s/%s): preparing files and job script for running ALE' % (qualified_og, n, len(og_to_process)))
            n += 1

            # get genomes on gene tree
            gene_gnm_set = set()
            gnm_to_gene_dict = dict()
            for each_gene in Tree(pwd_gene_tree_treefile).get_leaf_names():
                gene_gnm = each_gene.split('.gtdb')[0]
                gene_gnm_set.add(gene_gnm)
                if gene_gnm not in gnm_to_gene_dict:
                    gnm_to_gene_dict[gene_gnm] = {each_gene}
                else:
                    gnm_to_gene_dict[gene_gnm].add(each_gene)

            # subset genome tree
            genome_tree_leaf_set = Tree(genome_tree_file_rooted).get_leaf_names()
            gnms_in_both_trees = set(genome_tree_leaf_set).intersection(gene_gnm_set)
            gnm_tree_subset_str = subset_tree(genome_tree_file_rooted, gnms_in_both_trees, None)
            gnm_tree_subset_str_for_ale = gnm_tree_subset_str
            gnm_tree_subset_str_for_ale = gnm_tree_subset_str_for_ale.replace('GCA_', 'GCA').replace('GCF_', 'GCF')

            # write out genome tree subset
            with open(pwd_genome_tree_file_subset, 'w') as pwd_genome_tree_file_subset_handle:
                pwd_genome_tree_file_subset_handle.write(gnm_tree_subset_str)

            # write out genome tree subset for running ALE
            with open(pwd_genome_tree_file_subset_for_ale, 'w') as pwd_genome_tree_file_subset_for_ale_handle:
                pwd_genome_tree_file_subset_for_ale_handle.write(gnm_tree_subset_str_for_ale)

            # get genes to keep in gene tree
            gene_set_to_keep = set()
            for each_gnm in gnms_in_both_trees:
                gene_set_to_keep.update(gnm_to_gene_dict.get(each_gnm, set()))

            # subset gene_tree.treefile
            subset_tree(pwd_gene_tree_treefile, gene_set_to_keep, pwd_gene_tree_treefile_subset)

            # subset gene_tree.ufboot and rename leaves for running ALE
            pwd_gene_tree_ufboot_for_ale_handle = open(pwd_gene_tree_ufboot_for_ale, 'w')
            for each_gene_tree in open(pwd_gene_tree_ufboot):
                gene_tree_str = each_gene_tree.strip()
                gene_tree_str_subset_for_ale = subset_tree(gene_tree_str, gene_set_to_keep, None)
                gene_tree_str_subset_for_ale = gene_tree_str_subset_for_ale.replace('GCA_', 'GCA').replace('GCF_', 'GCF').replace('.gtdb', '')
                pwd_gene_tree_ufboot_for_ale_handle.write(gene_tree_str_subset_for_ale + '\n')
            pwd_gene_tree_ufboot_for_ale_handle.close()

            # get gene tree leaf name dict (for plot)
            leaf_name_dict = dict()
            for each_gene in Tree(pwd_gene_tree_treefile_subset).get_leaf_names():
                gene_id = each_gene
                gene_id_new = gene_id
                if '.gtdb' in gene_id_new:
                    gene_id_new = gene_id_new.replace('.gtdb', '')
                gene_genome = '_'.join(gene_id_new.split('_')[:-1])
                genome_pco = gnm_pco_dict[gene_genome]
                gene_id_with_taxon = '%s_%s' % (genome_pco, gene_id.split('_')[-1])
                leaf_name_dict[gene_id] = gene_id_with_taxon

            # plot and combine genome and gene trees
            #plot_tree(pwd_genome_tree_file_subset, 'Species tree (rooted)', gnm_pco_dict, dict(), align_leaf_name, show_scale, pwd_genome_tree_file_subset_png)
            #plot_tree(pwd_gene_tree_treefile_subset, 'Gene tree (unrooted)', leaf_name_dict, dict(), align_leaf_name, show_scale, pwd_gene_tree_treefile_subset_png)
            #plot_tree(pwd_genome_tree_file_subset, 'Species tree (rooted)', gnm_pco_dict, dict(), align_leaf_name, show_scale, pwd_genome_tree_file_subset_pdf)
            #plot_tree(pwd_gene_tree_treefile_subset, 'Gene tree (unrooted)', leaf_name_dict, dict(), align_leaf_name, show_scale, pwd_gene_tree_treefile_subset_pdf)

            #merge_image([pwd_genome_tree_file_subset_png, pwd_gene_tree_treefile_subset_png], pwd_combined_image)
            #merge_pdf(pwd_genome_tree_file_subset_pdf, pwd_gene_tree_treefile_subset_pdf, pwd_combined_pdf)
            #os.system('rm %s' % pwd_genome_tree_file_subset_png)
            #os.system('rm %s' % pwd_gene_tree_treefile_subset_png)

            # prepare job script for running ALE
            with open(pwd_js_ale, 'w') as pwd_js_ale_handle:
                obtain_ale_file_cmd = 'ALEobserve %s'                      % gene_tree_ufboot_for_ale
                reconciliation_cmd  = 'ALEml_undated %s %s.ale'            % (genome_tree_file_subset_for_ale, gene_tree_ufboot_for_ale)
                ale_splitter_cmd    = 'python3 %s -i %s.ale.uml_rec -sftr' % (ale_splitter_py, gene_tree_ufboot_for_ale)
                pwd_js_ale_handle.write('#!/bin/bash\n')
                pwd_js_ale_handle.write(obtain_ale_file_cmd + '\n')
                pwd_js_ale_handle.write(reconciliation_cmd + '\n')
                pwd_js_ale_handle.write(ale_splitter_cmd + '\n')

########################################################################################################################

if os.path.isdir(ale_hgt_plot_dir) is True:
    os.system('rm -r %s' % ale_hgt_plot_dir)
os.system('mkdir %s' % ale_hgt_plot_dir)
os.system('mkdir %s/without_HGT' % ale_hgt_plot_dir)

# parse ALE output
n = 1
for qualified_og in og_to_process:

    print('%s (%s/%s): Parsing ALE outputs' % (qualified_og, n, len(og_to_process)))
    n += 1
    gene_tree_treefile                                  = '%s.treefile'                                         % qualified_og
    genome_tree_file_subset_for_ale                     = '%s_genome_tree_for_ALE.treefile'                     % qualified_og
    genome_tree_file_subset_for_ale_png                 = '%s_genome_tree_for_ALE.treefile.png'                 % qualified_og
    gene_tree_ufboot_for_ale                            = '%s_for_ALE.ufboot'                                   % qualified_og
    uts_file                                            = '%s.ale.uTs'                                          % gene_tree_ufboot_for_ale
    uml_rec_file                                        = '%s.ale.uml_rec'                                      % gene_tree_ufboot_for_ale
    ale_formatted_gnm_tree                              = '%s_ALE_formatted_genome_tree.tree'                   % gene_tree_ufboot_for_ale
    ale_formatted_gnm_tree_with_len                     = '%s_ALE_formatted_genome_tree_with_len.tree'          % gene_tree_ufboot_for_ale
    ale_formatted_gnm_tree_with_len_prefixed            = '%s_ALE_formatted_genome_tree_with_len_prefixed.tree' % gene_tree_ufboot_for_ale
    ale_formatted_gnm_tree_with_len_prefixed_png        = '%s_genome_tree_with_HGT.pdf'                         % qualified_og
    itol_connection_txt                                 = '%s_iTOL_connection.txt'                              % qualified_og
    itol_label_txt                                      = '%s_iTOL_genome_pco.txt'                              % qualified_og
    gene_tree_itol_label_txt                            = '%s_iTOL_gene_pco.txt'                                % qualified_og
    gene_tree_treefile_subset_png                       = '%s_subset.treefile.pdf'                              % qualified_og
    combined_image_with_ale_hgts                        = '%s_combined_tree_with_HGTs.pdf'                      % qualified_og
    gene_tree_treefile_subset                           = '%s_subset.treefile'                                  % qualified_og
    gnm_tree_label_color_txt                            = '%s_genome_tree_label_color.txt'                      % qualified_og
    gene_tree_label_color_txt                           = '%s_gene_tree_label_color.txt'                        % qualified_og
    pwd_gene_tree_treefile_subset                       = '%s/%s/%s'                                            % (gene_tree_dir, qualified_og, gene_tree_treefile_subset)
    pwd_gene_tree_treefile_subset_png                   = '%s/%s/%s'                                            % (gene_tree_dir, qualified_og, gene_tree_treefile_subset_png)
    pwd_gene_tree_treefile                              = '%s/%s/%s'                                            % (gene_tree_dir, qualified_og, gene_tree_treefile)
    pwd_genome_tree_file_subset_for_ale                 = '%s/%s'                                               % (ale_op_dir, genome_tree_file_subset_for_ale)
    pwd_gene_tree_ufboot_for_ale                        = '%s/%s'                                               % (ale_wd, gene_tree_ufboot_for_ale)
    pwd_genome_tree_file_subset_for_ale_png             = '%s/%s'                                               % (ale_op_dir, genome_tree_file_subset_for_ale_png)
    pwd_itol_connection_txt                             = '%s/%s'                                               % (ale_hgt_plot_dir, itol_connection_txt)
    pwd_itol_label_txt                                  = '%s/%s'                                               % (ale_op_dir, itol_label_txt)
    pwd_gene_tree_itol_label_txt                        = '%s/%s'                                               % (ale_hgt_plot_dir, gene_tree_itol_label_txt)
    pwd_uts_file                                        = '%s/%s'                                               % (ale_op_dir, uts_file)
    pwd_uml_rec_file                                    = '%s/%s'                                               % (ale_op_dir, uml_rec_file)
    pwd_ale_formatted_gnm_tree                          = '%s/%s'                                               % (ale_op_dir, ale_formatted_gnm_tree)
    pwd_ale_formatted_gnm_tree_with_len                 = '%s/%s'                                               % (ale_op_dir, ale_formatted_gnm_tree_with_len)
    pwd_ale_formatted_gnm_tree_with_len_prefixed        = '%s/%s'                                               % (ale_op_dir, ale_formatted_gnm_tree_with_len_prefixed)
    pwd_ale_formatted_gnm_tree_with_len_prefixed_png    = '%s/%s'                                               % (ale_wd, ale_formatted_gnm_tree_with_len_prefixed_png)
    pwd_combined_image_with_ale_hgts                    = '%s/%s'                                               % (ale_hgt_plot_dir, combined_image_with_ale_hgts)
    pwd_gnm_tree_label_color_txt                        = '%s/%s'                                               % (ale_hgt_plot_dir, gnm_tree_label_color_txt)
    pwd_gene_tree_label_color_txt                       = '%s/%s'                                               % (ale_hgt_plot_dir, gene_tree_label_color_txt)

    # uts_to_itol_connections
    qualified_hgt_num = 0
    internal_node_to_leaf_dict = dict()
    paired_donor_to_recipient_leaf_dict = dict()
    if os.path.isfile(pwd_uts_file) is True:

        # write out ALE formatted genome tree
        renamed_genome_tree_str = open(pwd_uml_rec_file).readlines()[2].strip().split('\t')[1]
        with open(pwd_ale_formatted_gnm_tree, 'w') as ale_renamed_species_tree_handle:
            ale_renamed_species_tree_handle.write(renamed_genome_tree_str + '\n')

        qualified_hgt_num, internal_node_to_leaf_dict, paired_donor_to_recipient_leaf_dict = uts_to_itol_connections(pwd_genome_tree_file_subset_for_ale, pwd_ale_formatted_gnm_tree, interal_node_prefix, pwd_uts_file, hgt_freq_cutoff, ignore_leaf_hgt, ignore_vertical_hgt, donor_node_min_leaf_num, recipient_node_min_leaf_num, pwd_itol_connection_txt)

    else:
        print('%s: uTs file not found, you need to run ALE first!' % qualified_og)

    # combine_trees
    combine_trees(pwd_genome_tree_file_subset_for_ale, pwd_ale_formatted_gnm_tree, pwd_ale_formatted_gnm_tree_with_len)

    # prefix_internal_nodes of combined tree
    prefix_internal_nodes(pwd_ale_formatted_gnm_tree_with_len, interal_node_prefix, pwd_ale_formatted_gnm_tree_with_len_prefixed)

    g_to_dr_dict = dict()
    for each_d2r in paired_donor_to_recipient_leaf_dict:
        d_gene_list = paired_donor_to_recipient_leaf_dict[each_d2r][0]
        r_gene_list = paired_donor_to_recipient_leaf_dict[each_d2r][1]
        for each_d in d_gene_list:
            if each_d not in g_to_dr_dict:
                g_to_dr_dict[each_d] = {'d'}
            else:
                g_to_dr_dict[each_d].add('d')
        for each_r in r_gene_list:
            if each_r not in g_to_dr_dict:
                g_to_dr_dict[each_r] = {'r'}
            else:
                g_to_dr_dict[each_r].add('r')

    # write out iTOL label file for gene and genome tree
    pwd_itol_label_txt_handle  = open(pwd_itol_label_txt, 'w')
    pwd_itol_label_txt_handle.write('LABELS\nSEPARATOR TAB\n\nDATA\n')
    pwd_gene_tree_itol_label_txt_handle = open(pwd_gene_tree_itol_label_txt, 'w')
    pwd_gene_tree_itol_label_txt_handle.write('LABELS\nSEPARATOR TAB\n\nDATA\n')
    pwd_gene_tree_label_color_txt_handle = open(pwd_gene_tree_label_color_txt, 'w')
    pwd_gene_tree_label_color_txt_handle.write('DATASET_STYLE\nSEPARATOR TAB\nDATASET_LABEL\texample_style\nCOLOR\t#ffff00\n\nDATA\n')
    wrote_gnm_set = set()
    for each_gene in Tree(pwd_gene_tree_treefile).get_leaf_names():
        gene_gnm = each_gene.split('.gtdb')[0]
        gene_name_for_ale = gene_gnm
        gene_name_for_ale = gene_name_for_ale.replace('GCA_', 'GCA').replace('GCF_', 'GCF')
        gene_with_taxon = gnm_pco_dict[gene_gnm]

        gnm_dr = g_to_dr_dict.get(gene_name_for_ale, set())
        if gnm_dr == {'d'}:
            pwd_gene_tree_label_color_txt_handle.write('%s\tlabel\tnode\t%s\t1\tnormal\n' % (each_gene, d_color))
        elif gnm_dr == {'r'}:
            pwd_gene_tree_label_color_txt_handle.write('%s\tlabel\tnode\t%s\t1\tnormal\n' % (each_gene, r_color))
        elif len(gnm_dr) == 2:
            pwd_gene_tree_label_color_txt_handle.write('%s\tlabel\tnode\t%s\t1\tnormal\n' % (each_gene, '#FF7F50'))

        if gene_gnm not in wrote_gnm_set:
            pwd_itol_label_txt_handle.write('%s\t%s\n' % (gene_name_for_ale, gene_with_taxon))
            wrote_gnm_set.add(gene_gnm)
        pwd_gene_tree_itol_label_txt_handle.write('%s\t%s_%s\n' % (each_gene, gene_with_taxon, each_gene.split('_')[-1]))
    pwd_itol_label_txt_handle.close()
    pwd_gene_tree_itol_label_txt_handle.close()
    pwd_gene_tree_label_color_txt_handle.close()

    # write out gnm_tree_label_color_txt
    pwd_gnm_tree_label_color_txt_handle = open(pwd_gnm_tree_label_color_txt, 'w')
    pwd_gnm_tree_label_color_txt_handle.write('DATASET_STYLE\nSEPARATOR TAB\nDATASET_LABEL\texample_style\nCOLOR\t#ffff00\n\nDATA\n')
    for d2r in paired_donor_to_recipient_leaf_dict:
        pwd_gnm_tree_label_color_txt_handle.write('%s\tlabel\tclade\t%s\t1\tnormal\n' % (d2r.split('___')[0], d_color))
        pwd_gnm_tree_label_color_txt_handle.write('%s\tlabel\tclade\t%s\t1\tnormal\n' % (d2r.split('___')[1], r_color))
    pwd_gnm_tree_label_color_txt_handle.close()

    if qualified_hgt_num == 0:
        os.system('mv %s/%s* %s/without_HGT/' % (ale_hgt_plot_dir, qualified_og, ale_hgt_plot_dir))
    else:
        itol_tree(pwd_ale_formatted_gnm_tree_with_len_prefixed, [pwd_gnm_tree_label_color_txt, pwd_itol_label_txt, pwd_itol_connection_txt], project_name, API_key, display_mode, pwd_ale_formatted_gnm_tree_with_len_prefixed_png)
        itol_tree(pwd_gene_tree_treefile_subset, [pwd_gene_tree_itol_label_txt, pwd_gene_tree_label_color_txt], project_name, API_key, display_mode, pwd_gene_tree_treefile_subset_png)
        merge_pdf(pwd_ale_formatted_gnm_tree_with_len_prefixed_png, pwd_gene_tree_treefile_subset_png, pwd_combined_image_with_ale_hgts)

########################################################################################################################
