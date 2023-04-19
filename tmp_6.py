import os
og_to_process = []


qualified_og = ''
gene_tree_dir = ''
ale_wd = ''
ale_op_dir = ''
ale_hgt_plot_dir = ''
interal_node_prefix = ''
gnm_pco_dict = dict()
d_color = ''
r_color = ''
project_name = ''
API_key = ''
display_mode = ''
hgt_freq_cutoff = ''
ignore_leaf_hgt= ''
ignore_vertical_hgt= ''
donor_node_min_leaf_num= ''
recipient_node_min_leaf_num= ''


# parse ALE output

def parse_ale_op_worker(arg_list):

    qualified_og                = arg_list[0]
    gene_tree_dir               = arg_list[1]
    ale_wd                      = arg_list[2]
    ale_op_dir                  = arg_list[3]
    ale_hgt_plot_dir            = arg_list[4]
    interal_node_prefix         = arg_list[5]
    gnm_pco_dict                = arg_list[6]
    d_color                     = arg_list[7]
    r_color                     = arg_list[8]
    project_name                = arg_list[9]
    API_key                     = arg_list[10]
    display_mode                = arg_list[11]
    hgt_freq_cutoff             = arg_list[12]
    ignore_leaf_hgt             = arg_list[13]
    ignore_vertical_hgt         = arg_list[14]
    donor_node_min_leaf_num     = arg_list[15]
    recipient_node_min_leaf_num = arg_list[16]

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

