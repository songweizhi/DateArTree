
def prepare_ale_ip_worker(arg_list):

    qualified_og            = arg_list[0]
    gene_tree_dir           = arg_list[1]
    ale_wd                  = arg_list[2]
    genome_tree_file_rooted = arg_list[3]
    gnm_pco_dict            = arg_list[4]
    ale_splitter_py         = arg_list[5]

    genome_tree_file_subset             = '%s_genome_tree.treefile'             % qualified_og
    genome_tree_file_subset_for_ale     = '%s_genome_tree_for_ALE.treefile'     % qualified_og
    gene_tree_ufboot                    = '%s.ufboot'                           % qualified_og
    gene_tree_ufboot_for_ale            = '%s_for_ALE.ufboot'                   % qualified_og
    gene_tree_treefile                  = '%s.treefile'                         % qualified_og
    gene_tree_treefile_subset           = '%s_subset.treefile'                  % qualified_og
    js_ale                              = '/zjs_%s_ALE.sh'                      % qualified_og
    pwd_genome_tree_file_subset         = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, genome_tree_file_subset)
    pwd_genome_tree_file_subset_for_ale = '%s/%s'                               % (ale_wd, genome_tree_file_subset_for_ale)
    pwd_gene_tree_ufboot                = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, gene_tree_ufboot)
    pwd_gene_tree_ufboot_for_ale        = '%s/%s'                               % (ale_wd, gene_tree_ufboot_for_ale)
    pwd_gene_tree_treefile              = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, gene_tree_treefile)
    pwd_gene_tree_treefile_subset       = '%s/%s/%s'                            % (gene_tree_dir, qualified_og, gene_tree_treefile_subset)
    pwd_js_ale                          = '%s/%s'                               % (ale_wd, js_ale)


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

    # prepare job script for running ALE
    with open(pwd_js_ale, 'w') as pwd_js_ale_handle:
        obtain_ale_file_cmd = 'ALEobserve %s'                      % gene_tree_ufboot_for_ale
        reconciliation_cmd  = 'ALEml_undated %s %s.ale'            % (genome_tree_file_subset_for_ale, gene_tree_ufboot_for_ale)
        ale_splitter_cmd    = 'python3 %s -i %s.ale.uml_rec -sftr' % (ale_splitter_py, gene_tree_ufboot_for_ale)
        pwd_js_ale_handle.write('#!/bin/bash\n')
        pwd_js_ale_handle.write(obtain_ale_file_cmd + '\n')
        pwd_js_ale_handle.write(reconciliation_cmd + '\n')
        pwd_js_ale_handle.write(ale_splitter_cmd + '\n')
