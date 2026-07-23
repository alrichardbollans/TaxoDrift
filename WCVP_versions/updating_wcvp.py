import os

import pandas as pd
from wcvpy.wcvp_download import get_all_taxa, add_authors_to_col

from chaining_methods import compare_two_versions, chain_two_databases, get_direct_name_updates, compare_and_output_chained_and_direct_updates, \
    summarise_results, get_overrepresented_genera

repo_path = os.environ.get('KEWSCRATCHPATH')
this_repo_path = os.path.join(repo_path, 'TaxoDrift')
_output_path = os.path.join(this_repo_path, 'WCVP_versions', 'outputs')
_input_path = os.path.join(this_repo_path, 'WCVP_versions', 'inputs')

wcvp_version_order = ['v10', 'v11', 'v12', 'v13', 'v14', 'v15', 'v16']

if not os.path.isdir(_output_path):
    os.mkdir(_output_path)


def compare_all_pairs():
    v10_taxa, v11_taxa, v12_taxa, v13_taxa, v14_taxa, v15_taxa, v16_taxa = get_all_databases()

    compare_two_versions(v10_taxa, v11_taxa,
                         'v10', 'v11', _output_path)
    compare_two_versions(v10_taxa, v12_taxa,
                         'v10', 'v12', _output_path)
    compare_two_versions(v10_taxa, v13_taxa,
                         'v10', 'v13', _output_path)

    compare_two_versions(v11_taxa, v12_taxa,
                         'v11', 'v12', _output_path)
    compare_two_versions(v11_taxa, v13_taxa,
                         'v11', 'v13', _output_path)

    compare_two_versions(v12_taxa, v13_taxa, 'v12', 'v13', _output_path)

    ## 14
    compare_two_versions(v10_taxa, v14_taxa, 'v10', 'v14', _output_path)
    compare_two_versions(v11_taxa, v14_taxa, 'v11', 'v14', _output_path)
    compare_two_versions(v12_taxa, v14_taxa, 'v12', 'v14', _output_path)
    compare_two_versions(v13_taxa, v14_taxa, 'v13', 'v14', _output_path)

    ## 15
    compare_two_versions(v10_taxa, v15_taxa, 'v10', 'v15', _output_path)
    compare_two_versions(v11_taxa, v15_taxa, 'v11', 'v15', _output_path)
    compare_two_versions(v12_taxa, v15_taxa, 'v12', 'v15', _output_path)
    compare_two_versions(v13_taxa, v15_taxa, 'v13', 'v15', _output_path)
    compare_two_versions(v14_taxa, v15_taxa, 'v14', 'v15', _output_path)

    ## 16
    compare_two_versions(v10_taxa, v16_taxa, 'v10', 'v16', _output_path)
    compare_two_versions(v11_taxa, v16_taxa, 'v11', 'v16', _output_path)
    compare_two_versions(v12_taxa, v16_taxa, 'v12', 'v16', _output_path)
    compare_two_versions(v13_taxa, v16_taxa, 'v13', 'v16', _output_path)
    compare_two_versions(v14_taxa, v16_taxa, 'v14', 'v16', _output_path)
    compare_two_versions(v15_taxa, v16_taxa, 'v15', 'v16', _output_path)


def full_chain_results():
    # Note when chaining like this, in intermediary steps ambiguous/non resolving names may be dropped.
    # This may somewhat reflect real world situations but is optimistic about the chaining process
    out_dir = os.path.join('outputs', 'full_chain')
    v10_taxa, v11_taxa, v12_taxa, v13_taxa, v14_taxa, v15_taxa, v16_taxa = get_all_databases()

    def rename_columns_after_chaining(df, new_tag):
        df = df.rename(columns={f'{new_tag}_chained_accepted_name_w_author': 'accepted_name_w_author'})
        df = df[['taxon_name_w_authors', 'accepted_name_w_author']]
        return df

    # Start with 10 -> 11
    v10_11_chained = chain_two_databases(v10_taxa, v11_taxa, 'v10', 'v11', out_dir)
    v10_11_chained = rename_columns_after_chaining(v10_11_chained, 'v11')

    # Then chain -> 12
    v10_11_12_chained = chain_two_databases(v10_11_chained, v12_taxa, 'v10_11', 'v12', out_dir)
    v10_11_12_chained = rename_columns_after_chaining(v10_11_12_chained, 'v12')

    # Then -> 13
    v10_11_12_13_chained = chain_two_databases(v10_11_12_chained, v13_taxa, 'v10_11_12', 'v13', out_dir)
    v10_11_12_13_chained = rename_columns_after_chaining(v10_11_12_13_chained, 'v13')

    # Then -> 14
    v10_11_12_13_14_chained = chain_two_databases(v10_11_12_13_chained, v14_taxa, 'v10_11_12_13', 'v14', out_dir)
    v10_11_12_13_14_chained = rename_columns_after_chaining(v10_11_12_13_14_chained, 'v14')

    # Then -> 15
    v10_11_12_13_14_15_chained = chain_two_databases(v10_11_12_13_14_chained, v15_taxa, 'v10_11_12_13_14', 'v15', out_dir)
    v10_11_12_13_14_15_chained = rename_columns_after_chaining(v10_11_12_13_14_15_chained, 'v15')

    # Then -> 16
    v10_11_12_13_14_15_16_chained = chain_two_databases(v10_11_12_13_14_15_chained, v16_taxa, 'v10_11_12_13_14_15', 'v16', out_dir)


    direct_updated_records = get_direct_name_updates(v10_taxa, v16_taxa, 'v16', out_dir)
    results_df = compare_and_output_chained_and_direct_updates(v10_11_12_13_14_15_16_chained, direct_updated_records,
                                                               'v10_11_12_13_14_15', 'v16', out_dir)
    pass


def get_all_databases(do_summaries=False):
    # v10_taxa = get_all_taxa(version='10', output_csv=os.path.join(_input_path, 'v10_taxa.csv'))
    # v11_taxa = get_all_taxa(version='11', output_csv=os.path.join(_input_path, 'v11_taxa.csv'))
    # v12_taxa = get_all_taxa(version='12', output_csv=os.path.join(_input_path, 'v12_taxa.csv'))
    # v13_taxa = get_all_taxa(version='13', output_csv=os.path.join(_input_path, 'v13_taxa.csv'))
    # v14_taxa = get_all_taxa(version='14', output_csv=os.path.join(_input_path, 'v14_taxa.csv'))
    # v15_taxa = get_all_taxa(version='15', output_csv=os.path.join(_input_path, 'v15_taxa.csv'))
    # v16_taxa = get_all_taxa(get_new_version=True, output_csv=os.path.join(_input_path, 'v16_taxa.csv'))

    v10_taxa = pd.read_csv(os.path.join(_input_path, 'v10_taxa.csv'), index_col=0)
    v11_taxa = pd.read_csv(os.path.join(_input_path, 'v11_taxa.csv'), index_col=0)
    v12_taxa = pd.read_csv(os.path.join(_input_path, 'v12_taxa.csv'), index_col=0)
    v13_taxa = pd.read_csv(os.path.join(_input_path, 'v13_taxa.csv'), index_col=0)
    v14_taxa = pd.read_csv(os.path.join(_input_path, 'v14_taxa.csv'), index_col=0)
    v15_taxa = pd.read_csv(os.path.join(_input_path, 'v15_taxa.csv'), index_col=0)
    v16_taxa = pd.read_csv(os.path.join(_input_path, 'v16_taxa.csv'), index_col=0)

    v10_taxa['taxon_name_w_authors'] = add_authors_to_col(v10_taxa, 'taxon_name')
    v11_taxa['taxon_name_w_authors'] = add_authors_to_col(v11_taxa, 'taxon_name')
    v12_taxa['taxon_name_w_authors'] = add_authors_to_col(v12_taxa, 'taxon_name')
    v13_taxa['taxon_name_w_authors'] = add_authors_to_col(v13_taxa, 'taxon_name')
    v14_taxa['taxon_name_w_authors'] = add_authors_to_col(v14_taxa, 'taxon_name')
    v15_taxa['taxon_name_w_authors'] = add_authors_to_col(v15_taxa, 'taxon_name')
    v16_taxa['taxon_name_w_authors'] = add_authors_to_col(v16_taxa, 'taxon_name')

    if do_summaries:
        v10_taxa.describe(include='all').to_csv(os.path.join(_input_path, 'v10_taxa_summary.csv'))
        v11_taxa.describe(include='all').to_csv(os.path.join(_input_path, 'v11_taxa_summary.csv'))
        v12_taxa.describe(include='all').to_csv(os.path.join(_input_path, 'v12_taxa_summary.csv'))
        v13_taxa.describe(include='all').to_csv(os.path.join(_input_path, 'v13_taxa_summary.csv'))
        v14_taxa.describe(include='all').to_csv(os.path.join(_input_path, 'v14_taxa_summary.csv'))
        v15_taxa.describe(include='all').to_csv(os.path.join(_input_path, 'v15_taxa_summary.csv'))
        v16_taxa.describe(include='all').to_csv(os.path.join(_input_path, 'v16_taxa_summary.csv'))

    return v10_taxa, v11_taxa, v12_taxa, v13_taxa, v14_taxa, v15_taxa, v16_taxa


def main():
    get_all_databases(do_summaries=True)

    compare_all_pairs()
    full_chain_results()
    summarise_results(os.path.join(_output_path, f'full_chain'), f'v10_11_12_13_14_15_v16', old_tag='v10')
    for w in wcvp_version_order:
        for w2 in wcvp_version_order:
            try:
                summarise_results(os.path.join(_output_path, f'{w}_{w2}'), f'{w}_{w2}', old_tag=w)
            except:
                print(f'Could not summarise {w}, {w2}')

    # Genus results
    v10_taxa, v11_taxa, v12_taxa, v13_taxa, v14_taxa, v15_taxa, v16_taxa = get_all_databases()
    genus_counts = get_overrepresented_genera(_output_path, 'v10', 'v16', v10_taxa)
    print(genus_counts)
    genus_counts.to_csv(os.path.join(_output_path, f'v10_v16', 'genus_counts.csv'))


if __name__ == '__main__':
    main()
