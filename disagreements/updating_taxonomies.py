import os

import pandas as pd
from pkg_resources import resource_filename
from wcvpy.wcvp_download import get_all_taxa, add_authors_to_col
from wcvpy.wcvp_name_matching import get_genus_from_full_name

_output_path = resource_filename(__name__, 'outputs')

if not os.path.isdir(_output_path):
    os.mkdir(_output_path)


def get_accepted_name_from_record(record: pd.DataFrame, reported_name: str):
    if len(record.index) == 0:
        print(f'{reported_name} has no taxon record')
        raise ValueError
    else:
        accepted_names = record['accepted_name_w_author'].dropna().unique().tolist()

        if len(set(accepted_names)) > 1:
            print(f'{reported_name} resolves to more than one taxon record')
            print(record)
            print(accepted_names)
            raise ValueError
        if len(accepted_names) == 1:
            return accepted_names[0]
        else:
            return None


def compare_two_versions(v12_taxa: pd.DataFrame, v13_taxa: pd.DataFrame, old_tag: str, new_tag: str):
    # For all taxa with unique names (inc. author strings) in old taxon database
    # If the name resolves uniquely to a non-nan accepted name in both the old and new database
    # Find the accepted name resolution when the name is resolved first to the old taxonomy then the new taxonomy
    # If no such name exists, add to cases_that_cant_update
    # else, check if this name is the same as the name when directly resolved using the new database. If not, store in out_dict
    out_dir = os.path.join(_output_path,
                           '_'.join([old_tag, new_tag]))
    os.makedirs(out_dir, exist_ok=True)

    v12_taxa['taxon_name_w_authors'] = add_authors_to_col(v12_taxa, 'taxon_name')
    v13_taxa['taxon_name_w_authors'] = add_authors_to_col(v13_taxa, 'taxon_name')

    # Pick non-duplicated names
    unique_names = v12_taxa['taxon_name_w_authors'].dropna().unique().tolist()

    # Collect the records in the old and new database directly
    # relevant records
    v12_records = v12_taxa[v12_taxa['taxon_name_w_authors'].isin(unique_names)][
        ['taxon_name_w_authors', 'accepted_name_w_author']]
    v12_records = v12_records.dropna(subset=['accepted_name_w_author'])
    v12_records = v12_records.drop_duplicates(keep='first')
    v12_records = v12_records.drop_duplicates(subset=['taxon_name_w_authors'],
                                              keep=False)  # Ignore cases with multiple resolutions, as these are ambiguous anyway
    # I think Cubeba Raf. and Lamottea Pomel are examples but they are less interesting for this analysis

    v12_records = v12_records.rename(
        columns={'accepted_name_w_author': old_tag + '_accepted_name_w_author'})
    v12_records.describe(include='all').to_csv(os.path.join(out_dir, 'old_records_summary.csv'))

    # accepted names in old database
    relevant_v13_names_for_chaining = v12_records[old_tag + '_accepted_name_w_author'].dropna().unique().tolist()
    # names in new database where taxon name is accepted name in old database
    v13_records_for_chaining = v13_taxa[v13_taxa['taxon_name_w_authors'].isin(relevant_v13_names_for_chaining)][
        ['taxon_name_w_authors', 'accepted_name_w_author', 'accepted_species_w_author']]
    v13_records_for_chaining = v13_records_for_chaining.drop_duplicates(keep='first')
    v13_records_for_chaining = v13_records_for_chaining.drop_duplicates(subset=['taxon_name_w_authors'],
                                                                        keep=False)  # Ignore cases with multiple resolutions, as these are ambiguous anyway
    v13_records_for_chaining = v13_records_for_chaining.rename(
        columns={'taxon_name_w_authors': new_tag + '_taxon_name_w_authors', 'accepted_name_w_author': new_tag + '_chained_accepted_name_w_author',
                 'accepted_species_w_author': new_tag + '_chained_accepted_species_w_author'})
    # chained result where taxon_name_w_authors are resolved to v13_chained_accepted_name_w_author, mediated by old database
    chained_updated_records = pd.merge(v12_records, v13_records_for_chaining, left_on=old_tag + '_accepted_name_w_author',
                                       right_on=new_tag + '_taxon_name_w_authors')
    chained_updated_records.describe(include='all').to_csv(os.path.join(out_dir, 'chained_updated_records_summary.csv'))

    # relevant names in new database where taxon name is taxon name in old database
    v13_updated_records = v13_taxa[v13_taxa['taxon_name_w_authors'].isin(unique_names)][
        ['taxon_name_w_authors', 'accepted_name_w_author', 'accepted_species_w_author']]
    v13_updated_records = v13_updated_records.dropna(subset=['accepted_name_w_author'])
    v13_updated_records = v13_updated_records.drop_duplicates(keep='first')
    v13_updated_records = v13_updated_records.drop_duplicates(subset=['taxon_name_w_authors'],
                                                              keep=False)  # Ignore cases with multiple resolutions, as these are ambiguous anyway
    v13_updated_records = v13_updated_records.rename(
        columns={'accepted_name_w_author': new_tag + '_direct_accepted_name_w_author',
                 'accepted_species_w_author': new_tag + '_direct_accepted_species_w_author'})
    v13_updated_records.describe(include='all').to_csv(os.path.join(out_dir, 'direct_updated_records_summary.csv'))

    # direct result where taxon_name_w_authors are resolved directly to the new databse
    merged_df = pd.merge(chained_updated_records, v13_updated_records, on='taxon_name_w_authors')

    results_df = merged_df[merged_df[new_tag + '_direct_accepted_name_w_author'] != merged_df[new_tag + '_chained_accepted_name_w_author']]
    results_df = results_df.dropna(subset=[new_tag + '_direct_accepted_name_w_author'])

    # Add longer chains

    results_df.to_csv(os.path.join(out_dir, 'all_results.csv'))
    results_df.describe(include='all').to_csv(os.path.join(out_dir, 'all_results_summary.csv'))

    species_ambiguity_results = results_df.dropna(subset=[new_tag + '_direct_accepted_species_w_author'])
    species_ambiguity_results = species_ambiguity_results[
        species_ambiguity_results[new_tag + '_direct_accepted_species_w_author'] != species_ambiguity_results[
            new_tag + '_chained_accepted_species_w_author']]

    species_ambiguity_results.to_csv(os.path.join(out_dir, 'species_results.csv'))
    species_ambiguity_results.describe(include='all').to_csv(os.path.join(out_dir, 'species_results_summary.csv'))

    species_ambiguity_results[new_tag + '_chained_accepted_genus'] = species_ambiguity_results[new_tag + '_chained_accepted_species_w_author'].apply(
        get_genus_from_full_name)
    species_ambiguity_results[new_tag + '_direct_accepted_genus'] = species_ambiguity_results[new_tag + '_direct_accepted_species_w_author'].apply(
        get_genus_from_full_name)

    genus_ambiguity_results = species_ambiguity_results.dropna(subset=[new_tag + '_direct_accepted_genus'])
    genus_ambiguity_results = genus_ambiguity_results[
        genus_ambiguity_results[new_tag + '_chained_accepted_genus'] != genus_ambiguity_results[new_tag + '_direct_accepted_genus']]
    genus_ambiguity_results.to_csv(os.path.join(out_dir, 'genus_results.csv'))
    genus_ambiguity_results.describe(include='all').to_csv(os.path.join(out_dir, 'genus_results_summary.csv'))
    # cases_that_cant_update_df.to_csv(os.path.join(out_dir, 'v12_v13_cases_cant_update.csv'))
    # names_in_old_with_multiple_resolutions_df.to_csv(os.path.join(out_dir, 'v12_v13_names_in_old_with_multiple_resolutions.csv'))

    return results_df


def main():
    v12_taxa = get_all_taxa(version='12')
    v13_taxa = get_all_taxa(version=None)
    compare_two_versions(v12_taxa, v13_taxa, 'v12', 'v13')
    v10_taxa = get_all_taxa(version='10')
    compare_two_versions(v10_taxa, v13_taxa,
                         'v10', 'v13')
    v11_taxa = get_all_taxa(version='11')
    compare_two_versions(v11_taxa, v13_taxa,
                         'v11', 'v13')

    compare_two_versions(v10_taxa, v11_taxa,
                         'v10', 'v11')
    compare_two_versions(v11_taxa, v12_taxa,
                         'v11', 'v12')


if __name__ == '__main__':
    main()
