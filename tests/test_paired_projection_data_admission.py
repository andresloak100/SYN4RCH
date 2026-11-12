import json

from analysis import paired_projection_data_admission as A

SOFT = """^SERIES = GSE1
!Series_supplementary_file = ftp://x/suppl/GSE1_RAW.tar
^SAMPLE = GSM1
!Sample_title = LHA sample pool 1
!Sample_description = pooled live single-cell suspension generated from LHA tissue of two animals (male and female)
!Sample_characteristics_ch1 = strain: C57BL6/J
!Sample_supplementary_file_1 = ftp://x/GSM1_LHA_A.h5
"""

FILELIST = "#Archive/File\tName\tTime\tSize\tType\nFile\tGSM1_LHA_A.h5\tt\t56741529\tH5\n"


def nb(files, printed):
    return json.dumps({'cells': [
        {'cell_type': 'code', 'source': ["filenames=[" + ",".join(f"'{f}'" for f in files) + "]\n",
                                         "df['MouseSlideSlice']=1\n", "df.loc[0,'AP coord']=-1\n",
                                         "x(test_data=.5)\n"],
         'outputs': [{'text': [printed]}]}]})


def raw(repo_paths, nb_files=('AG09_C1_quant.xlsx', 'AG10_A2_quant.xlsx')):
    return {
        'rossi_geo_soft': SOFT, 'rossi_geo_filelist': FILELIST,
        'rossi_github_tree': json.dumps({'tree': [{'type': 'blob', 'path': p} for p in repo_paths]}),
        'rossi_github_repo': json.dumps({'license': None}),
        'rossi_hcr_notebook': nb(nb_files, 'Total Cells:  92866\nNumber of YFP cells:  816\n'
                                           'Number of Tomato cells:  857\n'),
        'rossi_pmc_xml': '<p>HCR C57BL/6J mice (n = 3)</p>', 'rossi_epmc_supp': 'not open access',
        'luo_pmc_xml': ('<p>Two cells were collected in one tube as a sample.</p>'
                        '<p>Among the four eGFP-labeled cells recovered from SuM-vDG mice</p>'
                        '<p>Cells with eGFP count values greater than 0 were considered eGFP-positive</p>'),
    }


def test_geo_and_filelist_parse():
    g = A.parse_geo_soft(SOFT)
    assert g['series_files'] == ['GSE1_RAW.tar']
    assert g['samples'][0]['gsm'] == 'GSM1' and 'two animals' in g['samples'][0]['description']
    assert A.parse_filelist(FILELIST)[0]['bytes'] == 56741529


def test_notebook_extracts_inputs_and_counts_only():
    n = A.parse_notebook(nb(['AG09_C1_quant.xlsx', 'AG11_C2_quant.xlsx'], 'Total Cells:  5\n'))
    assert n['input_files'] == ['AG09_C1_quant.xlsx', 'AG11_C2_quant.xlsx']
    assert n['hcr_mouse_ids_from_filenames'] == ['AG09', 'AG11']
    assert n['printed_total_cells'] == [5]
    assert n['svm_split_random_50pct_not_by_animal']


def test_rossi_rejected_when_hcr_tables_absent():
    r = A.audit(raw(['README.md', 'Fig3_S2_HCR analysis_Final.ipynb']))
    assert r['rossi2021_lha']['verdict'] == 'REJECT'
    ro = r['rossi2021_lha']
    assert ro['facts']['hcr_per_cell_counts_found_in_checked_sources'] is False
    assert 'NOT FOUND in the checked public sources' in ro['hard_blockers'][0]
    assert 'not a claim that the data are universally non-public' in ro['hard_blockers'][0]
    assert ro['unretrieved_sources']
    lu = r['luo2025_sum']
    assert lu['verdict'] == 'REJECT'
    assert lu['quotes']['patchseq_two_cells_per_tube']
    assert lu['facts']['snrna_label_detected_per_nucleus'] is True
    assert lu['quotes']['snrna_egfp_positives'] and lu['quotes']['snrna_egfp_threshold']
    assert not any('Epi-retro-seq' in b for b in lu['hard_blockers'])


def test_rna_pools_counted_from_geo_never_as_recovered_animals():
    f = A.audit(raw(['README.md']))['rossi2021_lha']['facts']
    assert f['rna_geo_pools_of_two_animals'] == 1          # one sample in the fixture
    assert 'GEO' in f['rna_pooling_sources'] and 'paper' in f['rna_pooling_sources']
    assert 'never six recovered animal IDs' in f['rna_pools_independence']
    assert 'rna_mice' not in f


def test_cached_mode_refuses_changed_raw(tmp_path, monkeypatch):
    monkeypatch.setattr(A, 'RAW', tmp_path)
    (tmp_path / 'x').write_bytes(b'abc')
    good = {'name': 'x', 'sha256': A.sha256_bytes(b'abc')}
    assert A.load_cached([good]) == {'x': 'abc'}
    import pytest
    with pytest.raises(ValueError):
        A.load_cached([{'name': 'x', 'sha256': '0' * 64}])


def test_public_hcr_tables_only_move_to_pending_never_admit():
    r = A.audit(raw(['data/AG09_C1_quant.xlsx']))
    assert r['rossi2021_lha']['verdict'] == 'PENDING'


def test_rna_deg_csv_in_repo_does_not_count_as_hcr():
    nbf = ('AG09_C1_quant.xlsx', 'YFPvTOM_DEGs_ALL_UPDATED210130.csv', 'HCR_data.csv')
    r = A.audit(raw(['out/YFPvTOM_DEGs_ALL_UPDATED210130.csv',
                     'Seurat_SCT_YFP_Cells_subset_genes_UPDATED210130.csv'], nb_files=nbf))
    assert r['rossi2021_lha']['verdict'] == 'REJECT'
    assert r['rossi2021_lha']['hcr_candidate_files_found'] == []


def test_hcr_candidate_names():
    assert A.is_hcr_candidate('AG10_A2_quant.xlsx') and A.is_hcr_candidate('HCR_data.csv')
    assert not A.is_hcr_candidate('HRC_seurat_deg_output.csv')
    assert not A.is_hcr_candidate('Correlate_logFC_SEQvHCR.csv')


def test_fetch_cap_is_50mb():
    assert A.MAX_BYTES == 50 * 1024 * 1024
