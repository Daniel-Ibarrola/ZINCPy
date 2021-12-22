from zincpy.zinc_client import ZincClient
import zincpy._private_tools.exceptions as exc
import pytest

def test_validate_filters():
    
    client = ZincClient()
    
    with pytest.raises(exc.InvalidFileFormatError, match="is not a valid fileformat"):
        client._validate_filters(fileformat="pdf")
    
    with pytest.raises(exc.InvalidAvailabilityError, match="is not a valid availability"):
        client._validate_filters(fileformat="smi", availability="instantaneous")
        
    with pytest.raises(exc.InvalidBioactiveError, match="is not a valid bioactivity"):
        client._validate_filters(fileformat="smi", bioactive="hello")
    
    with pytest.raises(exc.InvalidBiogenicError, match="is not a valid biogenic"):
        client._validate_filters(fileformat="smi", biogenic="hello")
        
    with pytest.raises(exc.InvalidReactivityError, match="is not a valid reactivity"):
        client._validate_filters(fileformat="smi", reactivity="explosive")
    
def test_append_filters_to_url():
    
    client = ZincClient()
    catalog_url = client._catalog_url + "/chembl27/substances"
    
    full_url = client._append_filters_to_url(catalog_url, "smi", count=1000)
    assert full_url == "https://zinc.docking.org/catalogs/chembl27/substances.smi?count=1000"
    
    full_url = client._append_filters_to_url(catalog_url, "smi", count=1000, availability="for-sale")
    assert full_url == "https://zinc.docking.org/catalogs/chembl27/substances/subsets/for-sale.smi?count=1000"
    
    full_url = client._append_filters_to_url(catalog_url, "smi", count=1000, availability="for-sale", reactivity="anodyne")
    assert full_url == "https://zinc.docking.org/catalogs/chembl27/substances/subsets/for-sale+anodyne.smi?count=1000"
    
    full_url = client._append_filters_to_url(catalog_url, "smi", count=1000, bioactive="fda")
    assert full_url == "https://zinc.docking.org/catalogs/chembl27/substances/subsets/fda.smi?count=1000"
    
    full_url = client._append_filters_to_url(catalog_url, "smi", count=1000, biogenic="metabolites")
    assert full_url == "https://zinc.docking.org/catalogs/chembl27/substances/subsets/metabolites.smi?count=1000"
    
    substances_url = client._substances_url
    
    full_url = client._append_filters_to_url(substances_url, "smi", count=1000)
    assert full_url == "https://zinc.docking.org/substances.smi?count=1000"
    
    full_url = client._append_filters_to_url(catalog_url, "smi", count=1000, availability="for-sale", reactivity="anodyne")
    assert full_url == "https://zinc.docking.org/substances/subsets/for-sale+anodyne.smi?count=1000"
    

def test_catalog_url():
    pass

def test_urls_for_tranches_2d():
    pass

def test_predifined_subset_2d_tranches():
    pass

def test_molw_and_logp_2d_tranches():
    pass

def test_predifined_subset_3d_urls():
    pass

def test_molw_and_logp_3d_urls():
    pass

@pytest.mark.parametrize("value,lower", [
    (230, True),
    (484, False),
    (600, True)
])
def test_discretize_values(value, lower):
    
    bins = [200, 250, 300, 325, 350, 375, 400, 425, 450, 500, 550]
    new_value = ZincClient.discretize_values(value=value, bins=bins, name="Test", lower=lower)

    if value == 230:
        assert new_value == 200
    elif value == 484:
        assert new_value == 500
    else:
        assert new_value == 550
        

@pytest.mark.parametrize("subset,mol_weight,logp,format", [
    ("Drug-Like", None, None, "smi"),
    (None, (250, 350), (-1, 1), "smi"),
    (None, (365, 415), (1.5, 2.25), "smi"),
    ("Drug-Like", None, None, "sdf"),
    (None, (200, 300), (-1, 2), "sdf"),
])
def test_download_ZINC2D_smiles(subset, mol_weight, logp, format):

    url_list = get_zinc_urls( 
        subset=subset,
        mw_range=mol_weight, 
        logp_range=logp,
        file_format=format,
        )
    
    if format == "smi":
        base_url = "http://files.docking.org/2D/"
        if subset == "Drug-like":
            assert len(url_list) == 90 * 4 * 2 
            assert url_list[0] == base_url + "BA/BAAA.smi"
            assert url_list[-1] == base_url + "JJ/JJEB.smi"
        elif mol_weight == (250, 350):
            assert len(url_list) == 12 * 4 * 2
            assert url_list[0] == base_url + "BA/BAAA.smi"
            assert url_list[-1] == base_url + "EC/ECEB.smi"
        elif mol_weight == (365, 415):
            assert len(url_list) == 12 * 4 * 2
            assert url_list[0] == base_url + "EC/ECAA.smi"
            assert url_list[-1] == base_url + "HE/HEEB.smi"
    else:
        base_url = "http://files.docking.org/3D/"
        if subset == "Drug-like":
            assert len(url_list) == 19420
            assert url_list[0] == base_url + "JJ/EDRP/JJEDRP.xaa.sdf.gz"
            assert url_list[-1] == base_url + "AB/AAMM/ABAAMM.xaa.sdf.gz"
        elif mol_weight == (200, 300):
            assert len(url_list) == 3720
            assert url_list[0] == base_url + "AA/AAML/AAAAML.xaa.sdf.gz"
            assert url_list[-1] == base_url + "DC/EDRP/DCEDRP.xaa.sdf.gz"