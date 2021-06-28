import pytest
import json
from CSVConvert import mcode_mapping, race_mapping, eliminate_duplicate
from MappingDict import mcode_dict, valid_mcode, possible_multi_val
import pandas as pd


df = pd.read_excel('data/pytest_data_v1.xlsx', dtype=str)  #read sheet from given data pathway
print(df)
df_drop_na = df.fillna('nan')
df_modified = df_drop_na.drop_duplicates(subset='Subject')
# df_modified = df_modified.fillna('nan')
df_mapped = mcode_mapping(df)

print (df_mapped['medication'])
print(df['THER_TX_NAME'])
counter = 0
for index, row in df_mapped.iterrows():
    print(type(row['gene_mutation']))

def count_array_length(mcode_item):
    counter = 0
    for index, row in df_drop_na.iterrows():
        counter += len(row[mcode_dict[mcode_item]].split(','))
    return counter

def count_nan(mcode_item):
    counter = 0
    for index, row in df_drop_na.iterrows():
        counter += row[mcode_dict[mcode_item]].count('nan')
    return counter

def count_mapped_array(mcode_item):
    counter = 0
    for index, row in df_mapped.iterrows():
        counter += len(row[mcode_item])
    return counter


def test_duplicate():
    assert any(df_mapped.duplicated(subset=['identifier'], keep=False)) == False

@pytest.mark.parametrize('mcode', valid_mcode)
def test_mapping(mcode):
    if mcode in df_mapped and mcode not in possible_multi_val and mcode!= 'race':
        check_mapping = df_mapped[mcode].eq(df_modified[mcode_dict[mcode]])
        assert check_mapping.all() == True, 'mapping incorrect'
    elif mcode in df_mapped and mcode in possible_multi_val:
        assert count_mapped_array(mcode) <= count_array_length(mcode), 'array fields have more values than expected'
        assert count_mapped_array(mcode) >= count_array_length(mcode) - count_nan(mcode), 'array fields have less values than expected'

@pytest.mark.parametrize('multival', possible_multi_val)
def test_array(multival):
    if multival in df_mapped:
        for index, row in df_mapped.iterrows():
            assert isinstance(row[multival], list), 'array field should be list datatype'
        





