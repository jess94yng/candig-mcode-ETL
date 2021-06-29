import pytest
import json
from CSVConvert import mcode_mapping, race_mapping, eliminate_duplicate, combine_sheets
from MappingDict import mcode_dict, valid_mcode, possible_multi_val
import pandas as pd

df = pd.read_excel('data/pytest_data_v2.xlsx', sheet_name=None, dtype=str, engine='openpyxl')  #read sheet from given data pathway
df_list = []
for page in df:
    df_list.append(mcode_mapping(df[page]))  #append all processed mcode dataframes to a list

mcode_dataframe = combine_sheets(df_list)   #combine sheets
mcode_dataframe = mcode_dataframe.fillna('nan')


# test for duplicated identifiers
def test_combine():
    assert any(mcode_dataframe.duplicated(subset=['identifier'], keep=False)) == False