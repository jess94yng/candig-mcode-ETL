# medidata_mCode_ETL
Convert medidata rave data to mCode model for Katsu ingestion with Python

## Set-up & Installation
Prerequisites: 
- [Python 3.6+](https://www.python.org/)
- [pip](https://github.com/pypa/pip/)

Dependencies
- [Pandas](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html)
  - `pip install pandas`
- [Biopython](https://biopython.org/wiki/Download)
  - `pip install biopython`
- [Httplib2](https://pypi.org/project/httplib2/)
  - `pip install httplib2`
- [xlrd](https://pypi.org/project/xlrd/)
  - `pip install xlrd`
- [openpyxl](https://pypi.org/project/openpyxl/)
  - `pip install openpyxl`

## Running from command line
`$ python CSVConvert.py --Data_path DATA_PATH --BioPortal_API_Key BIOPORTAL_API_KEY --Email EMAIL`

--Data_path: path to dataset to be converted to mcode data model

--BioPortal_API_Key: BioPortal API key found in BioPortal personal account settings

--Email: organizational email used to access NCBI clinvar API

## Testing
Continuous Integration is implemented through Pytest and Travis CI which runs when git pushes occur. 

To run tests manually, enter from command line `$ pytest`
