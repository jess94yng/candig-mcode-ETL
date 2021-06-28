import pytest
import json
from os import path

f = open('mCodePacket.json')
data = json.load(f)


# test overall data type
def test_data_type():
    assert isinstance(data, list), 'data must be array'
    assert all(isinstance(item, dict) for item in data), 'data array elements must be json object'


# function to check if element is in a data structure
def is_in_data(element, item):
    if element in item:
        return True
    return False


# function to check if element is a certain datatype
def is_datatype(element, datatype):
    if isinstance(element, datatype):
        return True
    return False


# test for presence of all required mcode data elements
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(name = 'test_mcode_fields')
def test_mcode_fields(patients):
    assert is_in_data('id', patients), 'required mcode id field missing'
    assert is_in_data('meta_data', patients), 'required mcode meta data field missing'
    assert is_in_data('subject', patients), 'required mcode subject field missing'
    assert is_in_data('genomics_report', patients), 'required mcode genomics report field missing'
    assert is_in_data('cancer_condition', patients), 'required mcode cancer condition field missing'
    assert is_in_data('cancer_related_procedures', patients), 'required mcode cancer related procedures field missing'
    assert is_in_data('medication_statement', patients), 'required mcode medication statement field missing'
    assert is_in_data('tumor_marker', patients), 'required mcode tumor marker field missing'


# test for correct data types of mcode data elements
@pytest.mark.parametrize('mcode_data_element', [
    'meta_data',
    'id',
    'subject',
    'genomics_report',
    'cancer_condition',
    'cancer_related_procedures',
    'medication_statement',
    'tumor_marker'
])
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(depends=['test_mcode_fields'])
def test_mcode_format(patients, mcode_data_element):
    if mcode_data_element == 'id' or mcode_data_element == 'date_of_death':
        assert is_datatype(patients[mcode_data_element], str), 'wrong mcode field data type'
    elif mcode_data_element == 'medication_statement' or mcode_data_element == 'tumor_marker':
        assert is_datatype(patients[mcode_data_element], list), 'wrong mcode field data type'
    else: 
        assert is_datatype(patients[mcode_data_element], dict), 'wrong mcode field data type'


# test for presence of all required subject data elements
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(name = 'test_subject_required_fields')
def test_subject_required_fields(patients):
    assert is_in_data('subject', patients), 'required subject field missing'
    assert is_in_data('id', patients['subject']), 'required subject id field missing'
    assert is_in_data('code', patients['subject']['comorbid_condition']), 'required subject comorbid condition field missing'
    assert is_in_data('communication_language', patients['subject']['extra_properties']), 'required subject communication language field missing'
    assert is_in_data('name', patients['subject']['extra_properties']), 'required subject name field missing'
    assert is_in_data('administrative_gender', patients['subject']['extra_properties']), 'required subject administrative gender field missing'

# test for correct data types of subject data elements
@pytest.mark.parametrize('subject_data_element', [
    'id',
    'date_of_birth',
    'sex',
    'taxonomy',
    'race',
    'ethnicity',
    'comorbid_condition',
    'ecog_performance_status',
    'karnofsky'
    'extra_properties'
])
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(depends=['test_subject_required_fields'])
def test_subject_format(patients, subject_data_element):
    if is_in_data(subject_data_element, patients['subject']):
        if subject_data_element == 'id' or subject_data_element == 'date_of_birth' or subject_data_element == 'sex' or subject_data_element == 'race' or subject_data_element == 'ethnicity':
            assert is_datatype(patients['subject'][subject_data_element], str), 'wrong subject field data type'
        else:
            assert is_datatype(patients['subject'][subject_data_element], dict), 'wrong subject field data type'


# test for presence of all required genomics report data elements
@pytest.mark.parametrize('genomics_data_element', [
    'id',
    'code',
    'issued'
])
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(name = 'test_genomics_report_required_fields')
def test_genomics_report_required_fields(patients, genomics_data_element):
    assert is_in_data('genomics_report', patients), 'required genomics report field missing'
    assert is_in_data(genomics_data_element, patients['genomics_report']), '1 or more required genomics report fields missing'
    if is_in_data('genetic_specimen', patients['genomics_report']):
        assert isinstance(patients['genomics_report']['genetic_specimen'], list), 'genetic specimen field must be array'
        for specimen_item in patients['genomics_report']['genetic_specimen']:
            assert is_in_data('specimen_type', specimen_item), 'required genomics report genetic specimen type field missing'


# test for correct data types of genomics report data elements
@pytest.mark.parametrize('report_data_element', [
    'id',
    'code',
    'performing_organization_name',
    'issued',
    'genetic_specimen',
    'genetic_variant',
    'genomic_region_studied',
    'extra_properties'
])
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(depends=['test_genomics_report_required_fields'])
def test_genomics_report_format(patients, report_data_element):
    if is_in_data(report_data_element, patients['genomics_report']):
        if report_data_element == 'id' or report_data_element == 'performing_organization_name' or report_data_element == 'issued':
            assert is_datatype(patients['genomics_report'][report_data_element], str), 'wrong genomics report field data type'
        elif report_data_element == 'genetic_specimen':
            assert is_datatype(patients['genomics_report'][report_data_element], list), 'wrong genomics report field data type'
        else: 
            assert is_datatype(patients['genomics_report'][report_data_element], dict), 'wrong genomics report field data type'


# test for presence of all required cancer condition data elements
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(name = 'test_cancer_condition_required_fields')
def test_cancer_condition_required_fields(patients):
    assert is_in_data('cancer_condition', patients), 'required  cancer condition field missing'
    assert is_in_data('code', patients['cancer_condition']), 'required cancer condition code field missing'


# test for correct data types of cancer condition data elements
@pytest.mark.parametrize('condition_data_element', [
    'id',
    'condition_type',
    'body_site',
    'code',
    'date_of_diagnosis',
    'histology_morphology_behavior',
    'extra_properties'
])
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(depends=['test_cancer_condition_required_fields'])
def test_cancer_condition_format(patients, condition_data_element):
    if is_in_data(condition_data_element, patients['cancer_condition']):
        if condition_data_element == 'id' or condition_data_element == 'condition_type' or condition_data_element == 'date_of_diagnosis':
            assert is_datatype(patients['cancer_condition'][condition_data_element], str), 'wrong cancer condition field data type'
        else: 
            assert is_datatype(patients['cancer_condition'][condition_data_element], dict), 'wrong cancer condition field data type'


# test for presence of all required cancer related procedure data elements
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(name = 'test_cancer_related_procedures_required_fields')
def test_cancer_related_procedures_required_fields(patients):
    assert is_in_data('cancer_related_procedures', patients), 'required cancer related procedures field missing'
    assert is_in_data('code',patients['cancer_related_procedures']), 'required cancer related procedures code field missing'


# test for correct data types of cancer related procedure data elements
@pytest.mark.parametrize('procedure_data_element', [
    'id',
    'procedure_type',
    'code',
    'body_site',
    'extra_properties'
])
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(depends=['test_cancer_related_procedures_required_fields'])
def test_cancer_procedure_format(patients, procedure_data_element):
    if is_in_data(procedure_data_element, patients['cancer_related_procedures']):
        if procedure_data_element == 'id' or procedure_data_element == 'procedure_type':
            assert is_datatype(patients['cancer_related_procedures'][procedure_data_element], str), 'wrong cancer related procedures field data type'
        elif procedure_data_element == 'body_site':
            assert is_datatype(patients['cancer_related_procedures'][procedure_data_element], list), 'wrong cancer related procedures field data type'
        else: 
            assert is_datatype(patients['cancer_related_procedures'][procedure_data_element], dict), 'wrong cancer related procedures field data type'


# test for presence of all required medication statement required fields
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(name = 'test_medication_statement_required_fields')
def test_medication_statement_required_fields(patients):
    assert is_in_data('medication_statement', patients), 'required medication statement field missing'


# test for correct data types of medication statement data elements
@pytest.mark.parametrize('medication_data_element', [
    'id',
    'medication_code',
    'extra_properties'
])
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(depends=['test_medication_statement_required_fields'])
def test_medication_format(patients, medication_data_element):
    if is_in_data(medication_data_element, patients['medication_statement']):
        if medication_data_element == 'id':
            assert is_datatype(patients['medication_statement'][medication_data_element], str), 'wrong medication id field data type'
        else:
            assert is_datatype(patients['medication_statement'][medication_data_element], dict), 'wrong medication code field data type'


# test for presence of all required tumor marker data elements
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(name = 'test_tumor_marker_required_fields')
def test_tumor_marker_required_fields(patients):
    assert is_in_data('tumor_marker', patients), 'required tumor marker field missing'
    assert isinstance(patients['tumor_marker'], list), 'tumor_marker field must be array'
    for tumor_item in patients['tumor_marker']:
        assert is_in_data('tumor_marker_code', tumor_item), 'required tumor marker code field missing'  


# test for correct data types of tumor marker data elements
@pytest.mark.parametrize('tumor_data_element', [
    'id',
    'individual',
    'tumor_marker_code',
    'tumor_marker_data_value',
    'extra_properties'
])
@pytest.mark.parametrize('patients', data)
@pytest.mark.dependency(depends=['test_tumor_marker_required_fields'])
def test_tumor_marker_format(patients, tumor_data_element):
    for tumor_item in patients['tumor_marker']:
        if is_in_data(tumor_data_element, tumor_item):
            if tumor_data_element == 'id' or tumor_data_element == 'individual':
                assert is_datatype(tumor_item[tumor_data_element], str), 'wrong tumor marker field data type'
            else: 
                assert is_datatype(tumor_item[tumor_data_element], dict), 'wrong tumor marker field data type'




