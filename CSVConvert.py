#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import json
from datetime import date
from MappingDict import race_mapping_dict, valid_mcode, mcode_dict, lab_results, vital_signs, possible_multi_val
from Ontologies import ontologies
from SchemaBase import return_mcodepacket, return_subject, return_genomics_report, return_genetic_region_studied, return_cancer_condition, return_cancer_related_procedures, return_tumor_marker
from IntegrateAPI import get_json, get_url, subtree_dict, hgnc_api, hgvs_api_esearch, hgvs_api_esummary
from SchemaGeneration import generate_genetic_specimen, generate_body_site, generate_medication, generate_tumor_marker
import argparse

# Set up argparse for command line inputs
# Data_path: path to dataset for ingestion/testing
# BioPortal_API_Key: found on BioPortal account, used to access BioPortal API
# Email: organizational email used to access NCBI clinvar API
parser = argparse.ArgumentParser()
parser.add_argument('--Data_path', type=str, required=True)
parser.add_argument('--BioPortal_API_Key', type=str, required=True)
parser.add_argument('--Email', type=str, required=True)
args = parser.parse_args()
dataset = args.Data_path
api_key = args.BioPortal_API_Key
email = args.Email


# test if value is null
def not_null(cell):
    if cell == 'nan':
        return False
    else:
        return True


# map race by sorting through columns of different race codes
def race_mapping(input_dataframe, ind, cur_row, mapping_list):
    for list_item in mcode_dict['race']:    #for each possible race
        if list_item in input_dataframe:    #check if in current dataframe
            if cur_row[list_item] != '0':       #check if value marked as '1' or possibly name of race written
                if cur_row[list_item] == '1':
                    mapping_list.append(('race', race_mapping_dict[list_item]))
                elif pd.notna(cur_row[list_item]):
                    mapping_list.append(('race', cur_row[list_item]))


#eliminate duplicate identifiers and group values from the same data element into a list in the dataframe
def eliminate_duplicate(input_dataframe):
    input_dataframe = input_dataframe.dropna(how='all')
    input_dataframe = input_dataframe.applymap(str)
    
    for index, row in input_dataframe.iterrows():   #for each element in the list of data elements that are expected to be in array form
        for multi_val_item in possible_multi_val:
            if multi_val_item in input_dataframe:   #check if data element exist in current dataframe
                row[multi_val_item] = row[multi_val_item].split(',')    #change string to appendable list

        if (input_dataframe.duplicated(subset=['identifier'], keep=False)[index]): #if found rows with duplicated identifier
            for multi_val_item in possible_multi_val:   #for each element in list of data elements that are expected to be in array form
                if multi_val_item in input_dataframe: #if element exists in current dataframe
                    for index1, row1 in input_dataframe.iterrows(): #find following repeated identifiers and append data elements to the list if they are not null
                        if index1 != index and row1['identifier'] == row['identifier'] and row1[multi_val_item] != 'nan':
                            row[multi_val_item].append(row1[multi_val_item])
    input_dataframe = input_dataframe.drop_duplicates(subset='identifier')  #drop all identifier duplicates after the first one
    return input_dataframe


# Process csv data, use mappings to return data frame with mCode data elements
def mcode_mapping(input_dataframe):

    mcode_csv_list = []

    for index, row in input_dataframe.iterrows():   #iterate through rows of dataframe
        data_element_list = []
        dict_item = {}
        for item in valid_mcode:    #iterate through data element
            if item == 'race':  #use race_mapping
                race_mapping(input_dataframe, index, row, data_element_list)
            elif mcode_dict[item] in input_dataframe:   #add a tuple with the data element and the corresponding value
                data_element_list.append((item, row[mcode_dict[item]]))

        for element in data_element_list:
            dict_item[element[0]] = element[1]  #turn list of tuples into dictionary

        mcode_csv_list.append(dict_item)

    df1 = pd.DataFrame(mcode_csv_list)  #transform dictionary into new mapped dataframe
    df1 = eliminate_duplicate(df1)  #eliminate duplicate identifiers
    return df1


# Combine dataframes from multiple sheets, delete any duplicate patients by merging data
def combine_sheets(list_dataframe):
    for ind in range(1, len(list_dataframe)):   #merge all dataframes in list
        list_dataframe[0] = list_dataframe[0].merge(
            list_dataframe[ind], how='outer')
    output_df = list_dataframe[0]
    for index, row in output_df.iterrows():     #compare across the newly merged sheets and eliminate repeated identifiers
        if(output_df.duplicated(subset=['identifier'], keep=False)[index]):
            for column in output_df:    #check all data elements
                notna_row = False
                if isinstance(row[column], list) or pd.notna(row[column]):  #if the field is a list or is a non-null single value, make not-null boolean true
                    notna_row = True
                for index1, row1 in output_df.iterrows(): #find repeated identifiers, replace identifier's first occurance field values if not null and first occurance field value is null
                    notna_row1 = False
                    if isinstance(row1[column], list) or pd.notna(row1[column]):
                        notna_row1 = True
                    if index1 > index and row1['identifier'] == row['identifier'] and notna_row == False and notna_row1 == True:
                        row[column] = row1[column]
    output_df = output_df.drop_duplicates(subset='identifier')
    return output_df


# Organize dataframe into mCode model
def mcodepacket_build(dataframe_item, counter_var):
    mcodepacket = return_mcodepacket(dataframe_item['identifier'])

    mcodepacket['subject'] = subject_build(dataframe_item)

    mcodepacket['genomics_report'] = genomics_report_build(dataframe_item, counter_var)

    mcodepacket['cancer_condition'] = cancer_condition_build(dataframe_item, counter_var)

    mcodepacket['cancer_related_procedures'] = cancer_related_procedures_build(dataframe_item, counter_var)

    mcodepacket['medication_statement'] = medication_statement_build(dataframe_item, counter_var)

    if dataframe_item.get('date_of_death', None) and not_null(dataframe_item['date_of_death']):
        mcodepacket['date_of_death'] = dataframe_item['date_of_death']

    mcodepacket['tumor_marker'] = tumor_marker_build(dataframe_item, counter_var)

    if dataframe_item.get('cancer_disease_status', None) and not_null(dataframe_item['cancer_disease_status']):
        if dataframe_item['cancer_disease_status'] in ontologies['cancer_disease_status']:
            cancer_disease_status = ontologies['cancer_disease_status'][dataframe_item['cancer_disease_status'].lower()]
            mcodepacket['cancer_disease_status'] = cancer_disease_status

    return (mcodepacket)


def subject_build(dataframe_item):
    subject = return_subject(dataframe_item['identifier'])

    if dataframe_item.get('date_of_birth', None) and not_null(dataframe_item['date_of_birth']):
        subject['date_of_birth'] = dataframe_item['date_of_birth']
    if dataframe_item.get('birth_sex', None) and not_null(dataframe_item['birth_sex']):
        subject['sex'] = ontologies['sex'][dataframe_item['birth_sex'].lower()]
    if dataframe_item.get('race', None) and not_null(dataframe_item['race']):
        subject['race'] = ontologies['race'][dataframe_item['race'].lower()]
    if dataframe_item.get('ethnicity', None) and not_null(dataframe_item['ethnicity']):
        subject['ethnicity'] = ontologies['ethnicity'][dataframe_item['ethnicity'].lower()]
    subject['comorbid_condition'] = {
        'code': {
            'id': 'SNOMED:103329007',
            'label': 'Not available'
        }
    }
    if dataframe_item.get('ecog_performance_status', None) and not_null(dataframe_item['ecog_performance_status']):
        subject['ecog_performance_status'] = ontologies['ecog_performance_status'][dataframe_item['ecog_performance_status']]
    subject['karnofsky_performance_status'] = {
        'id': 'KARNOFSKY: Not available',
        'label': 'Not available'
    }

    for vital_item in vital_signs:  #add all vital signs into subject extra properties
        if dataframe_item.get(vital_item, None) and not_null(dataframe_item[vital_item]):
            subject['extra_properties'][vital_item] = dataframe_item[vital_item]
    
    return subject


def genomics_report_build(dataframe_item, counter_var):
    genomics_report = return_genomics_report(counter_var)

    if dataframe_item.get('performing_organization_name', None) and not_null(dataframe_item['performing_organization_name']):
        genomics_report['performing_organization_name'] = dataframe_item['performing_organization_name']
    if dataframe_item.get('issued', None) and not_null(dataframe_item['issued']):
        genomics_report['issued'] = dataframe_item['issued']
    else:
        genomics_report['issued'] = 'Unknown'

    genetic_specimen_build(dataframe_item, genomics_report, counter_var)

    genetic_variant_build(dataframe_item, genomics_report, counter_var)

    genetic_region_studied = return_genetic_region_studied(counter_var)
    gene_region_add = False
    gene_region_add = genetic_mutation_build(dataframe_item, genetic_region_studied, gene_region_add)
    gene_region_add = gene_studied_build(dataframe_item, genetic_region_studied, gene_region_add)

    if gene_region_add:
        genomics_report['genetic_region_studied'] = genetic_region_studied

    return genomics_report


def genetic_specimen_build(dataframe_item, genomics_report, counter_var):
    genetic_specimen = []
    if dataframe_item.get('genetic_specimen', None) and not_null(dataframe_item['genetic_specimen']):
        temp_counter1 = 0
        for i in range(len(dataframe_item['genetic_specimen'])):    #for all array form data elements, iterate through list and append
            genetic_specimen.append({
                'id': str(counter_var) + '-' + str(temp_counter1),
                'collection_body': {
                    'id': generate_genetic_specimen(dataframe_item['collection_body_site'][i], api_key)[0],
                    'label': generate_genetic_specimen(dataframe_item['collection_body_site'][i], api_key)[1]
                }
            })
            temp_counter1 += 1

            if dataframe_item.get('specimen_type', None) and dataframe_item['specimen_type'] in ontologies['specimen_type']:
                genetic_specimen[-1]['specimen_type'] = ontologies['specimen_type'][dataframe_item['specimen_type'][i].lower()]
            else:
                genetic_specimen[-1]['specimen_type'] = {
                    'id': 'HL7:...',
                    'label': 'No suggested value'
                }
        genomics_report['genetic_specimen'] = genetic_specimen


def genetic_variant_build(dataframe_item, genomics_report, counter_var):
    if dataframe_item.get('genetic_variant', None) and not_null(dataframe_item['genetic_variant']):
        if dataframe_item['variant_data_value'] in ontologies['variant_data_value']:
            genomics_report['genetic_variant'] = {
                'id': str(counter_var),
                'data_value': ontologies['variant_data_value'][dataframe_item['variant_data_value'].lower()]
            }
        else:
            genomics_report['genetic_variant'] = {
                'id': str(counter_var),
                'data_value': {
                    'id': 'LOINC:LA4489-6',
                    'label': 'Unknown'
                }
            }


def genetic_mutation_build(dataframe_item, genetic_region_studied, genetic_region_add):
    if dataframe_item.get('gene_mutation', None) and not_null(dataframe_item['gene_mutation']):
        gene_region_add = True
        gene_mutation = []
        for i in range(len(dataframe_item['gene_mutation'])):
            result = hgvs_api_esearch(dataframe_item['gene_mutation'][i], email)
            if len(result['IdList']) > 0:
                gene_mutation.append({
                    'id': 'HGVS:' + result['IdList'][0],
                    'label': hgvs_api_esummary(result['IdList'][0], email)
                })
        genetic_region_studied['gene_mutation'] = gene_mutation
    return gene_region_add


def gene_studied_build(dataframe_item, genetic_region_studied, gene_region_add):
    if dataframe_item.get('gene_studied', None) and not_null(dataframe_item['gene_studied']):
        gene_region_add = True
        gene_studied = []
        for i in range(len(dataframe_item['gene_studied'])):
            result = hgnc_api(dataframe_item['gene_studied'][i])
            if result['response']['numFound'] >= 1:
                gene_studied.append({
                    'id': result['response']['docs'][0]['hgnc_id'],
                    'label': result['response']['docs'][0]['symbol']
                })
        genetic_region_studied['gene_studied'] = gene_studied
    return gene_region_add


def cancer_condition_build(dataframe_item, counter_var):
    cancer_condition = return_cancer_condition(counter_var)

    if dataframe_item.get('date_of_diagnosis', None) and not_null(dataframe_item['date_of_diagnosis']):
        cancer_condition['date_of_diagnosis'] = dataframe_item['date_of_diagnosis']
    if dataframe_item.get('history_morphology_behavior', None) and not_null(dataframe_item['history_morphology_behavior']):
        behavior_term = dataframe_item['history_morphology_behavior'].replace(
            ' ', '+')
        need_break = False
        for behavior_subtree in subtree_dict['history_morphology_behaviour']:   #several potential subtrees for ontology search
            term_colleciion = get_json(get_url('SNOMEDCT', behavior_term, False, subtree_dict['history_morphology_behaviour'][behavior_subtree]), api_key)['collection']
            if len(term_colleciion) > 0:    #only use ontologies when 1 or more concepts are found
                for list_item in term_colleciion:
                    if 'notation' in list_item:
                        identification = 'SNOMED:' + list_item['notation']
                        label = list_item['prefLabel']
                        need_break = True
                        break
            if need_break:
                break
            if behavior_subtree == 'subtreeThree':  #if on last subtree and still no results, return unknown
                identification = 'SNOMED:261665006'
                label = 'Unknown'
        cancer_condition['history_morphological_behavior'] = {
            'id': identification,
            'label': label
        }
    
    return cancer_condition

def cancer_related_procedures_build(dataframe_item, counter_var):
    cancer_related_procedures = return_cancer_related_procedures(counter_var)

    if dataframe_item.get('body_site', None) and not_null(dataframe_item['body_site']):
        bodySite = []
        for i in range(len(dataframe_item['body_site'])):
            bodySite.append({
                'id': generate_body_site(dataframe_item['body_site'][i], api_key)[0],
                'label': generate_body_site(dataframe_item['body_site'][i], api_key)[1]
            })

        cancer_related_procedures['body_site'] = bodySite
    
    return cancer_related_procedures

def medication_statement_build(dataframe_item, counter_var):
    medication_statement = []
    if dataframe_item.get('medication', None) and not_null(dataframe_item['medication']):
        temp_counter2 = 0
        for medicine in dataframe_item['medication']:
            medication_statement.append({
                'id': str(counter_var) + '-' + str(temp_counter2),
                'medication_code': {
                    'id': generate_medication(medicine, api_key)[0],
                    'label': generate_medication(medicine, api_key)[1]
                }
            })
            temp_counter2 += 1
    else:
        medication_statement.append({
            'id': str(counter_var),
            'medication_code': {
                'id': 'RxNorm:Unknown',
                'label': 'Unknown'
            }
        })
    
    return medication_statement

def tumor_marker_build(dataframe_item, counter_var):
    tumor_marker = []
    tumor_extra = {}
    extra = False
    for extra_item in lab_results:  #add lab_result data elements to tumor marker extra properties
        if dataframe_item.get(extra_item, None) and not_null(dataframe_item[extra_item]):
            tumor_extra[extra_item] = dataframe_item[extra_item]
            extra = True

    if dataframe_item.get('tumor_marker_data_value', None) and not_null(dataframe_item['tumor_marker_data_value']):
        temp_counter3 = 0
        for data in dataframe_item['tumor_marker_data_value']:
            tumor_marker.append(return_tumor_marker(
                counter_var, temp_counter3, dataframe_item['identifier'], True))
            tumor_marker = generate_tumor_marker(data, tumor_marker, api_key)

            if extra:
                tumor_marker[-1]['extra_properties'] = tumor_extra
            temp_counter3 += 1

    else:
        tumor_marker.append(return_tumor_marker(
            counter_var, 0, dataframe_item['identifier'], False))

        if extra:
            tumor_marker[-1]['extra_properties'] = tumor_extra

    return tumor_marker


def main():
    genomics_id = 1000  #set up unique identifiers for data elements in mcode data
    df_og = pd.read_excel(dataset, sheet_name=None, dtype=str)  #read sheet from given data pathway
    df_list = []
    for page in df_og:
        df_list.append(mcode_mapping(df_og[page]))  #append all processed mcode dataframes to a list

    mcode_dataframe = combine_sheets(df_list)   #combine sheets
    mcode_dataframe = mcode_dataframe.fillna('nan')
    mcode_list = []
    for index, row in mcode_dataframe.iterrows():   #build mcode data model from combined dataframe
        mcode_list.append(mcodepacket_build(row, genomics_id))
        genomics_id += 1
    print(mcode_list)
    with open('mCodePacket.json', 'w') as f:    #write to json file for ingestion
        json.dump(mcode_list, f, indent=4)


if __name__ == '__main__':
    main()
