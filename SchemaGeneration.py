#!/usr/bin/env python
# coding: utf-8
# Formulate identifications and labels for fields using codeable concepts and various ontololgies

from IntegrateAPI import get_json, get_url, subtree_dict, hgnc_api, hgvs_api_esearch, hgvs_api_esummary


# build genetic specimen schema with ontology search
def generate_genetic_specimen(term_str, api_key):
    str_list = term_str.split(' ')  #extract body site from procedure name
    need_break = False
    for term in str_list:
        term_collection = get_json(get_url('SNOMEDCT', term, True, subtree_dict['collection_body']), api_key)['collection']
        if len(term_collection) > 0:    #check if there are at least 1 concept
            for list_item in term_collection:  
                if 'notation' in list_item:
                    synonym_list = list_item['synonym'] 
                    for synonym in synonym_list:    #check if there are synonyms that match up with the term
                        if term.lower() == synonym.lower():
                            identification = 'SNOMED:' + list_item['notation']
                            label = list_item['prefLabel']
                            need_break = True
                            break
                if need_break:
                    break
        if need_break:
            break
        if term == str_list[-1]:    #else return unknown
            identification = 'SNOMED: 261665006'
            label = 'Unknown'

    id_label = [identification, label]
    return id_label


#find body site with ontology search
def generate_body_site(term_str, api_key): 
    if '/' in term_str:     #in case multiple body sites are present, break into individual terms
        body_part = term_str.split('/')[0]
    else:
        body_part = term_str
    body_part_new = body_part.replace(' ', '+')     #fill spaces with plus signs
    term_collection = get_json(get_url(
        'SNOMEDCT', body_part_new, False, subtree_dict['collection_body']), api_key)['collection']
    id_label = generate_id_label(
        term_collection, 'SNOMEDCT', 'SNOMED: 261665006')
    return id_label


#find medicine id and label based on ontology search
def generate_medication(med_item, api_key):
    med_string = med_item.replace(' ', '+')
    term_collection = get_json(get_url('RxNorm', med_string, False, ''), api_key)[
        'collection']
    id_label = generate_id_label(term_collection, 'RxNorm', 'RxNorm:Unknown')
    return id_label


#find tumor marker, used ontologies if needed
def generate_tumor_marker(tumor_item, tumor_list, api_key):
    if tumor_item.replace('.', '1').isdigit():
        tumor_list[-1]['tumor_marker_data_value'] = {
            'value': {
                'value': float(tumor_item),
                'comparator': '=',
            }
        }
    elif '/' in tumor_item:
        num = tumor_item.split('/')
        tumor_list[-1]['tumor_marker_data_value'] = {
            'value': {
                'numerator': {
                    'value': float(num[0]),
                    'comparator': '='
                },
                'denominator': {
                    'value': float(num[1]),
                    'comparator': '='
                }
            }
        }
    else:
        data_str = tumor_item.replace(' ', '+')
        term_collection = get_json(get_url('loinc', data_str, False, ''), api_key)[
            'collection']
        id_label = generate_id_label(
            term_collection, 'LOINC', 'LOINC:LA4489-6')
        tumor_list[-1]['tumor_marker_data_value'] = {
            'value': {
                'id': id_label[0],
                'label': id_label[1]
            }
        }
    return tumor_list


#generate id and label
def generate_id_label(collection, db_name, unknown_id):
    if len(collection) > 0:
        for list_item in collection:
            if 'notation' in list_item:
                identification = db_name + ':' + list_item['notation']
                label = list_item['prefLabel']
                break
            if list_item == collection[-1]:
                identification = unknown_id
                label = 'Unknown'
    else:
        identification = unknown_id
        label = 'Unknown'

    id_label = [identification, label]
    return id_label
