#!/usr/bin/env python
# coding: utf-8

# Create dictionaries with base required schema structures to use in mCode data organization
from datetime import datetime
from Ontologies import resources


#build mcodepacket base
def return_mcodepacket(identity):
    return {
        'id': identity,
        'meta_data': {
            'created': str(datetime.now()),
            'mcodepacket_schema_version': '1.0.0-RC3',
            'created_by': 'Ksenia Zaytseva',
            'submitted_by': 'Ksenia Zaytseva',
            'resources': [
                resources['NCBITaxon'],
                resources['SNOMED'],
                resources['RxNorm'],
                resources['LOINC'],
                resources['HL7'],
                resources['HGNC'],
                resources['HGVS'],
            ]
        }
    }


#build subject base
def return_subject(identity):
    return {
        'id': identity,
        'taxonomy': {
            'id': 'NCBITaxon:9606',
            'label': 'Homo sapiens'
        },
        'comorbid_condition': {
            'code': {
                'id': 'Not available',
                'label': 'Not available'
            }
        },
        'extra_properties': {
            'communication_language': 'Not available',
            'administrative_gender': 'Not available',
            'name': 'Not available'
        }
    }


#build genomics report base
def return_genomics_report(counter):
    return {
        'id': str(counter),
        'code': {
            'id': 'LOINC:LA7338-2',
            'label': 'Not available'
        }
    }


#build genetic region studied base
def return_genetic_region_studied(counter):
    return {
        'id': str(counter),
    }


#build cancer condition base
def return_cancer_condition(counter):
    return {
        'id': str(counter),
        'condition_type': 'primary',
        'code': {
            'id': 'SNOMED:103329007',
            'label': 'Not available'
        }
    }


#build cancer related procedures base
def return_cancer_related_procedures(counter):
    return {
        'id': str(counter),
        'procedure_type': 'radiation',
        #         code for surgical procedure required, but schema only allows one of radiation/surgical procedures
        'code': {
            'id': 'SNOMED:103329007',
            'label': 'Not available'
        }
    }


#build tumor marker base
def return_tumor_marker(big_counter, small_counter, identity, use_small):
    if use_small:
        return {
            'id': str(big_counter)+'-'+str(small_counter),
            'individual': identity,
            'tumor_marker_code': {
                'id': 'LOINC:LA7338-2',
                'label': 'Not available'
            }
        }
    else:
        return {
            'id': str(big_counter),
            'individual': identity,
            'tumor_marker_code': {
                'id': 'LOINC:LA7338-2',
                'label': 'Not available'
            }
        }
