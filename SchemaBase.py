#!/usr/bin/env python
# coding: utf-8

# In[1]:
# Create dictionaries with base required schema structures to use in mCode data organization
from datetime import datetime
from Ontologies import Resources

def returnMCodePacket(identity):
    return {
        'id': identity,
        "meta_data": {
            "created": str(datetime.now()),
            "mcodepacket_schema_version": "1.0.0-RC3",
            "created_by": "Ksenia Zaytseva",
            "submitted_by": "Ksenia Zaytseva",
            "resources": [
                Resources["NCBITaxon"],
                Resources["SNOMED"],
                Resources["RxNorm"],
                Resources["LOINC"],
                Resources["HL7"],
                Resources["HGNC"],
                Resources["HGVS"],
            ]
        }
    }
def returnSubject(identity):
    return {
        'id': identity,
        'taxonomy': {
            'id': 'NCBITaxon:9606',
            'label': 'Homo sapiens'
        },
        'comorbid_condition':{
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
def returnGenomicsReport(counter):
    return {
        'id': str(counter),
        'code': {
            'id': 'LOINC:LA7338-2',
            'label': 'Not available'
        }
    }
def returnGeneticRegionStudied(counter):
    return {
        'id': str(counter),
    }
def returnCancerCondition(counter):
    return {
        'id': str(counter),
        'condition_type': 'primary',
        'code': {
            'id': 'SNOMED:103329007',
            'label': 'Not available'
        }
    }

def returnCancerRelatedProcedures(counter):
    return {
        'id': str(counter),
        'procedure_type': 'radiation',
    #         code for surgical procedure required, but schema only allows one of radiation/surgical procedures
        'code': {
            'id': 'SNOMED:103329007',
            'label': 'Not available'
        }
    }

def returnTumorMarker(bigCounter, smallCounter, identity, useSmall):
    if useSmall:
        return {
            'id': str(bigCounter)+'-'+str(smallCounter),
            'individual': identity,
            'tumor_marker_code': {
                'id': 'LOINC:LA7338-2',
                'label': 'Not available'
            }
        }
    else:
        return {
            'id': str(bigCounter),
            'individual': identity,
            'tumor_marker_code': {
                'id': 'LOINC:LA7338-2',
                'label': 'Not available'
            }
        }

