#!/usr/bin/env python
# coding: utf-8

# Local ontologies
resources = {
    'NCBITaxon': {
        'name': 'NCBI Taxonomy OBO Edition',
        'version': '2020-11-02',
        'namespace_prefix': 'NCBITaxon',
        'id': 'NCBITaxon:2020-11-02',
        'iri_prefix': 'https://bioportal.bioontology.org/ontologies/NCBITAXON',
        'url': 'https://bioportal.bioontology.org/ontologies/NCBITAXON'
    },
    'SNOMED': {
        'name': 'SNOMED Clinical Terms',
        'version': '2021-04-03',
        'namespace_prefix': 'SNOMED',
        'id': 'SNOMED:2021-04-03',
        'iri_prefix': 'https://bioportal.bioontology.org/ontologies/SNOMEDCT',
        'url': 'https://bioportal.bioontology.org/ontologies/SNOMEDCT'
    },
    'RxNorm': {
        'name': 'RxNORM',
        'version': '2020-05-20',
        'namespace_prefix': 'RxNorm',
        'id': 'RxNorm:2020-05-20',
        'iri_prefix': 'https://bioportal.bioontology.org/ontologies/RXNORM',
        'url': 'https://bioportal.bioontology.org/ontologies/RXNORM'
    },
    'LOINC': {
        'name': 'Logical Observation Identifier Names and Codes',
        'version': '2021-04-03',
        'namespace_prefix': 'LOINC',
        'id': 'LOINC:2021-04-03',
        'iri_prefix': 'https://bioportal.bioontology.org/ontologies/LOINC',
        'url': 'https://bioportal.bioontology.org/ontologies/LOINC'
    },
    'HL7': {
        'name': 'Health Level Seven Terminologies',
        'version': '2021-03-07',
        'namespace_prefix': 'HL7',
        'id': 'HL7:2021-03-07',
        'iri_prefix': 'https://terminology.hl7.org/2.1.0/',
        'url': 'https://terminology.hl7.org/2.1.0/'
    },
    'HGNC': {
        'name': 'HUGO Gene Nomenclature Committee',
        'version': '2020-08-06',
        'namespace_prefix': 'HGNC',
        'id': 'HGNC:2020-08-06',
        'iri_prefix': 'http://rest.genenames.org/fetch/symbol/',
        'url': 'https://www.genenames.org/tools/search/'
    },
    'HGVS': {
        'name': 'Human Genome Variation Society',
        'version': '2020-05-01',
        'namespace_prefix': 'HGVS',
        'id': 'HGVS:2020-05-01',
        'iri_prefix': 'https://www.ncbi.nlm.nih.gov/clinvar/',
        'url': 'https://www.ncbi.nlm.nih.gov/clinvar/'
    }
}

ontologies = {
    'sex': {
        'female': 'FEMALE',
        'male': 'MALE',
        'not specified': 'UNKNOWN_SEX',
        'other': 'OTHER_SEX',
        'None': 'UNKNOWN_SEX'
    },
    'race': {
        'american indian or alaska native': '1002-5',
        'asian': '2028-9',
        'black or african american': '2054-5',
        'native hawaiian or pacific islander': '2076-8',
        'white': '2106-3',
        'unknown': 'UNK',
        'asked but no answer': 'ASKU'
    },
    'ethnicity': {
        'hispanic': '2135-2',
        'latino': '2135-2',
        'non-hispanic': '2186-5',
        'non-latino': '2186-5',
        'unknown': 'UNK',
        'not reported': 'ASKU'
    },
    'ecog_performance_status': {
        '0': {
            'id': 'ECOG: 0',
            'label': 'Fully active, able to carry on all pre-disease performance without restriction'
        },
        '1': {
            'id': 'ECOG: 1',
            'label': 'Restricted in physically strenuous activity but ambulatory and able to carry out work of a light or sedentary nature, e.g., light house work, office work'
        },
        '2': {
            'id': 'ECOG: 2',
            'label': 'Ambulatory and capable of all selfcare but unable to carry out any work activities; up and about more than 50% of waking hours'
        },
        '3': {
            'id': 'ECOG: 3',
            'label': 'Capable of only limited selfcare; confined to bed or chair more than 50% of waking hours'
        },
        '4': {
            'id': 'ECOG: 4',
            'label': 'Completely disabled; cannot carry on any selfcare; totally confined to bed or chair'
        },
        '5': {
            'id': 'ECOG: 5',
            'label': 'Dead'
        }
    },
    'specimen_type': {
        'amniotic fluid': {
            'id': 'HL7:AMN',
            'label': 'Amniotic fluid'
        },
        'bile Fluid': {
            'id': 'HL7:BIFL',
            'label': 'Bile Fluid'
        },
        'whole blood': {
            'id': 'HL7:BLD',
            'label': 'Whole blood'
        },
        'blood arterial': {
            'id': 'HL7:BLDA',
            'label': 'Blood arterial'
        },
        'cord blood': {
            'id': 'HL7:BLDCO',
            'label': 'Cord blood'
        },
        'blood venous': {
            'id': 'HL7:BLDV',
            'label': 'Blood venous'
        },
        'bone': {
            'id': 'HL7:BON',
            'label': 'Bone'
        },
        'convalescent serum': {
            'id': 'HL7:CSERU',
            'label': 'Serum, Convalescent'
        },
        'cerebral spinal fluid': {
            'id': 'HL7:CSF',
            'label': 'Cerebral spinal fluid'
        },
        'cervical mucus': {
            'id': 'HL7:CVM',
            'label': 'Cervical Mucus'
        },
        'cuodenal fluid': {
            'id': 'HL7:DUFL',
            'label': 'Duodenal fluid'
        },
        'fetal blood': {
            'id': 'HL7:FBLOOD',
            'label': 'Blood, Fetal'
        },
        'abdomen fluid': {
            'id': 'HL7:FGA',
            'label': 'Fluid, Abdomen'
        },
        'genital vaginal': {
            'id': 'HL7:GENV',
            'label': 'Genital vaginal'
        },
        'hydrocele fluid': {
            'id': 'HL7:HYDC',
            'label': 'Fluid, Hydrocele'
        },
        'joint fluid': {
            'id': 'HL7:JNTFLD',
            'label': 'Fluid, Joint'
        },
        'kidney fluid': {
            'id': 'HL7:KIDFLD',
            'label': 'Fluid, Kidney'
        },
        'lumbar sac fluid': {
            'id': 'HL7:LSAC',
            'label': 'Fluid, Lumbar Sac'
        },
        'marrow': {
            'id': 'HL7:MAR',
            'label': 'Marrow'
        },
        'pancreatic fluid': {
            'id': 'HL7:PAFL',
            'label': 'Pancreatic fluid'
        },
        'pericardial fluid': {
            'id': 'HL7:PCFL',
            'label': 'Fluid, Pericardial'
        },
        'placenta': {
            'id': 'HL7:PLC',
            'label': 'Placenta'
        },
        'pleural fluid': {
            'id': 'HL7:PLR',
            'label': 'Pleural fluid (thoracentesis fluid)'
        },
        'saliva': {
            'id': 'HL7:SAL',
            'label': 'Saliva'
        },
        'skin': {
            'id': 'HL7:SKN',
            'label': 'Skin'
        },
        'seminal fluid': {
            'id': 'HL7:SMN',
            'label': 'Seminal fluid'
        },
        'synovial fluid': {
            'id': 'HL7:SNV',
            'label': 'Fluid, synovial (Joint fluid)'
        },
        'sputum': {
            'id': 'HL7:SPT',
            'label': 'Sputum'
        },
        'tissue': {
            'id': 'HL7:TISS',
            'label': 'Tissue'
        },
        'vitreous fluid': {
            'id': 'HL7:VITF',
            'label': 'Vitreous Fluid'
        },
        'wound': {
            'id': 'HL7:WND',
            'label': 'Wound'
        }
    },
    'variant_data_value': {
        'positive': {
            'id': 'LOINC:LA9633-4',
            'label': 'Present'
        },
        'negative': {
            'id': 'LOINC:LA9634-2',
            'label': 'Absent'
        },
        'no call': {
            'id': 'LOINC:LA18198-4',
            'label': 'No call'
        },
        'indeterminate': {
            'id': 'LOINC:LA11884-6',
            'label': 'Indeterminate'
        }
    },
    'cancer_disease_status': {
        'not detected':
        {
            'id': 'SNOMED:260415000',
            'label': 'Not detected (qualifier)'
        },
        'improved': {
            'id': 'SNOMED:268910001',
            'label': 'Patient condition improved (finding)'
        },
        'stable': {
            'id': 'SNOMED:359746009',
            'label': 'Patient\'s condition stable(finding)'
        },
        'worsened': {
            'id': 'SNOMED:271299001',
            'label': 'Patient\'s condition worsened(finding)'
        },
        'undetermined': {
            'id': 'SNOMED:709137006',
            'label': 'Patient condition undetermined (finding)'
        }
    }
}
