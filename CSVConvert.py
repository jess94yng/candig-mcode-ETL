#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import json
from datetime import date
from MappingDict import raceMappingDict, validMCode, mCodeDict, biochemistry, labResults, vitalSigns, possibleMultiVal
from Ontologies import Ontologies
from SchemaBase import returnMCodePacket, returnSubject, returnGenomicsReport, returnGeneticRegionStudied, returnCancerCondition, returnCancerRelatedProcedures, returnTumorMarker

from IntegrateAPI import get_json, get_url, subtreeDict, hgncApi, hgvsApiESearch, hgvsApiESummary
from SchemaGeneration import generateGeneticSpecimen, generateBodySite, generateIdLabel, generateMedication, generateTumorMarker


# In[2]:


def notNull(cell):
    if cell == 'nan': return False
    else: return True


# In[3]:

# Process csv data, use mappings to return data frame with mCode data elements
def mCodeMapping(inputDataframe):

    mCodeCSVList = []

    for index, row in inputDataframe.iterrows():
        dataElementList = []
        dictItem = {}
        for item in validMCode:
            if isinstance(mCodeDict[item],list):
                if item == 'race':
                    for listItem in mCodeDict[item]:
                        if listItem in inputDataframe:
                            listItemIndex = mCodeDict[item].index(listItem)
                            if row[mCodeDict[item][listItemIndex]] != '0':
                                if row[mCodeDict[item][listItemIndex]] == '1':
                                    dataElementList.append((item,raceMappingDict[listItem]))
                                elif pd.notna(row[mCodeDict[item][listItemIndex]]):
                                    dataElementList.append((item,row[mCodeDict[item][listItemIndex]]))
                                
                else:
                    listItems = []
                    for listItem in mCodeDict[item]:
                        if listItem in inputDataframe:
                            listItemIndex = mCodeDict[item].index(listItem)
                            listItems.append(row[mCodeDict[item][listItemIndex]])
                    if listItems:
                        dataElementList.append((item,listItems))

            elif mCodeDict[item] in inputDataframe:
                dataElementList.append((item,row[mCodeDict[item]]))

        for element in dataElementList:
            dictItem[element[0]] = element[1]

        mCodeCSVList.append(dictItem)
    
    df1 = pd.DataFrame(mCodeCSVList)
    df1 = df1.dropna(how='all')
    df1 = df1.applymap(str)
    for index, row in df1.iterrows():
        toAdd = False
        for multiValItem in possibleMultiVal:
            if multiValItem in df1:
                if ',' in row[multiValItem]:
                    row[multiValItem] = row[multiValItem].split(',')
                    toAdd = True
        if(df1.duplicated(subset=['identifier'],keep = False)[index]):
            for multiValItem in possibleMultiVal:
                if multiValItem in df1:
                    newRow=[row[multiValItem]]
                    for index1, row1 in df1.iterrows():
                        if index1!=index and row1['identifier'] == row['identifier'] and row1[multiValItem] != 'nan':
                            newRow.append(row1[multiValItem])
                    if toAdd:
                        row[multiValItem].append(newRow)
                    else:
                        row[multiValItem] = newRow
    df1 = df1.drop_duplicates(subset = 'identifier')
    return df1


# In[6]:

# Combine dataframes from multiple sheets, delete any duplicate patients by merging data
def combineSheets(listDataframe):
    for ind in range (1, len(listDataframe)):
        listDataframe[0] = listDataframe[0].merge(listDataframe[ind], how='outer')
    outputDf = listDataframe[0]
    outputDf = outputDf.applymap(str)
    for index, row in outputDf.iterrows():
        if(outputDf.duplicated(subset=['identifier'],keep = False)[index]):
            for column in outputDf:    
                for index1, row1 in outputDf.iterrows():
                    if index1>index and row1['identifier'] == row['identifier'] and row[column] == 'nan' and row1[column] != 'nan':
                        row[column] = row1[column]
    outputDf = outputDf.drop_duplicates(subset = 'identifier')
    return outputDf


# In[7]:

# Organize dataframe into mCode model
def mCodePacketBuild(dataframeItem,counterVar):
    mCodePacket = returnMCodePacket(dataframeItem['identifier'])
    
    subject = returnSubject(dataframeItem['identifier'])
    
    if dataframeItem.get('date_of_birth', None) and notNull(dataframeItem['date_of_birth']):
        subject['date_of_birth'] = dataframeItem['date_of_birth']
    if dataframeItem.get('birth_sex', None) and notNull(dataframeItem['birth_sex']):
        subject['sex'] = Ontologies['sex'][dataframeItem['birth_sex'].lower()]
    if dataframeItem.get('race',None) and notNull(dataframeItem['race']):
        subject['race'] = Ontologies['race'][dataframeItem['race'].lower()]
    if dataframeItem.get('ethnicity',None) and notNull(dataframeItem['ethnicity']):
        subject['ethnicity'] = Ontologies['ethnicity'][dataframeItem['ethnicity'].lower()]
    subject['comorbid_condition'] = {
        'code':{
            'id': 'SNOMED:103329007',
            'label': 'Not available'
        }
    }
    if dataframeItem.get('ecog_performance_status',None) and notNull(dataframeItem['ecog_performance_status']):
        subject['ecog_performance_status'] = Ontologies['ecog_performance_status'][dataframeItem['ecog_performance_status']]
    subject['karnofsky_performance_status'] = {
        'id': 'KARNOFSKY: Not available',
        'label': 'Not available'
    }
    
    for vitalItem in vitalSigns:
        if dataframeItem.get(vitalItem, None) and notNull(dataframeItem[vitalItem]):
            subject['extra_properties'][vitalItem]=dataframeItem[vitalItem]

    mCodePacket['subject'] = subject

    genomicsReport = returnGenomicsReport(counterVar)
    
    geneticSpecimen = []
    if dataframeItem.get('performing_organization_name',None) and notNull(dataframeItem['performing_organization_name']):
        genomicsReport['performing_organization_name'] = dataframeItem['performing_organization_name']
    genomicsReport['issued'] = dataframeItem.get('issued',"Unknown")
    if dataframeItem.get('genetic_specimen',None) and notNull(dataframeItem['genetic_specimen']):
        if isinstance(dataframeItem['genetic_specimen'],list):
            tempCounter1 = 0
            for i in range(len(dataframeItem['genetic_specimen'])):
                geneticSpecimen.append({
                    'id': str(counterVar) + '-' + str(tempCounter1),
                    'collection_body': {
                        'id': generateGeneticSpecimen(dataframeItem['collection_body_site'][i])[0],
                        'label': generateGeneticSpecimen(dataframeItem['collection_body_site'][i])[1]
                    }

                })
                tempCounter1 += 1
                
                if dataframeItem.get('specimen_type',None) and dataframeItem['specimen_type'] in Ontologies['specimen_type']:
                    geneticSpecimen[-1]['specimen_type'] = Ontologies['specimen_type'][dataframeItem['specimen_type'][i].lower()]
                else:
                    geneticSpecimen[-1]['specimen_type'] = {
                            'id': 'HL7:...',
                            'label': 'No suggested value'
                    }
        else:
            geneticSpecimen.append({
                'id': str(counterVar),
                'collection_body': {
                    'id': generateGeneticSpecimen(dataframeItem['collection_body_site'])[0],
                    'label': generateGeneticSpecimen(dataframeItem['collection_body_site'])[1]
                }
            })
            if dataframeItem.get('specimen_type',None) and dataframeItem['specimen_type'] in Ontologies['specimen_type']:
                geneticSpecimen[-1]['specimen_type'] = Ontologies['specimen_type'][dataframeItem['specimen_type'].lower()]
            else:
                geneticSpecimen[-1]['specimen_type'] = {
                        'id': 'HL7:...',
                        'label': 'No suggested value'
                }
                
        genomicsReport['genetic_specimen'] = geneticSpecimen
    
    if dataframeItem.get('genetic_variant', None) and notNull(dataframeItem['genetic_variant']):
        if dataframeItem['variant_data_value'] in Ontologies['variant_data_value']:
            genomicsReport['genetic_variant'] = {
                'id': str(counterVar),
                'data_value': Ontologies['variant_data_value'][dataframeItem['variant_data_value'].lower()]
            }
        else: 
            genomicsReport['genetic_variant'] = {
                'id':str(counterVar),
                'data_value': {
                    'id': 'LOINC:LA4489-6',
                    'label': 'Unknown'
                }
            }

    geneticRegionStudied = returnGeneticRegionStudied(counterVar)
    geneRegionAdd=False

    if dataframeItem.get('gene_mutation',None) and notNull(dataframeItem['gene_mutation']):
        geneRegionAdd = True
        geneMutation = []
        if isinstance(dataframeItem['gene_mutation'],list):
            for i in range(len(dataframeItem['gene_mutation'])):
                result = hgvsApiESearch(dataframeItem['gene_mutation'][i])
                if len(result['IdList']) > 0:
                    geneMutation.append({
                        'id' : "HGVS:" + result['IdList'][0],
                        'label': hgvsApiESummary(result['IdList'][0])
                    })
        else:
            result = hgvsApiESearch(dataframeItem['gene_mutation'])
            if len(result['IdList']) > 0:
                geneMutation.append({
                    'id' : "HGVS:" + result['IdList'][0],
                    'label': hgvsApiESummary(result['IdList'][0])
                })
        geneticRegionStudied['gene_mutation'] = geneMutation

    if dataframeItem.get('gene_studied',None) and notNull(dataframeItem['gene_studied']):
        geneRegionAdd = True
        geneStudied = []
        if isinstance(dataframeItem['gene_studied'],list):
            for i in range(len(dataframeItem['gene_studied'])):
                result = hgncApi(dataframeItem['gene_studied'][i])
                if result['response']['numFound'] >= 1:
                    geneStudied.append({
                        'id' : result['response']['docs'][0]['hgnc_id'],
                        'label': result['response']['docs'][0]['symbol']
                    })
        else:
            result = hgncApi(dataframeItem['gene_studied'])
            if result['response']['numFound'] >= 1:
                geneStudied.append({
                    'id' : result['response']['docs'][0]['hgnc_id'],
                    'label': result['response']['docs'][0]['symbol']
                })
        geneticRegionStudied['gene_studied'] = geneStudied
                    
    if geneRegionAdd:
        genomicsReport['genetic_region_studied'] = geneticRegionStudied
    
    mCodePacket['genomics_report'] = genomicsReport
    
    cancerCondition = returnCancerCondition(counterVar)
    
    if dataframeItem.get('date_of_diagnosis',None) and notNull(dataframeItem['date_of_diagnosis']):
        cancerCondition['date_of_diagnosis'] = dataframeItem['date_of_diagnosis']
    if dataframeItem.get('history_morphology_behavior',None) and notNull(dataframeItem['history_morphology_behavior']):
        behaviorTerm = dataframeItem['history_morphology_behavior'].replace(' ','+')
        needBreak = False
        for behaviorSubtree in subtreeDict['history_morphology_behaviour']:
            termCollection = get_json(get_url('SNOMEDCT',behaviorTerm,False,subtreeDict['history_morphology_behaviour'][behaviorSubtree]))['collection']
            if len(termCollection) > 0:
                for listItem in termCollection:
                    if 'notation' in listItem:
                        identification = 'SNOMED:' + listItem['notation']
                        label = listItem['prefLabel']
                        needBreak = True
                        break
            if needBreak:
                break
            if behaviorSubtree == 'subtreeThree':
                identification = 'SNOMED:261665006'
                label = 'Unknown'
        cancerCondition['history_morphological_behavior'] = {
            'id': identification,
            'label': label
        }
    
    
    mCodePacket['cancer_condition'] = cancerCondition
    
    cancerRelatedProcedures = returnCancerRelatedProcedures(counterVar)
    
    if dataframeItem.get('body_site',None) and notNull(dataframeItem['body_site']):
        if isinstance(dataframeItem['body_site'],list):
            bodySite = []
            for i in range(len(dataframeItem['body_site'])):
                bodySite.append({
                    'id': generateBodySite(dataframeItem['body_site'][i])[0],
                    'label': generateBodySite(dataframeItem['body_site'][i])[1]
                })
        else:
            bodySite = [{
                'id': generateBodySite(dataframeItem['body_site'])[0],
                'label': generateBodySite(dataframeItem['body_site'])[1]
            }]
            
        cancerRelatedProcedures['body_site'] = bodySite
    
    mCodePacket['cancer_related_procedures'] = cancerRelatedProcedures
    
    medicationStatement = []
    if dataframeItem.get('medication',None) and notNull(dataframeItem['medication']):
        if isinstance(dataframeItem['medication'],list):
            tempCounter2 = 0
            for medicine in dataframeItem['medication']:
                medicationStatement.append({
                    'id': str(counterVar) + '-' + str(tempCounter2),
                    'medication_code': {

                        'id': generateMedication(medicine)[0],
                        'label': generateMedication(medicine)[1]
                    }
                })
                tempCounter2 += 1
        else:
            medicationStatement.append({
                'id': str(counterVar),
                'medication_code':{

                    'id': generateMedication(dataframeItem['medication'])[0],
                    'label': generateMedication(dataframeItem['medication'])[1]
                }
            })
    else:
        medicationStatement.append({
            'id': str(counterVar),
            'medication_code': {
                'id': 'RxNorm:Unknown',
                'label': 'Unknown' 
            }
        })
    
    mCodePacket['medication_statement'] = medicationStatement
    
    if dataframeItem.get('date_of_death',None) and notNull(dataframeItem['date_of_death']):
        mCodePacket['date_of_death'] = dataframeItem['date_of_death']
    
    tumorMarker = []
    tumorExtra={}
    extra=False
    for extraItem in labResults:
        if dataframeItem.get(extraItem,None) and notNull(dataframeItem[extraItem]):
                tumorExtra[extraItem] = dataframeItem[extraItem]
                extra=True

    if dataframeItem.get('tumor_marker_data_value',None) and notNull(dataframeItem['tumor_marker_data_value']):
        if isinstance(dataframeItem['tumor_marker_data_value'],list):
            tempCounter3 = 0
            for data in dataframeItem['tumor_marker_data_value']:
                tumorMarker.append(returnTumorMarker(counterVar, tempCounter3, dataframeItem['identifier'], True))
                tumorMarker = generateTumorMarker(data, tumorMarker)
                
                if extra:
                    tumorMarker[-1]['extra_properties']=tumorExtra
            tempCounter3 += 1
        else:
            tumorMarker.append(returnTumorMarker(counterVar, 0, dataframeItem['identifier'], False))
            tumorMarker = generateTumorMarker(dataframeItem['tumor_marker_data_value'], tumorMarker)
            
            if extra:
                tumorMarker[-1]['extra_properties'] = tumorExtra
            
    else:
        tumorMarker.append(returnTumorMarker(counterVar, 0, dataframeItem['identifier'], False))
        
        if extra:
            tumorMarker[-1]['extra_properties'] = tumorExtra
    
    mCodePacket['tumor_marker'] = tumorMarker
    
    if dataframeItem.get('cancer_disease_status',None) and notNull(dataframeItem['cancer_disease_status']):
        if dataframeItem['cancer_disease_status'] in Ontologies['cancer_disease_status']:
            cancerDiseaseStatus = Ontologies['cancer_disease_status'][dataframeItem['cancer_disease_status'].lower()]
            mCodePacket['cancer_disease_status'] = cancerDiseaseStatus
    
    return (mCodePacket)
    


# In[8]:


def main():
    genomicsID = 1000
#     mCodeDataframe = mCodeMapping('data/Synthetic_data_v2_radiology.csv')
    dfOg = pd.read_excel('data/synthetic_data_v2.xlsx', sheet_name=None, dtype=str)
    dfList=[]
    for page in dfOg:
        dfList.append(mCodeMapping(dfOg[page]))
    
    mCodeDataframe = combineSheets(dfList)
    mCodeList=[]
    for index, row in mCodeDataframe.iterrows():
        mCodeList.append(mCodePacketBuild(row,genomicsID))
        genomicsID+=1
    with open('mCodePacket.json', 'w') as f:
        json.dump(mCodeList, f, indent=4)
        

if __name__ == '__main__':
    main()

