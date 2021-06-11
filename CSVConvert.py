#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import json
from datetime import date
from MappingDict import raceMappingDict, validMCode, mCodeDict, biochemistry, labResults, vitalSigns, possibleMultiVal
from Ontologies import Ontologies
from SchemaBase import returnMCodePacket, returnSubject, returnGenomicsReport, returnGeneticRegionStudied, returnCancerCondition, returnCancerRelatedProcedures, returnTumorMarker

from IntegrateAPI import get_json, get_url, subtreeDict, hgncApi, hgvsApiESearch, hgvsApiESummary
from SchemaGeneration import generateGeneticSpecimen, generateBodySite, generateIdLabel, generateMedication, generateTumorMarker


def notNull(cell):
    if cell == 'nan': return False
    else: return True

def raceMapping(inputDataframe,ind,curRow,mappingList):
    for listItem in mCodeDict['race']:
        if listItem in inputDataframe:
            if curRow[listItem] != '0':
                if curRow[listItem] == '1':
                    mappingList.append(('race',raceMappingDict[listItem]))
                elif pd.notna(curRow[listItem]):
                    mappingList.append(('race',curRow[listItem]))

def eliminateDuplicate(inputDataframe):
    inputDataframe = inputDataframe.dropna(how='all')
    inputDataframe = inputDataframe.applymap(str)
    for index, row in inputDataframe.iterrows():
        for multiValItem in possibleMultiVal:
            if multiValItem in inputDataframe:
                row[multiValItem] = row[multiValItem].split(',')

        if(inputDataframe.duplicated(subset=['identifier'],keep = False)[index]):
            for multiValItem in possibleMultiVal:
                if multiValItem in inputDataframe:
                    newRow=[]
                    for index1, row1 in inputDataframe.iterrows():
                        if index1!=index and row1['identifier'] == row['identifier'] and row1[multiValItem] != 'nan':
                            row[multiValItem].append(row1[multiValItem])
    inputDataframe = inputDataframe.drop_duplicates(subset = 'identifier')
    return inputDataframe

# Process csv data, use mappings to return data frame with mCode data elements
def mCodeMapping(inputDataframe):

    mCodeCSVList = []

    for index, row in inputDataframe.iterrows():
        dataElementList = []
        dictItem = {}
        for item in validMCode:
            if item == 'race':
                raceMapping(inputDataframe, index, row, dataElementList)
            elif mCodeDict[item] in inputDataframe:
                dataElementList.append((item,row[mCodeDict[item]]))

        for element in dataElementList:
            dictItem[element[0]] = element[1]

        mCodeCSVList.append(dictItem)
    
    df1 = pd.DataFrame(mCodeCSVList)
    df1 = eliminateDuplicate(df1)
    return df1


# Combine dataframes from multiple sheets, delete any duplicate patients by merging data
def combineSheets(listDataframe):
    for ind in range (1, len(listDataframe)):
        listDataframe[0] = listDataframe[0].merge(listDataframe[ind], how='outer')
    outputDf = listDataframe[0]
    for index, row in outputDf.iterrows():
        if(outputDf.duplicated(subset=['identifier'],keep = False)[index]):
            for column in outputDf:
                notnaRow = False
                if isinstance (row[column],list) or pd.notna(row[column]): notnaRow = True
                for index1, row1 in outputDf.iterrows():
                    notnaRow1 = False
                    if isinstance (row1[column],list) or pd.notna(row1[column]): notnaRow1 = True
                    if index1>index and row1['identifier'] == row['identifier'] and notnaRow == False and notnaRow1 == True:
                        row[column] = row1[column]
    outputDf = outputDf.drop_duplicates(subset = 'identifier')
    return outputDf


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
    if dataframeItem.get('issued',None) and notNull(dataframeItem['issued']):
        genomicsReport['issued']=dataframeItem['issued']
    else: genomicsReport['issued']='Unknown'
    if dataframeItem.get('genetic_specimen',None) and notNull(dataframeItem['genetic_specimen']):
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
        for i in range(len(dataframeItem['gene_mutation'])):
            result = hgvsApiESearch(dataframeItem['gene_mutation'][i])
            if len(result['IdList']) > 0:
                geneMutation.append({
                    'id' : "HGVS:" + result['IdList'][0],
                    'label': hgvsApiESummary(result['IdList'][0])
                })
        geneticRegionStudied['gene_mutation'] = geneMutation

    if dataframeItem.get('gene_studied',None) and notNull(dataframeItem['gene_studied']):
        geneRegionAdd = True
        geneStudied = []
        for i in range(len(dataframeItem['gene_studied'])):
            result = hgncApi(dataframeItem['gene_studied'][i])
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
        bodySite = []
        for i in range(len(dataframeItem['body_site'])):
            bodySite.append({
                'id': generateBodySite(dataframeItem['body_site'][i])[0],
                'label': generateBodySite(dataframeItem['body_site'][i])[1]
            })
            
        cancerRelatedProcedures['body_site'] = bodySite
    
    mCodePacket['cancer_related_procedures'] = cancerRelatedProcedures
    
    medicationStatement = []
    if dataframeItem.get('medication',None) and notNull(dataframeItem['medication']):
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
        tempCounter3 = 0
        for data in dataframeItem['tumor_marker_data_value']:
            tumorMarker.append(returnTumorMarker(counterVar, tempCounter3, dataframeItem['identifier'], True))
            tumorMarker = generateTumorMarker(data, tumorMarker)

            if extra:
                tumorMarker[-1]['extra_properties']=tumorExtra
            tempCounter3 += 1
            
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
    


def main():
    genomicsID = 1000
    dfOg = pd.read_excel('data/test_data.xlsx', sheet_name=None, dtype=str)
    dfList=[]
    for page in dfOg:
        dfList.append(mCodeMapping(dfOg[page]))
    
    mCodeDataframe = combineSheets(dfList)
    mCodeDataframe = mCodeDataframe.fillna("nan")
    mCodeList=[]
    for index, row in mCodeDataframe.iterrows():
        mCodeList.append(mCodePacketBuild(row,genomicsID))
        genomicsID+=1
    print (mCodeList)
    with open('mCodePacket.json', 'w') as f:
        json.dump(mCodeList, f, indent=4)
        

if __name__ == '__main__':
    main()

