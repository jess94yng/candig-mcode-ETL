#!/usr/bin/env python
# coding: utf-8

# In[2]:
# Formulate identifications and labels for fields using codeable concepts and various ontololgies

from IntegrateAPI import get_json, get_url, subtreeDict, hgncApi, hgvsApiESearch, hgvsApiESummary


# In[3]:


def generateGeneticSpecimen(termStr,apiKey):
    strList=termStr.split(' ')
    needBreak=False
    for term in strList:
        termCollection=get_json(get_url("SNOMEDCT",term,True,subtreeDict['collection_body']),apiKey)['collection']
        if len(termCollection)>0:
            for listItem in termCollection:
                if 'notation' in listItem:
                    synonymList=listItem['synonym']
                    for synonym in synonymList:
                        if term.lower()==synonym.lower():
                            identification='SNOMED:'+ listItem['notation']
                            label=listItem['prefLabel']
                            needBreak=True
                            break
                if needBreak:
                    break
        if needBreak:
            break
        if term==strList[-1]:
            identification= "SNOMED: 261665006"
            label='Unknown'
    
    idLabel=[identification,label]
    return idLabel


# In[4]:


def generateBodySite(termStr,apiKey):
    if '/' in termStr:
        bodyPart = termStr.split('/')[0]
    else:
        bodyPart=termStr
    bodyPartNew=bodyPart.replace(' ','+')
    termCollection = get_json(get_url("SNOMEDCT",bodyPartNew,False,subtreeDict['collection_body']),apiKey)['collection']
    idLabel = generateIdLabel(termCollection, 'SNOMEDCT', "SNOMED: 261665006")
    return idLabel


# In[ ]:


def generateMedication(medItem,apiKey):
    medString = medItem.replace(' ','+')
    termCollection = get_json(get_url("RxNorm",medString,False,''),apiKey)['collection']
    idLabel = generateIdLabel(termCollection, 'RxNorm', 'RxNorm:Unknown')
    return idLabel
        


# In[ ]:


def generateTumorMarker(tumorItem, tumorList,apiKey):
    if tumorItem.replace('.','1').isdigit():
        tumorList[-1]['tumor_marker_data_value'] = {
            'value': {
                'value': float(tumorItem),
                'comparator': '=',
            }
        }
    elif '/' in tumorItem:
        num=tumorItem.split('/')
        tumorList[-1]['tumor_marker_data_value'] = {
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
        dataStr = tumorItem.replace(' ','+')
        termCollection = get_json(get_url("loinc",dataStr,False,''),apiKey)['collection']
        idLabel = generateIdLabel(termCollection, "LOINC", 'LOINC:LA4489-6')
        tumorList[-1]['tumor_marker_data_value'] = {
            'value': {
            'id': idLabel[0],
            'label': idLabel[1]
            }
        }
    return tumorList


# In[1]:


def generateIdLabel(collection, dbName, unknownId):
    if len(collection)>0:
        for listItem in collection:
            if 'notation' in listItem:
                identification = dbName + ':' + listItem['notation']
                label = listItem['prefLabel']
                break
            if listItem==collection[-1]:
                identification = unknownId
                label = 'Unknown'
    else:
        identification = unknownId
        label = 'Unknown'

    idLabel=[identification,label]
    return idLabel

