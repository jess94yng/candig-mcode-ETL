#!/usr/bin/env python
# coding: utf-8

# Functions for ontology search through APIs

import json
import urllib.request
import urllib.error
import urllib.parse
from pprint import pprint
from Bio import Entrez
import httplib2 as http

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse


REST_URL = 'http://data.bioontology.org'

#function for fetching BioPortal ontologies
def get_json(url, key):
    API_KEY = key
    opener = urllib.request.build_opener()
    opener.addheaders = [('Authorization', 'apikey token=' + API_KEY)]
    return json.loads(opener.open(url).read())


#building BioPortal ontology search urls
def get_url(ontology, term, exact_match, subtree_url):
    if exact_match:
        if subtree_url != '':
            url = REST_URL + '/search?q=' + term + '&ontology='+ontology+'&required_exact_match=true' + \
                '&include=notation,prefLabel,synonym' + '&subtree_root_id='+subtree_url
        else:
            url = REST_URL + '/search?q=' + term + '&ontology='+ontology + \
                '&required_exact_match=true'+'&include=notation,prefLabel,synonym'
    else:
        if subtree_url != '':
            url = REST_URL + '/search?q=' + term + '&ontology='+ontology + \
                '&include=notation,prefLabel,synonym' + '&subtree_root_id='+subtree_url
        else:
            url = REST_URL + '/search?q=' + term + '&ontology=' + \
                ontology+'&include=notation,prefLabel,synonym'

    return url


#specific subtrees for ontologies
subtree_dict = {
    'collection_body': 'http%3A%2F%2Fpurl.bioontology.org%2Fontology%2FSNOMEDCT%2F442083009',
    'history_morphology_behaviour': {
        'subtreeOne': 'http%3A%2F%2Fpurl.bioontology.org%2Fontology%2FSNOMEDCT%2F367651003',
        'subtreeTwo': 'http%3A%2F%2Fpurl.bioontology.org%2Fontology%2FSNOMEDCT%2F399919001',
        'subtreeThree': 'http%3A%2F%2Fpurl.bioontology.org%2Fontology%2FSNOMEDCT%2F399983006'
    }
}


#fetch id and label for gene studied
def hgnc_api(symbol):
    headers = {
        'Accept': 'application/json',
    }

    uri = 'http://rest.genenames.org'
    path = '/fetch/symbol/' + symbol

    target = urlparse(uri+path)
    method = 'GET'
    body = ''

    h = http.Http()

    response, content = h.request(
        target.geturl(),
        method,
        body,
        headers)

    if response['status'] == '200':
        # assume that content is a json reply
        # parse content with the json module
        data = json.loads(content)
        return data


#fetch id for gene mutation
def hgvs_api_esearch(toSearch, email):
    Entrez.email = email
    handle = Entrez.esearch(db='clinvar', term=toSearch, retmax=15)
    record = Entrez.read(handle)
    handle.close()
    return(record)


#fetch label for gene mutation
def hgvs_api_esummary(identifier, email):
    Entrez.email = email
    handle = Entrez.esummary(db='clinvar', id=identifier, retmode='xml')
    record = Entrez.read(handle)
    handle.close()
    return(record['DocumentSummarySet']['DocumentSummary'][0]['title'])
