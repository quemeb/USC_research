# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 21:56:49 2023

@author: quemeb
"""

import requests
import json

host = 'http://bioghost2.usc.edu:3403'

# Perform a count query first to get the total number of results
count_url = f"{host}/annoq-test/_count"

count_data = {
    "query": {
        "bool": {
            "filter": [
                {"term": {"chr": "2"}},
                {"range": {"pos": {"gte": 10, "lte": 20000}}}
            ]
        }
    }
}

count_response = requests.post(count_url, headers={'Content-Type': 'application/json'}, json=count_data)
count_json = count_response.json()
total_count = count_json.get('count', 0)

# Now perform the actual search query with the size set to the total count
search_url = f"{host}/annoq-test/_search?pretty"

search_data = {
    "size": total_count,
    "query": {
        "bool": {
            "filter": [
                {"term": {"chr": "20"}},
                {"range": {"pos": {"gte": 26188698, "lte": 26190631}}}
            ]
        }
    }
}

search_response = requests.post(search_url, headers={'Content-Type': 'application/json'}, json=search_data)
search_json = search_response.json()

# Parse the JSON response
# Note: The following is a placeholder, you will need to adjust the code according to the actual JSON structure.
gene_ids = [hit.get("_source", {}).get("SnpEff_ensembl_Gene_ID") for hit in search_json.get("hits", {}).get("hits", [])]

# Print only the Gene IDs and the length of the list
print(len(gene_ids))
