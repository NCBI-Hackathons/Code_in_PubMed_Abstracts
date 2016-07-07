from __future__ import print_function
from Bio import Entrez
import nltk
import time
import requests
from requests.exceptions import ConnectionError
import sys

DOMAIN_PARTS = [
    "github",
    "bitbucket",
    "dataverse",
    "figshare",
]
DATABASE = "pubmed"
EMAIL = "john.bradley@duke.edu"
# maximum records that can be returned
# expected to find much less than this...
NCBI_RETMAX = 100000
# How many abstracts to fetch before sleeping
SLEEP_AFTER_CNT = 100
# How long to sleep in seconds after fetching some ABSTRACTS to reduce load on NCBI servers.
SLEEP_AMT = 0.3
# How often to print out progress
PROG_CNT = 10

def is_source_code_url(domain_parts, token):
    for part in domain_parts:
        if part in token:
            if token != part:
                return is_url_valid(token)
    return False


def is_url_valid(url):
    try:
        r = requests.get(format_url(url))
        result = r.ok
        return result
    except ConnectionError as err:
        return False


def format_url(url):
    if not url.startswith("http"):
        if url.startswith("//"):
            return "http:" + url
        else:
            return "http://" + url
    return url


def print_urls_for_domain_parts(domain_parts, email):
    Entrez.email = EMAIL
    search_terms = " or ".join(domain_parts)
    print("Running esearch.")
    handle = Entrez.esearch(db=DATABASE, retmax=NCBI_RETMAX, term=search_terms)
    record = Entrez.read(handle)
    id_to_urls = {}
    cnt = 0
    print("Found", record['Count'], 'articles.')
    for id in record['IdList']:
        handle2 = Entrez.efetch(db=DATABASE, id=id, rettype="abstract", retmode="text")
        text = handle2.read().decode('utf-8').strip()
        tokens = nltk.word_tokenize(text)
        urls_in_abstract = [format_url(token) for token in tokens if is_source_code_url(domain_parts, token)]
        for url in urls_in_abstract:
            url_list = id_to_urls.get(id, None)
            if not url_list:
                url_list = []
                id_to_urls[id] = url_list
            url_list.append(url)
        handle2.close()
        cnt += 1
        if cnt % SLEEP_AFTER_CNT == 0:
            time.sleep(SLEEP_AMT)
        if cnt % PROG_CNT == 0:
            print("Fetched", cnt, "of", record['Count'])
    handle.close()
    return id_to_urls


outfile=sys.argv[1]
id_to_urls = print_urls_for_domain_parts(DOMAIN_PARTS, EMAIL)
print("Articles with code URLs", len(id_to_urls))
with open(outfile, 'w') as outfile:
    for id, urls in id_to_urls.items():
        for url in urls:
            outfile.write("{},{}\n".format(id, url))
