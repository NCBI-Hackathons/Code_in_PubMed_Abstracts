from __future__ import print_function
from Bio import Entrez
import nltk
import time

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


def is_source_code_url(domain_parts, token):
    for part in domain_parts:
        if part in token:
            return token != part
    return False


def print_urls_for_domain_parts(domain_parts, email):
    Entrez.email = EMAIL
    search_terms = " or ".join(domain_parts)
    handle = Entrez.esearch(db=DATABASE, retmax=NCBI_RETMAX, term=search_terms)
    record = Entrez.read(handle)
    cnt = 0
    for id in record['IdList']:
        handle2 = Entrez.efetch(db=DATABASE, id=id, rettype="abstract", retmode="text")
        text = handle2.read().decode('utf-8').strip()
        tokens = nltk.word_tokenize(text)
        gh_tokens = [token for token in tokens if is_source_code_url(domain_parts, token)]
        if gh_tokens:
            print("id", id, "urls", gh_tokens)
        handle2.close()
        cnt += 1
        if cnt == SLEEP_AFTER_CNT:
            time.sleep(SLEEP_AMT)
            cnt = 0
    handle.close()


print_urls_for_domain_parts(DOMAIN_PARTS, EMAIL)
