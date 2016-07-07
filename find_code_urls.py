from __future__ import print_function
from Bio import Entrez
import nltk
import requests
from requests.exceptions import ConnectionError
import sys

DOMAIN_PARTS = [
    "github",
    "bitbucket",
    "dataverse",
    "figshare",
    "bioconda",
    "omictools.com",
    "sourceforge.net",
    "bioinformatics.org",
    "bioinformatics.ca",
    "iubio.bio.indiana.edu",
    "bioweb.pasteur.fr",
    "vbio.tools",
    "scicrunch.org",
    "identifiers.org",
]
DATABASE = "pubmed"
EMAIL = "john.bradley@duke.edu"
# maximum records that can be returned
# expected to find much less than this...
NCBI_RETMAX = 100000
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
        print("BAD URL: {}".format(url))
        return False


def format_url(url):
    if not url.startswith("http"):
        if url.startswith("//"):
            return "http:" + url
        else:
            return "http://" + url
    return url


def find_urls_for_domain_parts(outfilename, domain_parts, email):
    Entrez.email = EMAIL
    search_terms = " or ".join(domain_parts)
    print("Running esearch.")
    handle = Entrez.esearch(db=DATABASE, retmax=NCBI_RETMAX, term=search_terms)
    record = Entrez.read(handle)
    handle.close()
    id_to_urls = {}
    cnt = 0
    print("Found", record['Count'], 'articles.')
    ids = ",".join(record['IdList'])
    num_urls = 0
    with open(outfilename, 'w') as outfile:
        handle2 = Entrez.efetch(db=DATABASE, id=ids, retmode="xml")
        records = Entrez.parse(handle2)
        for record in records:
            id = str(record['MedlineCitation']['PMID'])
            article = record['MedlineCitation']['Article']
            title = article['ArticleTitle'].encode('latin-1', 'replace')
            if article.get('Abstract'):
                try:
                    so = article['Abstract']['AbstractText'][0]
                    abstract_str = so.encode('latin-1', 'replace')
                    text = abstract_str
                    tokens = nltk.word_tokenize(text)
                    urls_in_abstract = [format_url(token) for token in tokens if is_source_code_url(domain_parts, token)]
                    for url in urls_in_abstract:
                        num_urls += 1
                        outfile.write("{}\t{}\t{}\n".format(id, title, url))
                except UnicodeDecodeError as err:
                    print(err)
            cnt += 1
            if cnt % PROG_CNT == 0:
                print(cnt)
        handle2.close()
    print("Urls in Articles=", num_urls)


find_urls_for_domain_parts(sys.argv[1], DOMAIN_PARTS, EMAIL)
