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
        if not result:
            print("Unable to connect to URL: {} result {}".format(url, r.status_code))
        return result
    except ConnectionError as err:
        print("Unable to connect to URL: {}".format(url))
        return False


def format_url(url):
    if not url.startswith("http"):
        if url.startswith("//"):
            return "http:" + url
        else:
            return "http://" + url
    return url


def search_for_articles(domain_parts):
    search_terms = " or ".join(domain_parts)
    print("Running esearch.")
    handle = Entrez.esearch(db=DATABASE, retmax=NCBI_RETMAX, term=search_terms)
    record = Entrez.read(handle)
    handle.close()
    count = record['Count']
    print("Found", count, 'articles.')
    return count, record['IdList']


def parse_abstract_data(record):
    pmid = str(record['MedlineCitation']['PMID'])
    article = record['MedlineCitation']['Article']
    title = article['ArticleTitle'].encode('latin-1', 'replace')
    abstract = None
    if article.get('Abstract'):
        so = article['Abstract']['AbstractText'][0]
        abstract = so.encode('latin-1', 'replace')
    return pmid, title, abstract


def find_urls_for_domain_parts(outfilename, domain_parts, email):
    Entrez.email = email
    num_found_articles, article_id_list = search_for_articles(domain_parts)
    cnt = 0
    ids = ",".join(article_id_list)
    num_urls = 0
    with open(outfilename, 'w') as outfile:
        handle = Entrez.efetch(db=DATABASE, id=ids, retmode="xml")
        records = Entrez.parse(handle)
        for record in records:
            pmid, title, abstract = parse_abstract_data(record)
            if abstract:
                try:
                    tokens = nltk.word_tokenize(abstract)
                    urls_in_abstract = [format_url(token) for token in tokens if is_source_code_url(domain_parts, token)]
                    for url in urls_in_abstract:
                        num_urls += 1
                        outfile.write("{}\t{}\t{}\n".format(pmid, title, url))
                except UnicodeDecodeError as err:
                    print("Failed to parse abstract for {}: {}".format(pmid, err))
            cnt += 1
            if cnt % PROG_CNT == 0:
                sys.stdout.write('\rProcessed {} of {}'.format(cnt, num_found_articles))
        handle.close()
    print("Urls in Articles=", num_urls)

if len(sys.argv) != 3:
    print("Usage: python find_code_urls.py <OUTPUT_FILENAME> <YOUR_EMAIL>")
else:
    find_urls_for_domain_parts(sys.argv[1], DOMAIN_PARTS, sys.argv[2])
