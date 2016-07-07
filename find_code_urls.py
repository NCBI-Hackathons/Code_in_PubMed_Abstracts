"""
Finds urls in pubmed articles and writes them out to a tsv file.
Checks for articles containing DOMAIN_PARTS below and verifies that it is a valid url.
"""

from __future__ import print_function
from Bio import Entrez
import nltk
import requests
from requests.exceptions import ConnectionError
import sys
import string
import re

# parts of the url that we should search for
DOMAIN_PARTS = [
    "github",
    "bitbucket",
    "dataverse",
    "figshare",
    "bioconda",
    "omictools",
    "sourceforge",
    "bioinformatics.org",
    "bioinformatics.ca",
    "iubio.bio.indiana.edu",
    "bioweb.pasteur.fr",
    "vbio.tools",
    "scicrunch",
    "identifiers.org",
]
DATABASE = "pubmed"

# maximum records that can be returned
# expected to find much less than this...
NCBI_RETMAX = 10000

# How often to print out progress
PROG_CNT = 10


def is_source_code_url(domain_parts, token):
    """
    Does the token contain one of our domain parts and resolves to a url.
    :param domain_parts: [str] domain part of url
    :param token: str possible code url token
    :return: Boolean: true if token is a source code url
    """
    for part in domain_parts:
        if part in token:
            if token != part:
                return is_url_valid(token)
    return False


def is_url_valid(url):
    """
    Is the specified URL valid? Fetches it to see if it works.
    :param url: str: url to check
    """
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
    """
    Add http:// or similar if necessary.
    :param url: url to format
    :return: formatted url
    """
    url = url.rstrip("/")
    if not url.startswith("http"):
        if url.startswith("//"):
            return "http:" + url
        else:
            return "http://" + url
    return url


def search_for_articles(domain_parts):
    """
    Search for articles containing domain_parts in pubmed.
    :param domain_parts: [str] list of terms to search for
    :return: int, list[ID]
    """
    search_terms = " or ".join(domain_parts)
    print("Running esearch.")
    handle = Entrez.esearch(db=DATABASE, retmax=NCBI_RETMAX, term=search_terms)
    record = Entrez.read(handle)
    handle.close()
    count = record['Count']
    print("Found", count, 'articles.')
    return count, record['IdList']


def clean_str(value):
    """
    Remove any non printable characters from the value.
    :param value: str value to be cleaned
    :return: str cleaned up value
    """
    return ''.join([c for c in value if c in string.printable])


def parse_article_data(record):
    """
    Pull out fields from NCBI article
    :param record:
    :return: pmid, title, abstract from the record
    """
    pmid = str(record['MedlineCitation']['PMID'])
    article = record['MedlineCitation']['Article']
    title = clean_str(article['ArticleTitle'])
    abstract = None
    if article.get('Abstract'):
        abstract = clean_str(article['Abstract']['AbstractText'][0])
    return pmid, title, abstract


def find_urls_in_abstract(domain_parts, abstract, pmid):
    """
    Find urls in the specified abstract by looking for domain parts.
    :param domain_parts: [str] list of terms to search for
    :param abstract: str text to search
    :param pmid: str id of the abstract
    :return: [str] urls we found
    """
    try:
        tokens = nltk.word_tokenize(abstract)
        urls_in_abstract = [format_url(token) for token in tokens if is_source_code_url(domain_parts, token)]
        return urls_in_abstract
    except UnicodeDecodeError as err:
        print("Failed to parse abstract for {}: {}".format(pmid, err))
        return []


def efetch_articles(article_id_list):
    """
    Fetch articles for a list of pmids
    :param article_id_list: [str] article pmids to fetch
    :return: handle to results
    """
    ids = ",".join(article_id_list)
    return Entrez.efetch(db=DATABASE, id=ids, retmode="xml")


def show_progress_message(cnt, num_found_articles):
    if cnt % PROG_CNT == 0:
        sys.stdout.write('\rProcessed {} of {}'.format(cnt, num_found_articles))


def find_urls_for_domain_parts(outfilename, domain_parts):
    """
    Write pmid, title, url records to outfilename for each article that contains domain_parts.
    :param outfilename: str: filename to write tsv into
    :param domain_parts: [str] list of terms to search for
    """
    num_found_articles, article_id_list = search_for_articles(domain_parts)
    cnt = 0
    found_urls = 0
    with open(outfilename, 'w') as outfile:
        handle = efetch_articles(article_id_list)
        records = Entrez.parse(handle)
        for record in records:
            pmid, title, abstract = parse_article_data(record)
            if abstract:
                urls_in_abstract = find_urls_in_abstract(domain_parts, abstract, pmid)
                for url in urls_in_abstract:
                    outfile.write("{}\t{}\t{}\n".format(pmid, title, url))
                    found_urls += 1
            cnt += 1
            show_progress_message(cnt, num_found_articles)
        handle.close()

    print("\nFound {} urls.".format(found_urls))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python find_code_urls.py <OUTPUT_FILENAME> <YOUR_EMAIL>")
    else:
        Entrez.email = sys.argv[2]
        find_urls_for_domain_parts(sys.argv[1], DOMAIN_PARTS)
