# Code_in_PubMed_Abstracts

## Current functionality

This tool searches for code cited in PubMed abstracts by searching for popular repositories such as GitHub and Bitbucket (full list below).  

It returns PMIDs and URLs where code is presumed to be.  

## Future functionality

It will check for valid URLs.  

## Current list of repositories

1. github
2. bitbucket
3. dataverse
4. figshare
5. bioconda
6. omictools
7. sourceforge
8. bioinformatics.org
9. bioinformatics.ca
10. iubio.bio.indiana.edu
11. bioweb.pasteur.fr
12. bio.tools
13. scicrunch.org
14. identifiers.org

* please feel free to suggest additional repos in issues!

## Notes

In its current manifestation, it will pull a maximum of 100,000 abstracts


## Setup
The tool should work in both python2 and python3. Using a [virtualenv](https://www.dabapps.com/blog/introduction-to-pip-and-virtualenv-python/) is strongly recommended.

```
$ pip install -r requirements.txt
```

You may need to install the Python development headers (Python.h). 
Consult your distribution for details. On Ubuntu, this is likely:

```
$ apt-get install -y libpython2.7-dev # or libpython3.5-dev 
```

On CentOS: 

```
$ yum install -y python-devel
```

Download nltk required files:

```
$ python -c 'import nltk ; nltk.download()'
```

Download `punkt` data files. By default this will end up in `$HOME/nltk_data`.

## Usage
```
python find_code_urls.py <OUTPUT_CSV_FILENAME> <YOUR_EMAIL>
```
Writes a tsv file containing a record per url we were able to find in pubmed articles.
Columns: pmid, article title, url

