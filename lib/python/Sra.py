#!/usr/bin/env python
""" A set of tools for parsing SRA database XML dumps. """
from xml.etree import ElementTree as ET
import re

from Bio import Entrez

# List of headers from the Run Info output by SRA website Save as -> File -> RunInfo
SraRunInfo_Headers = ['Run', 'ReleaseDate', 'LoadDate', 'spots', 'bases',
        'spots_with_mates', 'avgLength', 'size_MB', 'AssemblyName',
        'download_path', 'Experiment', 'LibraryName', 'LibraryStrategy',
        'LibrarySelection', 'LibrarySource', 'LibraryLayout', 'InsertSize',
        'InsertDev', 'Platform', 'Model', 'SRAStudy', 'BioProject',
        'Study_Pubmed_id', 'ProjectID', 'Sample', 'BioSample', 'SampleType',
        'TaxID', 'ScientificName', 'SampleName', 'g1k_pop_code', 'source',
        'g1k_analysis_group', 'Subject_ID', 'Sex', 'Disease', 'Tumor',
        'Affection_Status', 'Analyte_Type', 'Histological_Type', 'Body_Site',
        'CenterName', 'Submission', 'dbgap_study_accession', 'Consent',
        'RunHash', 'ReadHash']

# List of headers from the Run table output by SRA website Send to Run Table -> Download
SraRunTable_Headers = ['Assay_Type_s', 'AssemblyName_s', 'BioProject_s',
        'BioSample_s', 'Center_Name_s', 'Experiment_s', 'InsertSize_l',
        'LibraryLayout_s', 'LibrarySelection_s', 'LibrarySource_s',
        'Library_Name_s', 'LoadDate_s', 'MBases_l', 'MBytes_l', 'Platform_s',
        'ReleaseDate_s', 'Run_s', 'SRA_Sample_s', 'SRA_Study_s',
        'Sample_Name_s', 'Sex_s', 'genotype_s', 'source_name_s', 'strain_s',
        'tissue_s', 'Consent_s', 'Organism_s', 'g1k_analysis_group_s',
        'g1k_pop_code_s', 'source_s']

# Headers that Overlap in the two tables
InBoth_Headers = ['AssemblyName', 'BioProject', 'BioSample', 'Consent',
        'Experiment', 'LibraryLayout', 'LibrarySelection',
        'LibrarySource', 'LoadDate', 'Platform', 'ReleaseDate', 'Run',
        'Sex', 'g1k_analysis_group', 'g1k_pop_code', 'source']

class SraResultsTable(object):
    def __init__(self, fname):
        tree = ET.parse(fname)
        root = tree.getroot()
        self.experiment_package = root.getchildren()

#         self.rowList = []
# 
#         for package in experiment_package:
#             for run in package.findall('RUN_SET/RUN'):
#                 row = {}
#                 row['LoadDate'] = ''
#                 row['bases'] = run.get('total_bases')
#                 row['spots_with_mates'] = self.get_spots_with_mates(run)
#                 row['avgLength'] = self.get_avgLength(run)
#                 row['size_MB'] = ''
#                 row['AssemblyName'] = ''
#                 row['download_path'] = ''
#                 row['Experiment'] = package.find('EXPERIMENT/IDENTIFIERS/PRIMARY_ID').text
#                 row['LibraryName'] = package.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME').text
#                 row['LibraryStrategy'] = package.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY').text
#                 row['LibrarySelection'] = package.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SELECTION').text
#                 row['LibrarySource'] = package.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE').text
#                 row['LibraryLayout'] = self.get_LibraryLayout(package)
#                 row['InsertSize'] = self.get_InsertSize(package)
#                 row['InsertDev'] = ''
#                 row['Platform'] = self.get_Platform(package)
#                 row['Model'] = self.get_Model(package, row['Platform'])
#                 row['SRAStudy'] = package.find('STUDY/IDENTIFIERS/PRIMARY_ID').text
#                 row['BioProject'] = self.get_BioProject(package)
#                 row['Study_Pubmed_id'] = ''
#                 row['ProjectID'] = ''
#                 row['Sample'] = package.find('SAMPLE/IDENTIFIERS/PRIMARY_ID').text
#                 row['BioSample'] = self.get_BioSample(package)
#                 row['SampleType'] = ''
#                 row['TaxID'] = package.find('SAMPLE/SAMPLE_NAME/TAXON_ID').text
#                 row['ScientificName'] = package.find('SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME').text

    def wrap_try(func):
        def wrapper(self, *args):
            try:
                return func(self, *args)
            except:
                return 'UNDEFINED'
        return wrapper

    @wrap_try
    def get_Run(self, run):
        return run.find('IDENTIFIERS/PRIMARY_ID').text

    @wrap_try
    def get_ReleaseDate(self, run):
        return run.get('published').split(' ')[0]

    @wrap_try
    def get_LoadDate(self, package):
        return ''

    @wrap_try
    def get_spots(self, run):
        spots = run.get('total_spots') 
        if spots is not None:
            return spots
        else:
            return 0

    @wrap_try
    def get_bases(self, run):
        bases = run.get('total_bases') 
        if bases is not None:
            return bases
        else:
            return 0

    def get_spots_with_mates(self, run):
        try:
            reads = run.findall('Statistics/Read')
            mates = min(reads[0].get('count'), reads[1].get('count')) 
            return mates
        except:
            return 0

    def get_avgLength(self, run):
        try:
            if int(run.find('Statistics').get('nreads')) == 1:
                read = run.find('Statistics/Read')
                return int(read.get('average'))
            else:
                reads = run.findall('Statistics/Read')
                length = int(reads[0].get('average')) + int(reads[1].get('average'))
                return length
        except:
            return 0

    @wrap_try
    def get_size_MB(self, package):
        return ''

    @wrap_try
    def get_AssemblyName(self, package):
        return ''

    @wrap_try
    def get_download_path(self, run):
        return ''

    @wrap_try
    def get_Experiment(self, package):
        return package.find('EXPERIMENT/IDENTIFIERS/PRIMARY_ID').text

    @wrap_try
    def get_LibraryName(self, package):
        return package.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME').text

    @wrap_try
    def get_LibraryStrategy(self, package):
        return package.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY').text

    @wrap_try
    def get_LibrarySelection(self, package):
        return package.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SELECTION').text

    @wrap_try
    def get_LibrarySource(self, package):
        return package.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE').text

    @wrap_try
    def get_LibraryLayout(self, package):
        try:
            package.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PARIED').tag
            return 'PAIRED'
        except:
            package.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/SINGLE').tag
            return 'SINGLE'

    def get_InsertSize(self, package):
        try:
            length = package.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PARIED').get('NOMINAL_LENGTH')
            if length is not None:
                return length
        except:
            pass

        try:
            length = package.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/SINGLE').get('NOMINAL_LENGTH')
            if length is not None:
                return length
        except:
            pass

        return 0

    @wrap_try
    def get_InsertDev(self, package):
        return 0

    @wrap_try
    def get_Platform(self, package):
        platform = package.find('EXPERIMENT/PLATFORM')
        child = platform.getchildren()
        if len(child) == 1:
            return child[0].tag
        else:
            return 'TOO_MANY_CHILDREN_CHECK_XML'

    @wrap_try
    def get_Model(self, package):
        platform = package.find('EXPERIMENT/PLATFORM')
        child = platform.getchildren()
        if len(child) == 1:
            return child[0].find('INSTRUMENT_MODEL').text
        else:
            return 'TOO_MANY_CHILDREN_CHECK_XML'

    @wrap_try
    def get_SRAStudy(self, package):
        return package.find('STUDY/IDENTIFIERS/PRIMARY_ID').text

    @wrap_try
    def get_BioProject(self, package):
        external = package.find('STUDY/IDENTIFIERS/EXTERNAL_ID')
        if external.get('namespace') == 'BioProject':
            return external.text

    @wrap_try
    def get_Study_Pubmed_id(self, package):
        return ''

    @wrap_try
    def get_ProjectID(self, package):
        return ''

    @wrap_try
    def get_Sample(self, package):
        return package.find('SAMPLE/IDENTIFIERS/PRIMARY_ID').text

    @wrap_try
    def get_BioSample(self, package):
        external = package.find('SAMPLE/IDENTIFIERS/EXTERNAL_ID')
        if external.get('namespace') == 'BioSample':
            return external.text
        else:
            raise ValueError

    @wrap_try
    def get_SampleType(self, package):
        return ''

    @wrap_try
    def get_TaxID(self, package):
        return package.find('SAMPLE/SAMPLE_NAME/TAXON_ID').text

    @wrap_try
    def get_ScientificName(self, package):
        return package.find('SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME').text

    @wrap_try
    def get_SampleName(self, package):
        return package.find('SAMPLE/SAMPLE_NAME/COMMON_NAME').text

    @wrap_try
    def get_g1k_pop_code(self, package):
        pass

    @wrap_try
    def get_source(self, package):
        pass

    @wrap_try
    def get_g1k_analysis_group(self, package):
        pass

    @wrap_try
    def get_Subject_ID(self, package):
        pass

    @wrap_try
    def get_Sex(self, package):
        pass

    @wrap_try
    def get_Disease(self, package):
        pass

    @wrap_try
    def get_Tumor(self, package):
        pass

    @wrap_try
    def get_Affection_Status(self, package):
        pass

    @wrap_try
    def get_Analyte_Type(self, package):
        pass

    @wrap_try
    def get_Histological_Type(self, package):
        pass

    @wrap_try
    def get_Body_Site(self, package):
        pass

    @wrap_try
    def get_CenterName(self, package):
        pass

    @wrap_try
    def get_Submission(self, package):
        pass

    @wrap_try
    def get_dbgap_study_accession(self, package):
        pass

    @wrap_try
    def get_Consent(self, package):
        pass

    @wrap_try
    def get_RunHash(self, package):
        pass

    @wrap_try
    def get_ReadHash(self, package):
        pass

    def render(self):
        pass


def downloadSRA(count, webenv, query_key, batch_size=500, fname='../../output/sra_dump.xml'):
    """ Function to download SRA results in batches.
    
    SRA prefers that large queries are downloaded in batches. Here 
    I use a batch size of 500. Note because I am concatenating XML I do
    a little string maniuplation so that it is a valid XML.

    Parameters
    ----------

    count: str or int
        Number of quieries to download
    webenv: str
        Session cookie from eSearch
    query_key: str
        Key to access the query from eSearch

    batch_size: int
        Number of quieries to download at one time

    fname: str
        XML file name
    

    Example
    -------

    >>> handle = Entrez.esearch(db='sra', term='"Drosophila melanogaster"[Orgn]', 
            retmax=100, usehistory='y')
    >>> records = Entrez.read(handle)
    >>> webenv = records['WebEnv']
    >>> query_key = records['QueryKey']
    >>> Sra.downloadSRA(count=records['Count'], webenv=webenv, 
                        query_key=query_key, batch_size=10)

    """
    from urllib.error import HTTPError
    import time
    out_handle = open(fname, "w")
    out_handle.write('<?xml version="1.0" ?>\n<EXPERIMENT_PACKAGE_SET>\n')
    count = int(count)
    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)
        print("Going to download record %i to %i" % (start+1, end))
        attempt = 1
        while attempt <= 3:
            try:
                fetch_handle = Entrez.efetch(db="sra", retstart=start, retmax=batch_size,
                                             webenv=webenv, query_key=query_key)
                break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt)
                    attempt += 1
                    time.sleep(15)
                else:
                    raise
        data = fetch_handle.read()
        fetch_handle.close()
        data = re.sub('<\?xml version="1.0" \?>\n<EXPERIMENT_PACKAGE_SET>\n', '', data)
        data = re.sub('</EXPERIMENT_PACKAGE_SET>\n\n', '', data)
        out_handle.write(data)
    out_handle.write('</EXPERIMENT_PACKAGE_SET>')
    out_handle.close()
