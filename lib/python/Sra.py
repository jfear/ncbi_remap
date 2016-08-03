#!/usr/bin/env python
""" A set of tools for parsing SRA database XML dumps. """
from xml.etree import ElementTree as ET
import re

import pandas as pd

from Bio import Entrez

# List of headers from the Run Info output by SRA website Save as -> File -> RunInfo
SraRunInfo_Headers = ['Run', 'RunSecondary', 'ReleaseDate', 'LoadDate', 'spots', 'bases',
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

CV = {'sex': ['sex']}

class SraResultsTable(object):
    def __init__(self, fname):
        tree = ET.parse(fname)
        root = tree.getroot()
        self.experiment_package = root.getchildren()

    def build_rows(self):
        experiments = []
        for package in self.experiment_package:
            row = {}
            row['LoadDate'] = self.get_LoadDate(package)
            row['size_MB'] = self.get_size_MB(package)
            row['AssemblyName'] = self.get_AssemblyName(package)
            row['download_path'] = self.get_download_path(package)
            row['Experiment'] = self.get_Experiment(package)
            row['LibraryName'] = self.get_LibraryName(package)
            row['LibraryStrategy'] = self.get_LibraryStrategy(package)
            row['LibrarySelection'] = self.get_LibrarySelection(package)
            row['LibrarySource'] = self.get_LibrarySource(package)
            row['LibraryLayout'] = self.get_LibraryLayout(package)
            row['InsertSize'] = self.get_InsertSize(package)
            row['InsertDev'] = self.get_InsertDev(package)
            row['Platform'] = self.get_Platform(package)
            row['Model'] = self.get_Model(package)
            row['SRAStudy'] = self.get_SRAStudy(package)
            row['BioProject'] = self.get_BioProject(package)
            row['Study_Pubmed_id'] = self.get_Study_Pubmed_id(package)
            row['ProjectID'] = self.get_ProjectID(package)
            row['Sample'] = self.get_Sample(package)
            row['BioSample'] = self.get_BioSample(package)
            row['SampleType'] = self.get_SampleType(package)
            row['TaxID'] = self.get_TaxID(package)
            row['ScientificName'] = self.get_ScientificName(package)
            row['SampleName'] = self.get_SampleName(package)
            row['g1k_pop_code'] = self.get_g1k_pop_code(package)
            row['source'] = self.get_source(package)
            row['g1k_analysis_group'] = self.get_g1k_analysis_group(package)
            row['Subject_ID'] = self.get_Subject_ID(package)
            row['Sex'] = self.get_Sex(package)
            row['Disease'] = self.get_Disease(package)
            row['Tumor'] = self.get_Tumor(package)
            row['Affection_Status'] = self.get_Affection_Status(package)
            row['Analyte_Type'] = self.get_Analyte_Type(package)
            row['Histological_Type'] = self.get_Histological_Type(package)
            row['Body_Site'] = self.get_Body_Site(package)
            row['Submission'] = self.get_Submission(package)
            row['dbgap_study_accession'] = self.get_dbgap_study_accession(package)

            rows = []
            try:
                for run in package.findall('RUN_SET/RUN'):
                    current_run = {}
                    current_run['Run'] = self.get_Run(run)
                    current_run['RunSecondary'] = self.get_RunSecondary(run)
                    current_run['ReleaseDate'] = self.get_ReleaseDate(run)
                    current_run['spots'] = self.get_spots(run)
                    current_run['bases'] = self.get_bases(run)
                    current_run['spots_with_mates'] = self.get_spots_with_mates(run)
                    current_run['avgLength'] = self.get_avgLength(run)
                    current_run['CenterName'] = self.get_CenterName(run)
                    current_run['Consent'] = self.get_Consent(run)
                    current_run['RunHash'] = self.get_RunHash(run)
                    current_run['ReadHash'] = self.get_ReadHash(run)
                    current_run.update(row)
                    rows.append(pd.Series(current_run))
                experiments.extend(rows)
            except:
                pass

        return pd.concat(experiments, axis=1).T[SraRunInfo_Headers]


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

    def get_RunSecondary(self, run):
        try:
            secondary = run.findall('IDENTIFIERS/SECONDARY_ID')
            return ';'.join([x.text for x in secondary])
        except:
            return ''

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
        # NOTE: It looks like ERR samples use COMMON_NAME as their SampleName.
        # I think it is still better to use the alias because it is a more
        # descriptive name
        return package.find('SAMPLE').get('alias')

    @wrap_try
    def get_g1k_pop_code(self, package):
        return ''

    @wrap_try
    def get_source(self, package):
        return ''

    @wrap_try
    def get_g1k_analysis_group(self, package):
        return ''

    @wrap_try
    def get_Subject_ID(self, package):
        return ''

    def get_Sex(self, package):
        try:
            attribs = package.findall('SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')
            for attrib in attribs:
                if attrib.tag.lower() in CV['sex']:
                    return attrib.value
        except:
            pass

        return ''

    @wrap_try
    def get_Disease(self, package):
        return ''

    @wrap_try
    def get_Tumor(self, package):
        return 'no'

    @wrap_try
    def get_Affection_Status(self, package):
        return ''

    @wrap_try
    def get_Analyte_Type(self, package):
        return ''

    @wrap_try
    def get_Histological_Type(self, package):
        return ''

    @wrap_try
    def get_Body_Site(self, package):
        try:
            attribs = package.findall('SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')
            for attrib in attribs:
                if attrib.tag.lower() == 'body_site':
                    return attrib.value
        except:
            pass

        return ''

    def get_CenterName(self, run):
        try:
            return run.get('run_center')
        except:
            pass
        return ''

    @wrap_try
    def get_Submission(self, package):
        return package.find('SUBMISSION/IDENTIFIERS/PRIMARY_ID').text

    @wrap_try
    def get_dbgap_study_accession(self, package):
        return ''

    @wrap_try
    def get_Consent(self, run):
        if run.get('is_public') == 'true':
            return 'public'
        else:
            return 'not public'

    @wrap_try
    def get_RunHash(self, package):
        return ''

    @wrap_try
    def get_ReadHash(self, package):
        return ''

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
