from unittest import TestCase
import Sra

fname = './test_xml.xml'


class TestSraResultsTableRowQuery(TestCase):
    def setUp(self):
        self.sra = Sra.SraResultsTable(fname)
        self.package = self.sra.experiment_package[0]
        self.run = self.package.find('RUN_SET/RUN')

    def test_get_Run(self):
        print('Run', self.sra.get_Run(self.run))

    def test_get_ReleaseDate(self):
        self.sra.get_ReleaseDate(self.run)

    def test_get_LoadDate(self):
        self.sra.get_LoadDate(self.package)

    def test_get_spots(self):
        self.sra.get_spots(self.run)

    def test_get_bases(self):
        self.sra.get_bases(self.run)

    def test_get_spots_with_mates(self):
        self.sra.get_spots_with_mates(self.run)

    def test_get_avgLength(self):
        self.sra.get_avgLength(self.run)

    def test_get_size_MB(self):
        self.sra.get_size_MB(self.package)

    def test_get_AssemblyName(self):
        self.sra.get_AssemblyName(self.package)

    def test_get_download_path(self):
        self.sra.get_download_path(self.package)

    def test_get_Experiment(self):
        self.sra.get_Experiment(self.package)

    def test_get_LibraryName(self):
        self.sra.get_LibraryName(self.package)

    def test_get_LibraryStrategy(self):
        self.sra.get_LibraryStrategy(self.package)

    def test_get_LibrarySelection(self):
        self.sra.get_LibrarySelection(self.package)

    def test_get_LibrarySource(self):
        self.sra.get_LibrarySource(self.package)

    def test_get_LibraryLayout(self):
        self.sra.get_LibraryLayout(self.package)

    def test_get_InsertSize(self):
        self.sra.get_InsertSize(self.package)

    def test_get_InsertDev(self):
        self.sra.get_InsertDev(self.package)

    def test_get_Platform(self):
        self.sra.get_Platform(self.package)

    def test_get_Model(self):
        self.sra.get_Model(self.package)

    def test_get_SRAStudy(self):
        self.sra.get_SRAStudy(self.package)

    def test_get_BioProject(self):
        self.sra.get_BioProject(self.package)

    def test_get_Study_Pubmed_id(self):
        self.sra.get_Study_Pubmed_id(self.package)

    def test_get_ProjectID(self):
        self.sra.get_ProjectID(self.package)

    def test_get_Sample(self):
        self.sra.get_Sample(self.package)

    def test_get_BioSample(self):
        self.sra.get_BioSample(self.package)

    def test_get_SampleType(self):
        self.sra.get_SampleType(self.package)

    def test_get_TaxID(self):
        self.sra.get_TaxID(self.package)

    def test_get_ScientificName(self):
        self.sra.get_ScientificName(self.package)

    def test_get_SampleName(self):
        print(self.sra.get_SampleName(self.package))

    def test_get_g1k_pop_code(self):
        pass

    def test_get_source(self):
        pass

    def test_get_g1k_analysis_group(self):
        pass

    def test_get_Subject_ID(self):
        pass

    def test_get_Sex(self):
        pass

    def test_get_Disease(self):
        pass

    def test_get_Tumor(self):
        pass

    def test_get_Affection_Status(self):
        pass

    def test_get_Analyte_Type(self):
        pass

    def test_get_Histological_Type(self):
        pass

    def test_get_Body_Site(self):
        pass

    def test_get_CenterName(self):
        pass

    def test_get_Submission(self):
        pass

    def test_get_dbgap_study_accession(self):
        pass

    def test_get_Consent(self):
        pass

    def test_get_RunHash(self):
        pass

    def test_get_ReadHash(self):
        pass

if __name__ == '__main__':
    import unittest
    unittest.main()







