import sys
import os
import unittest
import py_compile
import vcf

from pvactools.lib.proximal_variant import ProximalVariant
from .test_utils import *

class ProximalVariantTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.executable_dir = os.path.join(pvactools_directory(), 'pvactools', 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'proximal_variant.py')
        cls.test_data_dir  = os.path.join(pvactools_directory(), 'tests', 'test_data', 'proximal_variant')
        proximal_variant_vcf_path = os.path.join(cls.test_data_dir, 'input.vcf.gz')
        cls.klass = ProximalVariant(proximal_variant_vcf_path, False, 30)
        somatic_vcf_path = os.path.join(cls.test_data_dir, 'somatic.vcf.gz')
        cls.somatic_vcf_fh = open(somatic_vcf_path, 'rb')
        cls.somatic_vcf_reader = vcf.Reader(cls.somatic_vcf_fh)

    @classmethod
    def tearDownClass(cls):
        cls.klass.fh.close()
        cls.somatic_vcf_fh.close()

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_no_potential_proximal_variants(self):
        somatic_variant = next(self.somatic_vcf_reader.fetch('chr1', 16006133, 16006134))
        (phased_somatic_variant, potential_proximal_variants) = self.klass.find_phased_somatic_variant_and_potential_proximal_variants(somatic_variant, "T", "ENST00000329454")

        self.assertTrue(phased_somatic_variant)
        self.assertFalse(potential_proximal_variants)

    def test_found_potential_proximal_variants(self):
        somatic_variant = next(self.somatic_vcf_reader.fetch('chr2', 227893862, 227893863))
        (phased_somatic_variant, potential_proximal_variants) = self.klass.find_phased_somatic_variant_and_potential_proximal_variants(somatic_variant, "T", "ENST00000309931")

        self.assertTrue(phased_somatic_variant)
        self.assertTrue(potential_proximal_variants)

    def test_found_actual_proximal_variant(self):
        somatic_variant = next(self.somatic_vcf_reader.fetch('chr2', 227893862, 227893863))
        proximal_variants = self.klass.extract(somatic_variant, "T", "ENST00000309931")
        self.assertTrue(proximal_variants)

    def test_combine_conflicting_variants(self):
        self.assertEqual(self.klass.combine_conflicting_variants(["ttC/ttA", "tTc/tAc"]), '*')
        self.assertEqual(self.klass.combine_conflicting_variants(["Cca/Tca", "cCa/cTa"]), 'L')

