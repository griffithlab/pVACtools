.. image:: ../images/pVACseq_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACseq logo

.. raw:: html

   <style> .large {font-size: 110%; font-weight: bold} </style>
   <style> .large-code {font-size: 110%; font-family: monospace} </style>


Frequently Asked Questions
==========================

.. role:: large
.. role:: large-code

:large:`What type of variants does pVACseq support?`

pVACseq makes predictions for all transcripts of a variant that were annotated
as ``missense_variant``, ``inframe_insertion``, ``inframe_deletion``,
inframe ``protein_altering_variant``, or
``frameshift_variant`` by VEP as long as the transcript was not also annotated
as ``start_lost``. In addition, pVACseq only includes variants that were
called as homozygous or heterozygous variant. Variants that were not called
in the sample specified are skipped (determined by examining the ``GT`` genotype
field in the VCF). In addition, some variants might be skipped in cases where
the VEP annotation does not contain protein position information.

:large:`My pVACseq command has been running for a long time. Why is
that?`

The rate-limiting factor in running pVACseq is the number of calls that are
made to the IEDB software for binding score predictions.

.. note::

  It is generally faster to make IEDB calls using a local install of IEDB than
  using the IEDB web API. It is, therefore, recommended to use a local IEDB
  install for any in-depth analysis. You should either install IEDB locally yourself
  or use the pvactools docker image that includes it.

There are a number of factors that determine the number of IEDB calls to be made:

- Number of variants in your VCF

  pVACseq will make predictions for each missense, inframe insertion,
  inframe deletion, protein altering, and frameshift variant in your VCF.

  **Speedup suggestion**: Split the VCF into smaller subsets and process each one
  individually, in parallel.

- Number of transcripts for each variant

  pVACseq will make predictions for each transcript of a supported variant
  individually. The number of transcripts for each variant depends on how VEP was
  run when the VCF was annotated.

  **Speedup suggestion**: Use the ``--pick`` option when running VEP to
  annotate each variant with the top transcript only.

- The ``--fasta-size`` parameter value

  pVACseq takes an input VCF and creates a wildtype and a mutant
  FASTA for each transcript. The number of FASTA entries that get submitted
  to IEDB at a time is limited by the ``--fasta-size`` parameter in order
  to reduce the load on the IEDB servers. The smaller the FASTA size, the
  more calls have to be made to IEDB.

  **Speedup suggestion**: When using a local IEDB install, increase the size
  of this parameter.

- Number of prediction algorithms, epitope lengths, and HLA-alleles

  One call to IEDB is made for each combination of these parameters for each chunk
  of FASTA sequences. That means, for example, when 8 prediction
  algorithms, 4 epitope lengths (8-11), and 6 HLA-alleles are chosen, 7*4*6=192 calls
  to IEDB have to be made for each chunk of FASTA.

  **Speedup suggestion**: Reduce the number of prediction algorithms,
  epitope lengths, and/or HLA-alleles to the ones that will be the most
  meaningful for your analysis. For example, the NetMHCcons method is
  already a consensus method between NetMHC, NetMHCpan, and PickPocket.
  If NetMHCcons is chosen, you may want to omit the underlying prediction
  methods. Likewise, if you want to run NetMHC, NetMHCpan, and PickPocket
  individually, you may want to skip NetMHCcons.

- ``--downstream-sequence-length`` parameter value

  This parameter determines how many amino acids of the downstream sequence
  after a frameshift mutation will be included in the wildtype FASTA sequence.
  The shorter the downstream sequence length, the lower the number of epitopes
  that IEDB needs to make binding predictions for.

  **Speedup suggestion**: Reduce the value of this parameter.

- ``-t`` parameter value

  This parameter determines the number of threads pvacseq will use for parallel 
  processing.

  **Speedup suggestion**: Use a host with multiple cores and sufficient memory and 
  use a larger number of threads.

:large:`My pVACseq all_epitopes output file does not contain entries for all of the
variants in the VCF. Why is that?`

There are a variety of reasons why a variant in your VCF might've been skipped
by pVACseq.

pVACseq only makes predictions for those transcripts of a variant that were annotated
as ``missense_variant``, ``inframe_insertion``, ``inframe_deletion``,
inframe ``protein_altering_variant``, or
``frameshift_variant`` by VEP as long as the transcript was not also annotated
as ``start_lost``.

In addition, pVACseq only includes variants that were
called as homozygous or heterozygous variant. Variants that were not called
in the sample specified are skipped (determined by examining the ``GT`` genotype
field in the VCF).

Some variants might be skipped in cases where
the VEP annotation does not contain protein position information.

If the ``--pass-only`` flag is set, variants with a filter in the VCF FILTER
column will be skipped.

Variants may also be skipped if the resulting peptide sequence includes an
unsupported amino acid or other character. Only ``A``, ``R``, ``N``, ``D``, ``C``, ``E``,
``Q``, ``G``, ``H``, ``I``, ``L``, ``K``, ``M``, ``F``, ``P``, ``S``, ``T``, ``W``, ``Y``,
and ``V`` are supported.

Lastly, a variant may be skipped if the resulting peptide sequence would be
too short for the chosen epitope length(s).

:large:`My pVACseq all_epitopes output file does not contain entries for all of the
alleles I chose. Why is that?`

There could be a few reasons why the pVACseq output does not contain
predictions for alleles:

- The alleles you picked might have not been compatible with the prediction algorithm and/or epitope lengths chosen. In that case no calls for that allele would've been made and a status message would've printed to the screen.

- It could be that all epitope predictions for some alleles got filtered out. You can check the ``<sample_name>.all_epitopes.tsv`` file to see all called epitopes before filtering.

:large:`Why are some values in the` :large-code:`WT Epitope Seq` :large:`column` :large-code:`NA` :large:`?`

Not all mutant epitope sequences will have a corresponding wildtype epitope sequence. 
This occurs when the mutant epitope sequence is novel and a comparison is therefore not
meaningful. For example:

- An epitope in the downstream portion of a frameshift might not have a corresponding wildtype epitope at the same position at all. The epitope is completely novel.

- An epitope that overlaps an inframe indel or multinucleotide polymorphism (MNP) might have a large number of amino acids that are different from the wildtype epitope at the corresponding position. If less than half of the amino acids between the mutant epitope sequence and the corresponding wildtype sequence match, the corresponding wildtype sequence in the report is set to ``NA``.

:large:`What filters are applied during a pVACseq run?`

By default we filter the neoepitopes on their binding score. If readcount
and/or expression annotations are available in the VCF we also filter on the depth, VAF,
and gene/trancript FPKM. In addition, candidates where the mutant epitope sequence is the
same as the wildtype epitope sequence will also be filtered out (i.e., they
don't overlap the mutation). pVACseq also filters on the transcript support
level, if the ``--tsl`` option was chosen during VEP annotation. Lastly, the
top score filter will pick the best epitope for each variant.

:large:`How can I see all of the candidate epitopes without any filters
applied?`

The ``<sample_name>.all_epitopes.tsv`` will contain all of the epitopes predicted
before filters are applied.

:large:`Why have some of my epitopes been filtered out even though the` :large-code:`Best MT Score` :large:`is below 500?`

By default, the binding filter will be applied to the ``Median MT Score``
column. This is the median score value among all chosen prediction algorithms.
The ``Best MT Score`` column shows the lowest score among all
chosen prediction algorithms. To change this behavior and apply the binding
filter to the ``Best MT Score`` column you may set the ``--top-score-metric``
parameter to ``lowest``.

:large:`Why are entries with` :large-code:`NA` :large:`in the`
:large-code:`VAF` :large:`and` :large-code:`depth` :large:`columns not
filtered?`

We do not filter out ``NA`` entries for depth and VAF since there is not
enough information to determine whether the cutoff has been met one way or another.

:large:`Why do some of my epitopes have no score predictions for certain prediction methods?`

Not all prediction methods support all epitope lengths or all alleles. To see
a list of supported alleles for a prediction method you may use the
``pvacseq valid_alleles`` :ref:`command <valid_alleles>`. For more details on
each algorithm refer to the IEDB MHC `Class I <http://tools.iedb.org/mhci/help/#Method>`_
and `Class II <http://tools.iedb.org/mhcii/help/#Method>`_ documentation.

:large:`How is pVACseq licensed?`

pVACseq is licensed under the open source license `NBSD 3-Clause Clear License
<https://spdx.org/licenses/BSD-3-Clause-Clear.html>`_.

:large:`How do I cite pVACseq?`

Jasreet Hundal+, Susanna Kiwala+, Joshua McMichael, Christopher A Miller,
Alexander T Wollam, Huiming Xia, Connor J Liu, Sidi Zhao, Yang-Yang Feng,
Aaron P Graubert, Amber Z Wollam, Jonas Neichin, Megan Neveau, Jason Walker,
William E Gillanders, Elaine R Mardis, Obi L Griffith, Malachi Griffith.
`pVACtools: a computational toolkit to select and visualize cancer
neoantigens <https://doi.org/10.1101/501817>`_. (+)equal contribution.
bioRxiv 501817; doi: https://doi.org/10.1101/501817

Jasreet Hundal, Susanna Kiwala, Yang-Yang Feng, Connor J. Liu, Ramaswamy Govindan, 
William C. Chapman, Ravindra Uppaluri, S. Joshua Swamidass, Obi L. Griffith, 
Elaine R. Mardis, and Malachi Griffith. `Accounting for proximal variants improves 
neoantigen prediction <https://www.nature.com/articles/s41588-018-0283-9>`_. Nature Genetics. 
2018, DOI: 10.1038/s41588-018-0283-9. PMID: `30510237 <https://www.ncbi.nlm.nih.gov/pubmed/30510237>`_.

Jasreet Hundal, Beatriz M. Carreno, Allegra A. Petti, Gerald P. Linette, Obi
L. Griffith, Elaine R. Mardis, and Malachi Griffith. `pVACseq: A genome-guided
in silico approach to identifying tumor neoantigens <http://www.genomemedicine.com/content/8/1/11>`_. Genome Medicine. 2016,
8:11, DOI: 10.1186/s13073-016-0264-5. PMID: `26825632
<http://www.ncbi.nlm.nih.gov/pubmed/26825632>`_.
