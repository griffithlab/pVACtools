.. image:: images/pVACvector_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACvector logo

pVACvector
====================================

pVACvector is designed to aid specifically in the construction of
DNA vector based personalized cancer vaccines. It takes as input either a pVACseq output 
tsv file or a FASTA file containing peptide sequences and returns a peptide ordering that
minimizes the effects of junctional epitopes (that may create novel peptides)
between the sequences. It does this by using the core pVACseq services to
predict the binding scores for each junctional peptide.

Running pVACvector with spacer amino acid sequences may help eliminate junctional
epitopes. The list of spacers to be tested is specified using the ``--spacers``
parameter. Peptide combinations without a spacer can be tested by including
``None`` in the list of spacers. The default spacer amino acid sequences are
"None", "AAY", "HHHH", "GGS", "GPGPG", "HHAA", "AAL", "HH", "HHC", "HHH", "HHHD",
"HHL", "HHHC". Peptide junctions are tested with each spacer in the order that
they are specified. If a tested spacers results in a valid junction without any
well-binding junction epitopes, that junction will not be tested with any
other spacers, even if a different spacer could potentially result in better
junction scores. This reduces runtime. If a tested spacer for a junction doesn't
yield a valid junction (i.e., there are well-binding junction epitopes) the junction
is tested wtih the next spacer in the input list.

If, after testing all spacers, no valid path is found, clipped versions of
peptides are tested by removing leading and/or trailing amino acids and
constructing junctions with the clipped peptides. The maximum number of amino
acids to clip is controlled by the ``--max-clip-length`` argument.

In some cases, the (core) neoantigen candidate of a peptide sequence may be located
toward the beginning or end of the sequence. In these cases, clipping may
accidentially remove amino acids of the core neoantigen. To prevent this, the
``--max-clip-length`` should be set to the shortest number of flanking amino
acids of any of the peptides to include in the vector. Alternatively, pVACvector also
supports specifying the core neoantigen in the FASTA header when using a FASTA
file as the input to pVACvector. If the core neoantigens for each sequence are specified in the
input FASTA file, pVACvector will not clip into these neoantigens, even if the
flanking sequence is smaller than the ``--max-clip-length``. The core neoantigen should
be specified like so:

.. code-block:: none

    >Peptide1 {"Best Peptide": "LYYSYGLLHI"}
    WLYYSYGLLHIYGSGGYALYF

In this example ``Peptide1`` is the ID of the sequence, ``LYYSYGLLHI`` is
the core neoantigen candidate, and ``WLYYSYGLLHIYGSGGYALYF`` is the peptide
sequence to include in the vector. The Best Peptide information will already
be included in the FASTA headers if the FASTA file is created by using the ``pvacseq
generate_protein_fasta`` command in conjunction with an aggregated report TSV
as the ``--input-tsv`` parameter.

If no solution is found after testing all spacers and after clipping peptides, pVACvector
will attempt to find a partial solution by excluding peptide sequences. The
number of peptide sequences that are allowed to be removed is controlled via
the ``--allow-n-peptide-exclusion`` parameter. Partial solutions will be
written to their own result subdirectory. The subdirectory name reflects which
peptide(s) were removed from the partial solution.

The final vaccine ordering is achieved through a simulated annealing procedure that
returns a near-optimal solution, when one exists.

.. toctree::
   :glob:

   pvacvector/prerequisites
   pvacvector/getting_started
   pvacvector/run
   pvacvector/additional_commands
   pvacvector/output_files
