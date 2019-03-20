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
predict the binding scores for each junctional peptide separated by a spacer amino acid
sequence that may help to eliminate junctional epitopes. The list of spacers
to be tested is specified using the ``--spacers`` parameter. Peptide combinations without a
spacer can be tested by including ``None`` in the list of spacers.

Peptide junctions are tested with
each spacer in the order that they are specified. If a valid peptide ordering
is found that doesn't result in any well-binding junction epitopes, that
ordering is returned. No other spacer are tested, even if they could
potentially result in better junction scores. This reduces runtime.
If no valid path is found, the next spacer in the input list is tested.
The default spacer amino acid sequences are "None", "AAY", "HHHH", "GGS", "GPGPG", "HHAA",
"AAL", "HH", "HHC", "HHH", "HHHD", "HHL", "HHHC".

The final vaccine ordering is
achieved through a simulated annealing procedure that returns a near-optimal
solution, when one exists.

.. toctree::
   :glob:

   pvacvector/prerequisites
   pvacvector/getting_started
   pvacvector/run
   pvacvector/additional_commands
   pvacvector/output_files
