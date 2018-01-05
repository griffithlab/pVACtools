.. image:: images/pVACvector_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACvector logo

pVACvector
====================================

pVACvector is designed to aid specifically in the construction of
DNA-based cancer vaccines. It takes as input either a pVACseq output tsv file
or a fasta file containing peptide sequences and returns an ordering that
minimizes the effects of junctional epitopes (that may create novel peptides)
between the sequences. It does this by using the core pVACseq services to
predict the binding scores for each junctional peptide. It also tests
junctions with spacer amino acid sequences that may help to reduce reactivity.
These spacer amino acid sequences can be "HH", "HHC", "HHH", "HHHD", "HHHC",
"AAY", "HHHH", "HHAA", "HHL" or  "AAL. The final vaccine ordering is
achieved through a simulated annealing procedure that returns a near-optimal
solution, when one exists.

.. toctree::
   :glob:

   pvacvector/prerequisites
   pvacvector/getting_started
   pvacvector/run
   pvacvector/additional_commands
   pvacvector/output_files
