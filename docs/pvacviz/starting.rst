.. image:: ../images/pVACviz_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACviz logo

.. _pvacviz_starting:

Starting Processes
==================

pVACviz provides a helpful form for specifying all of the parameters for a pVACseq process, as an alternative to constructing these commands and executing them via the command line.


Populating and Submitting the Start Form
----------------------------------------
The form is divided into two sections: required and optional parameters. To be submitted, all required fields must be filled and validated. Optional parameters are pre-filled with sensible defaults - the same defaults that would be applied when submitting processes via the command line.

The form provides feedback as to which fields remain to be filled and validated, both with a red highlight and message around fields in question, and at the bottom of the form a list of incomplete or invalid fields.

Once all the required form fields are completed with valid values, the Submit Process button is activated. Clicking on this button will submit the process to pVACapi. Clicking the Reset button will restore the form to its initial pristine state.

Notes
-----
* The Input VCF and Phased Proximal Variant VCF fields require the selection of VCF files. These selectors list all VCF files within the :ref:`/input folder <pvacviz_directories>` found within the ~/pVACseq folder located in the user's home folder.

* The alleles selector only shows alleles relevant to selected prediction algorithms. Choose prediction algorithms to enable and populate the alleles selector.
