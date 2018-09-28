.. image:: ../images/pVACviz_logo_trans-bg_sm_v4b.png
    :align: right
    :alt: pVACviz logo

pVACapi Directories
===================

The pVACapi process that supports pVACviz features creates several directories in the user's home folder in ~/pVACseq/. These directories are used to hold input files to pVACseq processes, results files for visualization, archives and exported projects.

archive
---------
pVACviz provides an archive function within its :ref:`Manage section <pvacviz_managing>`. When processes are archived, they are placed in this archive folder.

export
--------
pVACviz provides an export function within its :ref:`Manage section <pvacviz_managing>`. When processes are exported, they are placed in this export folder.

input
-------
The pVACviz :ref:`Start form <pvacviz_starting_processes>` has Input VCF and Phased Proximal Variant fields that accept VCF files. The selectors for these fields list all relevant files placed within the `~/input` directory. You may sort these files into directories of any depth and the selectors will keep them grouped by directory.

visualize
---------
The :ref:`Visualize feature <pvacviz_visualizing>` allows users to visualize any pVACseq result VCF files. Any pVACseq VCF file placed in this /visualize folder will be displayed on the Visualize page in the right column. Directory structures will be preserved so that users may group files in whatever way they wish. 
