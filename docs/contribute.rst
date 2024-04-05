How To Contribute
=================

pVACtools is an open-source project. We enthusiastically welcome outside contributions.

Making Code Contributions
_________________________

The `pVACtools source code <https://github.com/griffithlab/pVACtools>`_ is hosted on GitHub.
Code changes are made via pull requests to the main repository.
The `First Contributions
<https://github.com/firstcontributions/first-contributions>`_ project has a
neat guide of all the steps involved. Please follow the steps in this guide
when making code contributions to pVACtools, with the modifications as noted
below.

Development Branches
--------------------

The pVACtools project uses ``staging`` and
``hotfix`` branches to manage major/minor and patch releases, respectively.
Depending on the type of your contribution you will
need to branch your change branch off of either ``staging`` for new features
or ``hotfix`` for bugfixes and use one these two branches as the base branch
of your pull request. Pull requests should never be made to the ``master``
branch.

Testing
-------

To test code changes, users will first need to install pVACtools from their
local checkout of the pVACtools respository. From inside the main repository
folder run:

.. code-block:: none

   pip install -e .

Please ensure that all tests are passing by running:

.. code-block:: none

   TEST_FLAG=1 python -m unittest discover

Depending on what changes were made, test failures might be expected because,
for example, the data in output files has changed as part of the code changes.
In those cases, please update the affected test data files, as necessary, after
confirming that the updated outputs are correct.

When adding new features or fixing a bug, we ask that you add new test(s) to
cover these changes. Test files can be found under the ``tests`` directory and
generally correspond to the name of the pVACtools module files that are being
tested. Test data can be found under ``tests/test_data`` with subfolders
corresponding to each test file. The names of tests being added should
generally describe the feature of bug that is being tested.

If new test input files need to be added, we ask that they are as small in size
as possible. For VCF files, this can generally be achieved by removing
optional header lines and limiting the data to one or more entries that reflect
the issue being addressed. Please ensure that all PHI is removed from your
test data by, for example, using dummy sample names.

After making a GitHub pull request for your changes, GitHub will also run
these tests again on different versions of Python. All tests need to pass on
GitHub before we will consider merging a pull request. The coverage test may
be optional if the decrease in coverage percentage is small enough.

Milestones
----------

The pVACtools projects uses milestones to assign issues to specific releases.
Each milestone is named with a version number, which indicates the version of
the upcoming release. Upcoming milestones/versions can be found `here
<https://github.com/griffithlab/pVACtools/milestones>`_. An issue in an open
milestone may be closed if the respective feature has been addressed (i.e.
added via a merged pull request), even though the release for the open
milestone has not been made yet. Therefore, a closed issue does not indicate
that a feature is already available in a release yet. We urge users to
also search closed issues for any problem they may encounter and to check
whether a desired feature addition has already been made but isn't released
yet.


Updating documentation pages
____________________________

If you wish to make a suggested fix/improvement to the documentation for 
pVACtools, all documentation exists in the github pVACtools repo in the 
**docs** sub directory. The pVACtools docs are built using sphinx and published 
using readthedocs. You can either install sphinx locally on your system or 
you can use the pVACtools docker image. Sphinx allows you to build a local 
copy of the docs and see how the build looks in your browser before you 
commit/push and submit a pull request.  The following instructions detail 
how to run sphinx from the pVACtools docker container and view changes in 
the brower on your host system.

The following instructions assume: (1) that you already have Docker installed 
and running on your sytem and (2) that you have a clone of the pVACtools 
github repository.

The following ``docker run`` command publishes port 8000 in the docker 
image so that we can access services on this port from the host machine. 

It also mounts the location of the pVACtools git repo clone on the host machine
to a location within the docker image. "~/git/pVACtools/" in this example should 
be updated to the location of your clone of the repo on your system.

After launching the Docker image, we will install sphinx and its dependencies
using ``pip install``. 

Finally, we will use ``sphinx-autobuild`` to build the docs website. 
Once sphinx completes the build, you can view the docs website in a web browser 
on your system. As you make edits to the docs files (e.g. .rst files) the updates 
should automatically build and be viewable. 

.. code-block:: none

    docker pull griffithlab/pvactools
    docker run -p 8000:8000 -v ~/git/pVACtools/:/opt/git/pVACtools -it griffithlab/pvactools 
    pip install sphinx sphinx-autobuild sphinx-fontawesome sphinxcontrib-images sphinxjp.themes.basicstrap sphinxcontrib.programoutput
    cd /opt/git/pVACtools/docs/
    sphinx-autobuild --host 0.0.0.0 --port 8000 ./ ./_build/html

To view the updated docs site in your web browser go to: http://0.0.0.0:8000/. Once 
you are satisfied, submit a pull request as described above.

