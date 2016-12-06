.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Usage in Publications
~~~~~~~~~~~~~~~~~~~~~

If you use CliMT to perform research, your publication is a valuable resource
for others looking to learn the ways they can leverage CliMT's capabilities.
If you have used CliMT in a publication, please let us know so we can add it to
the list.

Presenting CliMT to Others
~~~~~~~~~~~~~~~~~~~~~~~~~~

CliMT is meant to be an accessible, community-driven model. You can help the
community of users grow and be more effective in many ways, such as:

* Running a workshop
* Offering to be a resource for others to ask questions
* Presenting research that uses CliMT

If you or someone you know is contributing to the CliMT community by presenting
it or assisting others with the model, please let us know so we can add that
person to the contributors list.

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/CliMT/climt/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
and "help wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

CliMT could always use more documentation. You could:

* Clean up or add to the official CliMT docs and docstrings.
* Write useful and clear examples that are missing from the examples folder.
* Create a Jupyter notebook that uses CliMT and share it with others.
* Prepare reproducible model scripts to distribute with a paper using CliMT.
* Anything else that communicates useful information about CliMT.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/CliMT/climt/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `climt` for local development.

1. Fork the `climt` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/climt.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv climt
    $ cd climt/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the tests, including testing other Python versions with tox::

    $ flake8 climt tests
    $ python setup.py test or py.test
    $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 2.7, 3.4 and 3.5. Check
   https://travis-ci.org/CliMT/climt/pull_requests
   and make sure that the tests pass for all supported Python versions.

Style
-----

In the CliMT code, we follow PEP 8 style guidelines (tested by flake8). You can
test style by running "tox -e flake8" from the root directory of the repository.
There are some exceptions to PEP 8:

* All lines should be shorter than 80 characters. However, lines
  longer than this are permissible if this increases readability (particularly
  for lines representing complicated equations).
* Space should be assigned around arithmetic operators in a way that maximizes
  readability. For some cases, this may mean not including whitespace around
  certain operations to make the separation of terms clearer,
  e.g. "Cp*T + g*z + Lv*q".
* While state dictionary keys are full and verbose, within components they may
  be assigned to shorter names if it makes the code clearer.
* We can take advantage of known scientific abbreviations for quantities within
  components (e.g. "T" for "air_temperature") even thought they do not follow
  pothole_case.

Tips
----

To run a subset of tests::

$ py.test tests.test_timestepping

