opencast-bio
============

A biological data mining project as part of a Masters degree at the University of Edinburgh.

-----------------------------------

There are three parts to this project.
Firstly, we're trying to replicate and build on previous work in the area of protein interaction prediction.
Secondly, applying various different supervised classification algorithms to this data as a machine learning exercise.
Thirdly, Using this classification algorithm to build a weighted graph to see if the weighted graph can improve the performance of a community detection algorithm used at the University of Edinburgh.

For detailed information about the proposed project, look at the [proposal][].

# Using the code

If you would like to repeat the weighting of the edges and community detection(and you have access to the data required), from start to finish:

1. Run the [Training and Test set generation notebook][traintest]
2. Run the [Classifier Training notebook][classtrain]
3. Run the [Community Detection notebook][comdet]

This assumes you have access to the extracted features.
If you would like to repeat those, the notebooks required are referenced in the report's appendices.

## ocbio.extract

The dedicated code that was written as part of this project.
This is mostly useful for accessing the various pre-extracted features, if you have access to the data.
Instructions on how to use it in general can be found in the [ocbio.extract usage notes][ocbio].
The other pieces of Python code in ocbio are not intended to be used directly - they are simply imported in `extract.py`.

# Required packages

All that's required for most of this code is pylab and Scikit-learn. For a full list of all possible requirements and their versions look at the notebook on [required packages][reqnotes]

# Progress

This project has completed.
The report is available in the repository.

[opencastwiki]: https://github.com/ggray1729/opencast-bio/wiki
[proposal]: https://github.com/ggray1729/opencast-bio/raw/master/proposal/proposal4.pdf
[proteinlist]: https://github.com/ggray1729/opencast-bio/wiki/Full-protein-list
[sourcelist]: https://github.com/ggray1729/opencast-bio/wiki/Feature-extraction
[reqnotes]: http://nbviewer.ipython.org/github/ggray1729/opencast-bio/blob/master/notebooks/Required%20packages.ipynb
[traintest]: http://nbviewer.ipython.org/github/ggray1729/opencast-bio/blob/master/notebooks/Training%20and%20test%20featureset%20generation.ipynb
[classtrain]: http://nbviewer.ipython.org/github/ggray1729/opencast-bio/blob/master/notebooks/Classifier%20Training.ipynb
[comdet]: http://nbviewer.ipython.org/github/ggray1729/opencast-bio/blob/master/notebooks/Community%20Detection.ipynb
[ocbio]: https://github.com/ggray1729/opencast-bio/blob/master/notebooks/ocbio.extract%20usage%20notes.ipynb
