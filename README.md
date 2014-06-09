oencast-bio
============

A biological data mining project which is part of a Masters degree at the University of Edinburgh.

-----------------------------------

There are three parts to this project.
Firstly, we're trying to replicate and build on previous work in the area of protein interaction prediction.
Secondly, applying various different supervised classification algorithms to this data as a machine learning exercise.
Thirdly, Using this classification algorithm to build a weighted graph to see if the weighted graph can improve the performance of a community detection algorithm used at the University of Edinburgh.

For detailed information about the proposed project, look at the [proposal][].

## Progress

To look at detailed information about the progress of this project look at the [wiki][opencast wiki].

At this point we are trying to build feature vectors for the [full protein list][proteinlist].
We have a [long list of data sources][sourcelist] which we are working through trying to quickly extract features.

To do list:

1. ~~Full list of proteins of interest~~
2. Feature vectors
3. Supervised classification code
4. Community detection analysis

## Combining data sources

Protein interaction prediction is a binary prediction task so involves gathering features from many different biological data sources and constructing positive and negative feature vectors.
Picking the different data sources from what's available depends on knowledge of the biology and what might work well as a predictor.
Luckily, many papers have been published on this subject.

[opencastwiki]: 
[proposal]: 
[proteinlist]: 
[sourcelist]:
