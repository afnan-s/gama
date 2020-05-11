# gama: a Genetic Approach to Maximize Clustering Criteria in R

We presented an R package to perform hard partitional clustering guided by an user-specified cluster validation criterion. The algorithm obtains high cluster validation indices when applied to datasets who contains superellipsoid clusters. The algorithm is capable of estimate the number of partitions for a given dataset by an automatic inference of the elbow in WCSSE graph or by using a broad search in 24 cluster validation criteria. The package brings six different built-in datasets for experimentation, two of them are in-house datasets collected from real execution of distributed machine learning algorithms on Spark clusters. The others are well-known datasets used in the benchmark of clustering problems.

This version contains modifications of the original:
- Random solution generation more suited to large sparse datasets.
- Option to recieve the distance matrix (for efficiency)