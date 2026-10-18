# Data — Higgs-Twitter (SNAP)

The review-revision scripts read the public SNAP Higgs-Twitter dataset.
Download both files into a local `data/` directory at the repository root:

```
mkdir -p data
wget https://snap.stanford.edu/data/higgs-twitter.zip
unzip higgs-twitter.zip higgs-activity_time.txt.gz higgs-retweet_network.edgelist.gz -d data/
# rename to the expected names:
mv data/higgs-activity_time.txt.gz data/higgs-activity.txt.gz
mv data/higgs-retweet_network.edgelist.gz data/higgs-retweet.edgelist.gz
```

Cite: A. De Domenico et al., "The anatomy of a scientific rumor,"
Scientific Reports 3, 2980 (2013). Dataset: SNAP, Stanford University.
The data files themselves are NOT redistributed in this repository.
