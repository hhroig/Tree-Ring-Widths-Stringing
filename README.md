# Tree-Ring-Widths-Stringing
Data illustration: stringing data on tree-ring-widths.
Author: Harold Antonio Hernández Roig
June, 2019

## General Description:

We devote this work to mimic Section 4. Data Illustrations / 4.1 Stringing of Tree-Ring Widhts in [1].

We focus our attention on the concept of \emph{stringing} applied to \emph{high-dimensional (HD) data}, the $n \ll p$ scenario. Moreover, classical tools as Functional Principal Component Analysis (FPCA), recovers an important role in the analysis and further summary of the stringed HD data.

Stringing is based on the assumption that the predictors are correlated or have some sort of neighborhood relationship. Said so, a smooth stochastic process could be constructed from the data. Another view is that we actually observe predictors with a randomly permuted order, but after reordering (stringing), we observe realizations of a smooth function.

Through stringing, we can map the HD predictor vectors to the infinite-dimensional function space. Once we have this, it is possible to analyze the resulting smooth random functions with FPCA. Particularly, we address the usage of R package: fdapace [2], which implements Principal Analysis through Conditional Expectation (PACE) [3]. This package also provides an implementation of stringing functionality.

## Data

Data is freely available at: ftp://ftp.ncdc.noaa.gov/pub/data/paleo/treering/measurements/northamerica/usa/ca645.rwl

It consists on tree-ring widhts measurements taken from period 1932-1976, at Mary Ranch, Santa Clara, California. Main contributed by D.W. Stahle and R.D. Griffin.

We compare measurments to rainy season precipitation (defined as precipitation from previous December to April of the current year). Data was obtained at WestMap: Climate Analysis & Mapping Toolbox (check https://cefa.dri.edu/Westmap/Westmap_home.php).

## References

[1] K Chen, H G Müller, and J Wang. Stringing High-Dimensional Data for Functional Analysis. Journal of the American Statistical Association, 106(493):275–284, 2011.

[2] Xiongtao Dai, Pantelis Z. Hadjipantelis, Kynghee Han, Hao Ji, Shu-Chin Lin, Hans-Georg
Mueller, and Jane-Ling Wang. fdapace: Functional Data Analysis and Empirical Dynamics, 2018. R package version 0.4.0.

[3] Fang Yao, Hans Georg Müller, and Jane Ling Wang. Functional data analysis for sparse
longitudinal data. Journal of the American Statistical Association, 100(470):577–590, 2005.
