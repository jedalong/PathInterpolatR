---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# PathInterpolatR

The package PathInterpolatR provides access to some simple functions for path interpolation (e.g., from GPS tracking data). Most notably, it is associated with the methods described in:

Long, J.A. (2016) Kinematic interpolation of movement data. International Journal of
Geographical Information Science. 30(5): 854-868.

## Installation

You can install PathInterpolatR from github with:


```r
# install.packages("devtools")
devtools::install_github("jedalong/PathInterpolatR")
```

## Quick Demo

Below a quick demo shows how the package can be used with the kinematic interpolation function.


```r
library(PathInterpolatR)
data(contrived)
xyt <- contrived
###times for interpolation
t.slice <- c(1.5,2,2.5,3,3.5,4,4.5,5,5.5)
### add velocities for kinematic analysis
xytvv <- cbind(xyt[2:3,],rbind(c(0,3),c(3,0)))
a <- kinematic(xytvv,t.slice)
plot(xyt[,1],xyt[,2],pch=20,asp=1)
points(a[,1],a[,2])
```

![plot of chunk example](README-example-1.png)

--- End ---
