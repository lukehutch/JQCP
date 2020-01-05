# JQCP

Java port of Liu and Theobold's [Quaternion Characteristic Polynomial (QCP)](https://theobald.brandeis.edu/qcp/) method for determining the minimum RMSD between two structures and for determining the optimal least-squares rotation matrix.

This is a direct port from the C code to Java, with the exception that the out parameters were moved to the end of the parameter list for each method.

This repository also contains a second port of the code to the [JOML](https://github.com/JOML-CI/JOML) Java linear algebra library.
