# fastbndpolycube

Take a surface triangular mesh along with a text file containing a integer (a flagging) per triangle (0 to 5, -> {+X, -X, +Y, -Y, +Z, -Z}) and returns a mesh where triangular are planar in their flagged dimension. The rest does its best to catch-up in a least-square way. 

The code is quite fast, and is a nice use case of Disjoint-set to produce polycubes. 

# Use CMake to build the project:
```sh
git clone --recurse-submodules https://github.com/fprotais/fastbndpolycube &&
cd fastbndpolycube &&
mkdir build &&
cd build &&
cmake .. &&
make -j 
```

# Running the code :

```sh
./fastpolycube ../boundary.obj ../labeling.txt result.obj
```
For the supported mesh formats, see ultimaille. 