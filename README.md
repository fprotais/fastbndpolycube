# fastbndpolycube

Take a surface triangular mesh along with a text file containing an integer (a flagging) per triangle (0 to 5, -> {+X, -X, +Y, -Y, +Z, -Z}) and returns a mesh where triangles are planar in their flagged dimension. The rest does its best to catch-up in a least-square way. 

The code is quite fast, and is a nice use case of [Disjoint-set](https://en.wikipedia.org/wiki/Disjoint-set_data_structure) to produce polycubes. 

I added [amgcl](https://github.com/ddemidov/amgcl), meaning that you need [boost](https://www.boost.org/) to compile the project. You can only compile main.cpp though, which doesn't require amgcl. Amgcl under the right tuning looks faster, while openNL is less troublesome to set-up. 

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
For the supported mesh formats, see [ultimaille](https://github.com/ssloy/ultimaille). 
