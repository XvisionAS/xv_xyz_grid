# What is it ?
This tool takes a XYZ file and transform it into an OBJ file.

# How does it works ?
* It reads the XYZ file into a list of points
* It then remove duplicate
* It then use "triangle" to generate a mesh from the points
* If you decide to simplify the mesh, it will then generate a new set of points by ray tracing to the just generated mesh, generating a grid of points, and then use "triangle" to generate a new mesh from this set of points.
* It generate the OBJ file.

# How to build ?
This tool use premake to generate the necessary files for compiling. Just go to the folder "Platforms", peek the right "premake5" executable, and :
* for windows : premake5 vs2015
* for osx : premake5.osx gmake
* for linux : premake5.linux gmake

Then go to the "build" folder, and launch either **make** or **visual studio** from the generated files.

# Futur ?
This tool is a simple proof of concept for now ( needed to be able to display a heightmap readed from XYZ ). At some point I will need this tool to generate a complete tile set with LOD, and the export should be a TIFF or a custom binary format ( TBD ).

# Can I use it ?
Of course, just drop me a line. Note that the XYZ *parser* is far from handling all the cases ( values must be separated by spaces or tabs, no header ). Other than that, it is also not very memory efficient, so big datasets might not work ( biggest one I tried was about  500megs ).
