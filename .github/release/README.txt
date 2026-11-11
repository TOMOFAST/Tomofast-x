Tomofast-x pre-built binaries
=============================

This is a standalone build of Tomofast-x. It needs no compiler, MPI or other
libraries: everything it uses is in this folder. It runs as one process and uses
OpenMP threads for parallelism (built with "make USE_MPI=0").

Runs on macOS 11 (Big Sur) or later, and on 64-bit Windows 10 or later.

For the source code, documentation and more examples, see
  https://github.com/TOMOFAST/Tomofast-x
  https://github.com/TOMOFAST/Tomofast-manual


Running
-------

macOS (Terminal, from this folder):

    export OMP_NUM_THREADS=4
    ./tomofastx -p parfiles/Parfile_mansf_slice.txt

Windows (Command Prompt, from this folder):

    set OMP_NUM_THREADS=4
    tomofastx.exe -p parfiles\Parfile_mansf_slice.txt

Windows (PowerShell, from this folder):

    $env:OMP_NUM_THREADS = 4
    .\tomofastx.exe -p parfiles\Parfile_mansf_slice.txt

Set OMP_NUM_THREADS to the number of CPU cores you want to use. If it is not set,
all cores are used.

The example writes its results to output/mansf_slice/. When it finishes, the log
ends with "THE END." Open output/mansf_slice/Paraview/grav_final_model3D_full.vtk
in ParaView to see the final model.

Paths in a parameter file are relative to the folder you run the program from.

If you get a crash on a large model, increase the thread stack size, for example
OMP_STACKSIZE=512M (and on macOS also run "ulimit -s unlimited" first).


macOS security warning
----------------------

These binaries are not notarised by Apple. If macOS says the program "cannot be
opened" or is "damaged", remove the download quarantine flag once, from this folder:

    xattr -dr com.apple.quarantine .

Use the "macos-arm64" download on Apple Silicon Macs (M1 and later) and the
"macos-x86_64" download on Intel Macs.


Contents
--------

    tomofastx / tomofastx.exe   the program
    lib/                        runtime libraries (macOS only)
    parfiles/                   example parameter file
    data/                       example input data
    LICENSE                     MIT license
