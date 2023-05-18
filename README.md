Custom physics in Godot
=======================

This is a quick overview for the "custom physics" module, which facilitates implementing
custom physics functionality. The module is simply an adapted version of the "Bullet" module
that ships with Godot.

You should hopefully only need to edit the `custom_dynamics_world.cpp` file, as well as any other files
that you might want to add yourself. In particular, you should **not** have to edit any of the files in the
`godot_module_files` folder.

Compiling Godot from source
---------------------------
The way native modules work in godot is to include them directly in its source code and rebuilt the project.

Before doing anything with this module, you need to make sure that you are able to compile Godot 3.2.1. from sources.
To do so, go to the Releases page for godot's github page and download the source code for version 3.2.1:

https://github.com/godotengine/godot/releases

Thankfully, building godot is very straightforward under normal circumstances. In short, you need only install
the `scons` build system and run

```
scons platform=windows target=release_debug -j4
```

in your terminal of choice for Windows, and

```
scons platform=x11 target=release_debug -j 4
```

for Linux platforms. Substitute platform=osx for Mac OSX. The above command uses 4 jobs for compilation.
Depending on the number of cores you have available, you may want to change this number.

--- Windows recommendations

scons may be installed through Python's package manager pip. However, you may find that this does not
necessarily give you a "scons" executable that you can call from your terminal. Our recommended approach
is to use Anaconda (or miniconda, a stripped down version) to manage your Python packages
(http://www.anaconda.com). After installing Anaconda (or miniconda), you can run the Anaconda shell and run

```
conda install scons
```

after which you should be able to run `scons` from inside your terminal (at least through the Anaconda shell).

Compilation
-----------

Copy the `customphysics` folder into the `modules` folder in the source distribution of Godot 3.2.1.
This effectively adds a `customphysics` module to Godot. When you (re)compile Godot again, files
from the `customphysics` module should be listed in the output. That is, you should only need to compile
Godot like usual (see up above).

Once you have compiled Godot with the `customphysics` module, you should run the binary
(e.g. `bin/godot.x11.opt.tools.64` for 64-bit Linux) and check that "CustomPhysics" is a
possible choice for the `Physics > 3d > Physics Engine` property (next to "Bullet", "GodotPhysics"
and "DEFAULT").

Faster re-compilation
---------------------

When running the above `scons` command, the build system is forced to check for changes in the entire
source tree of Godot, which may take some time. To speed up incremental rebuilds on some platforms,
you can choose to only rebuild the `customphysics` module as a shared library:

```
scons customphysics platform=x11 customphysics_dylib=yes target=release_debug -j 4
```

NOTE: This only works on Linux and possibly Mac OS (with some tweaks) platforms. In particular,
it does not seem to work on Windows.

Note that if you happen to make any changes of a more "major" nature, you might occassionally have to
do a "full" rebuild again.
