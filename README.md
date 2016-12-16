HyVis
=====
(c) 2015 NYU CAL

Install
-------

You'll need QT. Get it from [here](http://www.qt.io/download/).

Clone this repository into a directory of your choice. When the projects build, by default another folder will be created in the parent directory to the source code.
```bash
$ cd hyvisParent
$ git clone git@github.com:NYU-CAL/HyVis.git
```

Within the HyVis directory, open the file `HyVis.pro` and configure for your QT installation.

Installing Without the QT IDE
-----------------------------

The most important part of the QT install is the qmake command. If you have a working qmake, you can compile HyVis purely on the command line.  From the `src/` directory:

```bash
$ qmake -makefile HyVis.pro
$ make
```

The `qmake` command may generate some linking warnings, these are probably fine.

The result (on a Mac) will be a `HyVis.app/` which can be launched as an Application from Finder and contains a binary (for command line use) in `HyVis.app/Contents/MacOS/` which can be added to your PATH.

Switching between Disco and Jet Mode
------------------------------------

Eventually this will be automatic, for now you must edit the `filedata` variable declaration in `viewer.h` (currently on line 97).  For Disco this variable should be declared as `DiscoData filedata`, for Jet `JetData filedata`.



Building QT from source on OSX 10.6.8 
-------------------------------------

This was necessary to get a working QT distribution on an OSX 10.6.8 MacBook Pro which has various technical issues and not great C++11 support.  

Download the source of Qt 5.0.2. Compile as:
```bash
$ ./configure -prefix $PWD/qtbase -opensource -no-c++11 -debug-and-release -nomake examples -nomake demos
$ make -j 4
```

Add `qtbase/bin` to the `PATH`.

