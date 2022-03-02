# OPRA


## WHAT IS IT?
OPRA (OTF-based Phase Retrieval Analysis) is a phase retrieval
tool that uses phase diversity. That is, you feed it a few images
-in between which some *phase diversity* has been introduced,
namely some focus- and it returns the corresponding phase.
OPRA specificity includes:

* Instead of minimizing the distance model/measurement in the
  image plane, it does it in the OTF plane.
* It can use several basis of modes to reconstruct the phase:
  Disk harmonic (the default), Karhuenen-Loeve or Zernike.
* It can also work within the yao framework, which allow the user
  to define more complex conditions/systems. In particular, from
  v1.5, one can use a set of images to reconstruct the phase
  aberration tomography using a set of modes at various altitudes.

OPRA was developed within the framework on the GEMINI MCAO system.


## AUTHORS
OPRA, Francois Rigaut (c) 2010-2011

Authors:
Francois Rigaut    frigaut@gemini.edu
Benoit Neichel     bneichel@gemini.edu
Damien Gratadour   damien.gratadour@obspm.fr


## LICENSE
This program is free software; you can redistribute it and/or  modify it
under the terms of the GNU General Public License  as  published  by the
Free Software Foundation; either version 2 of the License,  or  (at your
option) any later version.

This program is distributed in the hope  that  it  will  be  useful, but
WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
General Public License for more details (to receive a  copy  of  the GNU
General Public License, write to the Free Software Foundation, Inc., 675
Mass Ave, Cambridge, MA 02139, USA).


## INSTALLATION

OPRA is just a set of interpreted yorick files. You can easily install from
source (it's all source anyway) using the usual:

    yorick -batch make.i
    [sudo] make install

Note that if you want to use the yao framework (in particular, the tomography
mode requires it), you will need yao installed.


## TEST IT / USE IT

There are a few test cases in the examples directory. Look at the README
file in this directory for more details.

    yorick -i opra_test1.i

* opra_test1.i Is the simplest test case. Using disk harmonics.
* opra_test2.i uses yao, so you'll need yao installed. Same images
  and same phase as opra_test1.i.

For the test cases 3,4, and 5 (tomographic cases), you need to
generate the data set before running the analysis:

    yorick -i gen_test3_images.i
    yorick -i opra_test3.i

(obviously replace 3 by 4 or 5 to run the other cases).
test3 is with 2DMs. test4 is with 3DMs. test5 adds noise to test4.

* opra_test6.i uses images from the Gemini NIRI imager.