/* pyk.i
 * 2 way communication interface to python (useful for GUIs)
 * 
 * $Id: pyk.i,v 1.4 2009-04-01 06:43:07 paumard Exp $
 * 
 * Author: Francois Rigaut
 * Written May 2007
 * last revision/addition: 2007
 *
 * Copyright (c) 2007, Francois Rigaut
 *
 * This program is free software; you can redistribute it and/or  modify it
 * under the terms of the GNU General Public License  as  published  by the
 * Free Software Foundation; either version 2 of the License,  or  (at your
 * option) any later version.
 *
 * This program is distributed in the hope  that  it  will  be  useful, but
 * WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
 * MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
 * General Public License for more details (to receive a  copy  of  the GNU
 * General Public License, write to the Free Software Foundation, Inc., 675
 * Mass Ave, Cambridge, MA 02139, USA).
 *
   Inspired from tyk.i (tcltk interface)
*/

Y_PYTHON = get_env("Y_PYTHON");
if (noneof(Y_PYTHON))                                                   \
  Y_PYTHON="./:"+Y_USER+":"+pathform(_(Y_USER,Y_SITES,Y_SITE)+"python/");


func pyk(py_command)
/* DOCUMENT pyk, py_command
 *       or value = pyk(py_command)
 *
 *   send PY_COMMAND to python front-end.  If the wish front-end is not
 *   already running, pyk starts it.
 *   In the second form, pyk suspends the yorick interpreted caller
 *   until python invokes a pyk_resume command returing a value.
 *
 * SEE ALSO: pyk_debug, pyk_resume
 */
{
  require,"pathfun.i";
  
  if (is_void(_pyk_proc)) {
    // start python through pyk.py (includes python initialization):
    pykpy = find_in_path("pyk.py",takefirst=1);
    if (pykpy==[]) pykpy = find_in_path("../python/pyk.py",takefirst=1);
    if (pykpy==[]) error,"Can't find pyk.py";
    if (pyk_debug) write,format="Init python, loading pyk.py, found at "+pykpy;
    _pyk_proc = spawn(pykpy, _pyk_callback);
  }
  if (is_void(py_command)) return;

  if (pyk_debug) write,format="to python: %s\n",py_command;

  block = !am_subroutine();
  /* if blocking, first wait for python to synchronize */
  if (block) {
    if (_pyk_sync) {
      _pyk_sync = 0;
      error, "already blocked waiting for python";
    }
    _pyk_sync = 1;
    _pyk_proc, "pyk_sync()\n";
    suspend;
    if (_pyk_sync) {
      _pyk_sync = 0;
      error, "python front-end shut down unexpectedly";
    }
  }

  /* send the command to python */
  if (strpart(py_command,0:0) != "\n") py_command += "\n";
  _pyk_proc, py_command;

  /* block for response if requested
   * note that python knows from previous pyk_sync command that it needs
   * to respond with a pyk_resume command to resume yorick
   * - it is free to respond with multiple lines if necessary to
   *   build a complicated return value
   */
  if (block) {
    suspend;
    value = _pyk_value;
    _pyk_value = [];
    return value;
  }
}

local pyk_debug;
/* DOCUMENT pyk_debug = 1
 *       or pyk_debug = []
 *
 *   If pyk_debug is non-nil and non-zero, print message traffic to and
 *   from wish.  This is useful for debugging py-yorick interaction.
 *
 * SEE ALSO: pyk
 */


func pyk_import(filename)
/* DOCUMENT pyk_import, filename
 *
 *   source FILENAME containing python program into python front-end.
 *   If python front-end is not already running, pyk_source starts it.
 *
 * SEE ALSO: pyk
 */
{
  if (strpart(filename,1:1)=="/") { // absolute path
    dirs="";
  } else {
    dirs = Y_PYTHON;
    //    dirs = [".","~/yorick","~/Yorick",Y_SITE+"i",Y_SITE+"contrib",
    //            Y_SITE+"i0",Y_HOME+"lib"]+"/";
  }

  for (i=1;i<=numberof(dirs);i++) {
    if (open(dirs(i)+filename,"r",1)) {
      pyk, "import " + dirs(i)+filename;
      if (pyk_debug) {
        write,format="Loading %s, found in %s\n",filename,dirs(i);
      }
      return;
    }
  }
  error,"Can't find "+filename;
}


func pyk_resume(value)
/* DOCUMENT pyk_resume value
 *
 *   This function must be invoked by the python front-end; it should
 *   never be called by yorick programs.
 *
 *   When python receives a pyk_sync command, it must eventually
 *   respond with a pyk_resume in order to unblock the pyk function.
 *
 * SEE ALSO: pyk
 */
{
  extern _pyk_value;
  _pyk_value = value;
  resume;
}

func pyk_set(&v1,x1,&v2,x2,&v3,x3,&v4,x4,&v5,x5,&v6,x6,&v7,x7,&v8,x8)
/* DOCUMENT pyk_set var1 val1 var2 val2 ...
 *
 *   This function is designed to be invoked by the python front-end;
 *   it is not useful for yorick programs.
 *
 *   Equivalent to
 *     var1=val1; var2=val2; ...
 *   Handles at most 8 var/val pairs.
 *   As a special case, if given an odd number of arguments, pyk_set
 *   sets the final var to [], e.g.-
 *     pyk_set var1 12.34 var2
 *   is equivalent to
 *     var1=12.34; var2=[];
 *
 * SEE ALSO: pyk
 */
{
  v1 = x1;
  if (is_void(x1)) return; else v2 = x2;
  if (is_void(x2)) return; else v3 = x3;
  if (is_void(x3)) return; else v4 = x4;
  if (is_void(x4)) return; else v5 = x5;
  if (is_void(x5)) return; else v6 = x6;
  if (is_void(x6)) return; else v7 = x7;
  if (is_void(x7)) return; else v8 = x8;
}

_pyk_sync = 0;
_pyk_linebuf = string(0);

func _pyk_callback(line)
{
  extern _pyk_proc;
  if (!line) {
    _pyk_proc = [];
    if (_pyk_sync) resume;
    if (pyk_debug) write, "from python -> <python terminated>";
    return;
  }
  /* must be prepared for python output to dribble back a fraction of
   * a line at a time, or multiple lines at a time
   * _pyk_linebuf holds the most recent incomplete line,
   *   assuming the the remainder will arrive in future callbacks
   */
  _pyk_linebuf += line;
  selist = strword(_pyk_linebuf, "\n", 256);
  line = strpart(_pyk_linebuf, selist);
  line = line(where(line));
  n = numberof(line);
  if (n && selist(2*n)==strlen(_pyk_linebuf)) {
    /* final character of input not \n, store fragment in _pyk_linebuf */
    _pyk_linebuf = line(0);
    if (n==1) return;
    line = line(1:-1);
  } else {
    _pyk_linebuf = string(0);
  }
  strtrim, line;
  line = line(where(strlen(line)));

  if (pyk_debug) write, "from python:", line;

  /* record whether final line is a sync message, remove it */
  mask =(line != "-s+y-n+c-+p-y+k-");
  synchronized = !allof(mask);
  line = line(where(mask));
  nofline = numberof(line);

  /* parse and execute yorick command lines */
  for (i=1 ; i<=nofline ; i++) funcdef(line(i));

  /* permit pyk function to continue after successful synchronization */
  if (synchronized && _pyk_sync) {
    _pyk_sync = 0;
    resume;
  }
}
