#!/usr/bin/env python
# authors: Matthieu Bec, Francois Rigaut
# opra.py

import gtk
import gtk.glade
import sys
import gobject
import os, fcntl, errno

class opra:

   def destroy(self, wdg, data=None):
      self.py2yo('opra_quit')
#      gtk.main_quit()

   def __init__(self,opratop,ntab,dpi,is_yao):
      self.opratop = opratop

      self.glade = gtk.glade.XML(os.path.join(self.opratop,'opra_gui.glade'))
      self.window = self.glade.get_widget('window1')
      # handle destroy event
      if (self.window):
         self.window.connect('destroy', self.destroy)
      # autoconnect to callback functions
      self.glade.signal_autoconnect(self)

      # set stdin non blocking, this will prevent readline to block
      fd = sys.stdin.fileno()
      flags = fcntl.fcntl(fd, fcntl.F_GETFL)
      fcntl.fcntl(fd, fcntl.F_SETFL, flags | os.O_NONBLOCK)

      # add stdin to the event loop (yorick input pipe by spawn)
      gobject.io_add_watch(sys.stdin,gobject.IO_IN|gobject.IO_HUP,self.yo2py,None)

      # Create a new notebook, place the position of the tabs
      table = gtk.Table(3,6,False)
      self.glade.get_widget('tabs').add(table)
      self.notebook = gtk.Notebook()
      self.notebook.set_tab_pos(gtk.POS_TOP)
      table.attach(self.notebook, 0,6,0,1)
      self.notebook.show()
      self.show_tabs = True
      self.show_border = True

      swid = ""
      for i in range(ntab):
         darea = gtk.DrawingArea()
         darea.set_size_request(6*dpi,6*dpi)
         bufferl = " %d " % (i+1)
         darea.show()
         label = gtk.Label(bufferl)
         self.notebook.append_page(darea, label)
         # table.show()
         swid = swid+str(darea.window.xid)+" ";

      if is_yao:
         darea = gtk.DrawingArea()
         darea.set_size_request(6*dpi,6*dpi)
         bufferl = "DM"
         darea.show()
         label = gtk.Label(bufferl)
         self.notebook.append_page(darea, label)
         swid = swid+str(darea.window.xid)+" ";

         darea = gtk.DrawingArea()
         darea.set_size_request(6*dpi,6*dpi)
         bufferl = "Geometry"
         darea.show()
         label = gtk.Label(bufferl)
         self.notebook.append_page(darea, label)
         swid = swid+str(darea.window.xid)+" ";

      table.show()

      self.py2yo('opra_win_init %s' % swid)
      # run
      gtk.main()

   done_init = 0

   def area_expose_cb(self, area, event):
      # self.style = self.area.get_style()
      # self.gc = self.style.fg_gc[gtk.STATE_NORMAL]
      self.draw_pixmap(0, 0)
      return True

   def on_about_activate(self,wdg):
      dialog = self.glade.get_widget('aboutdialog')
      dialog.run()
      dialog.hide()

   def on_quit_activate(self,*args):
      self.py2yo('opra_quit')
      raise SystemExit

   def on_tabs_prev_clicked(self,wdg):
      self.notebook.prev_page()

   def on_tabs_next_clicked(self,wdg):
      self.notebook.next_page()

   #
   # Yorick to Python Wrapper Functions
   #

   def on_unzoom_clicked(self,wdg):
      self.py2yo('opra_stop')

   def on_next_stage_clicked(self,wdg):
      self.py2yo('next_stage')

   def on_window1_map_event(self,wdg,*args):
      pass
      # drawingarea = self.glade.get_widget('drawingarea1')
      # mwid1 = drawingarea.window.xid;
      # self.py2yo('opra_win_init %d' % mwid1)

   def on_quit_requested(self):
      raise SystemExit

   #
   # minimal wrapper for yorick/python communication
   #
   def yo2py_flush(self):
      sys.stdin.flush()

   def py2yo(self,msg):
      # sends string command to yorick's eval
      sys.stdout.write(msg+'\n')
      sys.stdout.flush()

   def yo2py(self,cb_condition,*args):
      if cb_condition == gobject.IO_HUP:
         raise SystemExit, "lost pipe to yorick"
      # handles string command from yorick
      # note: inidividual message needs to end with \n for proper ungarbling
      while 1:
         try:
            msg = sys.stdin.readline()
            # sys.stderr.write('received from python: %s\n' % msg)
            if msg=="":
               return False
            msg = "self."+msg
            exec(msg)
         except IOError, e:
            if e.errno == errno.EAGAIN:
               # the pipe's empty, good
               break
            # else bomb out
            raise SystemExit, "yo2py unexpected IOError:" + str(e)
         except Exception, ee:
            raise SystemExit, "yo2py unexpected Exception:" + str(ee)
      return True


if len(sys.argv) != 5:
   print 'Usage: opra.py path_to_glade ntab dpi is_yao'
   raise SystemExit

opratop = str(sys.argv[1])
ntab = long(sys.argv[2])
dpi  = long(sys.argv[3])
is_yao  = long(sys.argv[4])
top  = opra(opratop,ntab,dpi,is_yao)
