# these values filled in by    yorick -batch make.i
Y_MAKEDIR=/usr/lib/yorick/2.2
Y_EXE=/usr/lib/yorick/2.2/bin/yorick
Y_EXE_PKGS=
Y_EXE_HOME=/usr/lib/yorick/2.2
Y_EXE_SITE=/usr/share/yorick/2.2

#
# !! THIS IS NOT A PLUGIN !!
# This is a package made of several interpreted
# include file. This makefile is just used to install,
# uninstall it or build the distribution tar file.

# ------------------------------------------------ macros for this package

# used for distribution
PKG_NAME = opra
# include files for this package
PKG_I=opra.i opra_lmfit.i opra_utils.i opra_caller.i opra_libkl.i opra_libdh.i opra_structs.i opra_gui.i

# autoload file for this package, if any
PKG_I_START =

# override macros Makepkg sets for rules and other macros
# Y_HOME and Y_SITE in Make.cfg may not be correct (e.g.- relocatable)
Y_HOME=$(Y_EXE_HOME)
Y_SITE=$(Y_EXE_SITE)

# include $(Y_MAKEDIR)/Make.cfg
DEST_Y_SITE=$(DESTDIR)$(Y_SITE)
DEST_Y_HOME=$(DESTDIR)$(Y_HOME)

# ------------------------------------- targets and rules for this package

build:
	@echo "Nothing to build. This is not a plugin"
	@echo "other targets: install, uninstall, clean"
	@echo "for maintainers: package, distpkg"

clean:
	-rm -rf pkg *~

install:
	mkdir -p $(DEST_Y_SITE)/i
	mkdir -p $(DEST_Y_SITE)/g
	mkdir -p $(DEST_Y_SITE)/python
	mkdir -p $(DEST_Y_SITE)/glade
	cp -p $(PKG_I) $(DEST_Y_SITE)/i/
	cp -p opra.gs $(DEST_Y_SITE)/g/
	cp -p opra_gui.py $(DEST_Y_SITE)/python/.
	cp -p opra_gui.glade $(DEST_Y_SITE)/glade/.

uninstall:
	-cd $(DEST_Y_SITE)/i; rm $(PKG_I)
	-cd $(DEST_Y_SITE)/g; rm opra.gs

