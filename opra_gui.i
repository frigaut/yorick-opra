/* pykex.i
   main function to call the pygtk GUI to pykex.
   syntax: yorick -i pykex.i
   Authors: F.Rigaut, May 2007
*/

require,"pyk.i";

//
// this function is called by python after the main GUI window has been mapped:
//
func opra_win_init(xid,..)  // main display window.
{
  xids = [xid];
  while (nid=next_arg()) grow,xids,nid;
  xids;

  //#################################
  // create graphic windows
  //#################################
  winnum = opp.winnum;
  for (n=1;n<=opp.npos;n++) {
    winn = winnum+n-1;
    if (window_exists(winn)) { winkill,winn;}
    window,winn,wait=0,style="opra.gs",dpi=fdpi,width=0,height=0,parent=xids(n);
    if (pal) palette,pal;
    for (i=1;i<=3;i++) { plsys,i; limits,square=1; }
    if (opp.modes_type=="yao") {
      if (window_exists(winnum+opp.npos)) winkill,winnum+opp.npos;
      window,winnum+opp.npos,style="nobox.gs",width=0,height=0,\
        dpi=fdpi,wait=0,parent=xids(opp.npos+1);
      if (window_exists(winnum+opp.npos)) winkill,winnum+opp.npos+1;
      window,winnum+opp.npos+1,width=0,height=0,\
        dpi=fdpi,wait=0,parent=xids(opp.npos+2);
      plp,wfs.gspos(2,),wfs.gspos(1,),symbol='*',color="red";
      // for (i=1;i<=numberof(wfs.gspos(1,));i++) \
      for (i=1;i<=nwfs;i++) \
        plt,swrite(format="  %d",i),wfs.gspos(1,i),wfs.gspos(2,i),tosys=1;
      plmargin;
      redraw;
    }
  }
  sem_give,semkey,0;
  pyk,"done_init = 1";
}

//
// wrapper functions
//

func next_stage(void)
{
  shm_write,shmkey,"next_stage",&([1]);
}

func opra_stop(void)
{
  shm_write,shmkey,"stop",&([1]);
}

func opra_quit(void)
{
  shm_write,shmkey,"stop",&([1]);
  shm_write,shmkey,"quit?",&([1]);
  write,format="%s\n","Quit queued";
}


func opra_plots(void)
{
  read_data_from_master,a,op,opp,mircube;
  for (n=1;n<=opp.npos;n++) {
    opra_info_and_plots,a,opp,op(,n),yao=(opp.modes_type=="yao"),star=n,noprint=1;
  }
  if (opp.modes_type=="yao") {
    window,opp.winnum+opp.npos;
    tv,mircube(,*),square=1;
  }
}

// restore original quit behavior:
quit = original_quit;

func opra_gui_poll(void)
{
  // check if we have to quit:
  if (shm_read(shmkey,"quit?")(1)) {
    write,format="%s\n","Fork quitting";
    pyk,"on_quit_requested()";
    shm_write,shmkey,"quit!",&([1]);
    quit;
  }

  // check if we have to do plots:
  if (shm_read(shmkey,"do plots")(1)) {
    shm_write,shmkey,"do plots",&([0]);
    status = opra_plots();
  }

  after,0.05,opra_gui_poll;
}
opra_gui_poll;

func fork_quit(void)
{
  _pyk_proc = [];
  quit;
}
after_error = fork_quit;

// main stuff:
pyk_debug=0;

// build command
python_exec = find_in_path("../python/opra_gui.py",takefirst=1);
path_to_glade = dirname(find_in_path("../glade/opra_gui.glade",takefirst=1))+"/";
pyk_cmd=[python_exec,path_to_glade,swrite(format="%d",opp.npos), \
         swrite(format="%d",fdpi),swrite(format="%d",(opp.modes_type=="yao"))];

// spawn it and attach to _pyk_callback (see pyk.i):
_pyk_proc = spawn(pyk_cmd, _pyk_callback);

