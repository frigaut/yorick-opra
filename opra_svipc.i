require,"svipc.i";
require,"opra_pyk.i";

func opra_svipc_init(dpi)
{
  extern pyk_debug, _yorick_gui_proc;
  extern am_master;
  extern shmkey,semkey;

  randomize;
  shmkey = long(random()*1e7);
  semkey = long(random()*1e7);
  // shmkey = 0x0badcafe;
  // semkey = 0x0badbeef;
  write,format="shmkey = %#x\n",shmkey;
  write,format="semkey = %#x\n",semkey;

  // init SHM and SEM stuff:
  shm_init,shmkey,slots=20;
  sem_init,semkey,nums=20;
  shm_write,shmkey,"stop",&([0]);
  shm_write,shmkey,"next_stage",&([0]);
  shm_write,shmkey,"quit?",&([0]);
  shm_write,shmkey,"quit!",&([0]);
  shm_write,shmkey,"do plots",&([0]);
  shm_write,shmkey,"opra_structs",&([0]);

  // spawn GUI yorick process
  // build command
  pyk_cmd=["yorick","-q","-i","opra_gui.i",swrite(format="%d",dpi),
    swrite(format="%d",shmkey),swrite(format="%d",semkey)]
  // spawn it and attach to _pyk_callback (see pyk.i):
    _yorick_gui_proc = spawn(pyk_cmd, _pyk_callback);

  extern quit;
  quit = opra_quit;
  set_idler,poll_quit;
  return;
}


func write_yao_struct_to_shm(wfs,dm)
{
  shm_free,shmkey,"yao_structs";
  s = vsave("wfs",wfs,"dm",dm);
  shm_write,shmkey,"yao_structs",&s;
}


func read_yao_struct_from_shm(&wfs,&dm)
{
  restore,openb(shm_read(shmkey,"yao_structs"));
}


func write_opra_struct_to_shm(a,op,opp,mircube)
{
  if (mircube==[]) mircube=[0];
  shm_free,shmkey,"opra_structs";
  s = vsave("op",op,"opp",opp,"a",a,"mircube",mircube,\
    "lmfititer_pass",lmfititer_pass,"lmfit_itmax", \
    lmfit_itmax,"lmfititer",lmfititer);
  shm_write,shmkey,"opra_structs",&s;
}


func read_opra_struct_from_shm(&a,&op,&opp,&mircube)
{
  extern lmfititer_pass,lmfit_itmax,lmfititer;
  restore,openb(shm_read(shmkey,"opra_structs"));
}
