func opra_svipc_init(void)
{
  extern am_master;
  require,"svipc.i";
  extern shmkey,semkey;
  shmkey = 0x00abacab;
  semkey = 0x0aabacab;
  randomize;
  shmkey = long(random()*1e7);
  semkey = long(random()*1e7);
  write,format="shmkey = %#x\n",shmkey;
  write,format="semkey = %#x\n",semkey;
  shm_init,shmkey,slots=20;
  sem_init,semkey,nums=20;
  shm_write,shmkey,"stop",&([0]);
  shm_write,shmkey,"next_stage",&([0]);
  shm_write,shmkey,"quit?",&([0]);
  shm_write,shmkey,"quit!",&([0]);
  shm_write,shmkey,"do plots",&([0]);
  shm_write,shmkey,"opra_structs",&([0]);
  if (fork()==0) { // I'm the child
    am_master = 0;
    // Start GUI
    extern fdpi; fdpi=dpi;
    images = [];
    require,"opra_gui.i";
    exit;
  } else {         // I'm the parent
    am_master = 1;
    extern quit;
    quit = opra_quit;
    set_idler,poll_quit;
    return;
  }
}


func write_data_to_fork(a,op,opp,mircube)
{
  shm_free,shmkey,"opra_structs";
  s = vsave("op",op,"opp",opp,"a",a,"mircube",mircube,\
      "lmfititer_pass",lmfititer_pass,"lmfit_itmax", \
      lmfit_itmax,"lmfititer",lmfititer);
  shm_write,shmkey,"opra_structs",&s;
}

func read_data_from_master(&a,&op,&opp,&mircube)
{
  extern lmfititer_pass,lmfit_itmax,lmfititer;
  restore,openb(shm_read(shmkey,"opra_structs"));
}
