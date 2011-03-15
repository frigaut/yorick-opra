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
  shm_init,shmkey,slots=20;
  sem_init,semkey,nums=20;
  shm_write,shmkey,"stop",&([0]);
  shm_write,shmkey,"next_stage",&([0]);
  shm_write,shmkey,"quit?",&([0]);
  shm_write,shmkey,"quit!",&([0]);
  shm_write,shmkey,"do plots",&([0]);
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
  opv =  array(char,sizeof(op));
  mem_copy,mem_base(opv),op;
  oppv = array(char,sizeof(opp));
  mem_copy,mem_base(oppv),opp;
  actionv   = strchar(opp.action);
  av = a_struct2vec(a);
  shm_free,shmkey,"opv";
  shm_free,shmkey,"oppv";
  shm_free,shmkey,"actionv";
  shm_free,shmkey,"av";
  shm_free,shmkey,"lmfitit";
  shm_free,shmkey,"modes";
  shm_write,shmkey,"opv",&opv;
  shm_write,shmkey,"oppv",&oppv;
  shm_write,shmkey,"actionv",&actionv;
  shm_write,shmkey,"av",&av;
  shm_write,shmkey,"mircube",&mircube;
  if (*opp.modes!=[]) shm_write,shmkey,"modes",opp.modes;
  else shm_write,shmkey,"modes",&([0]);
  shm_write,shmkey,"lmfitit",&([lmfititer_pass,lmfit_itmax,lmfititer]);
 }

func read_data_from_master(&a,&op,&opp,&mircube)
{
  extern lmfititer_pass,lmfit_itmax,lmfititer,emodes;
  opv  = shm_read(shmkey,"opv");
  oppv = shm_read(shmkey,"oppv");
  mem_copy,mem_base(op),opv;
  mem_copy,mem_base(opp),oppv;
  opp.action = strchar(shm_read(shmkey,"actionv"));
  av = shm_read(shmkey,"av");
  a  = a_vec2struct(av,opp);
  mircube = shm_read(shmkey,"mircube");
  tmp = shm_read(shmkey,"lmfitit");
  lmfititer_pass = tmp(1);
  lmfit_itmax = tmp(2);
  lmfititer   = tmp(3);
  // I don't know what's going on with opp.modes but I get seg violation
  // if I try to opp.modes = &emodes. here or up in calling routine.
  // hence I have to pass it in extern.
  emodes = shm_read(shmkey,"modes");
}
