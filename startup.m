%%***********************************************************************
%%
%%
%%***********************************************************************

  restoredefaultpath;
  addpath(genpath(pwd),path);  
%%
  set(0,'defaultaxesfontsize',14);
  set(0,'defaultlinemarkersize',6);
  set(0,'defaulttextfontsize',14);
%%***********************************************************************

run(sprintf('%s/startup.m', userpath))

A = shrikhande(4);
B = lineGraph(completeBipartite(4, 4));
