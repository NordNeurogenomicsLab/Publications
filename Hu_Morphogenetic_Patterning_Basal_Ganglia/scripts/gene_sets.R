# Gene sets transcribed from the August 2025 analysis section of the methods notebook.
# The neuronal panel emphasizes specification, subtype identity, and maturation rather
# than the most variable genes dominated by unrelated programs. It is an embedding
# feature set, not a declaration that every selected cell is a projection neuron.
# The progenitor panel includes neuroepithelial/radial-glial, neurogenic, signaling,
# metabolic, and cell-cycle programs to resolve developmental progenitor states.
# Preserve the source list (including duplicate symbols) in config/ for inspection.
# unique() prevents duplicate features from receiving artificial extra weight.
progenitor_genes <- c(
  "Adamts1", "Adamts5", "Adcy2", "Aff2", "Akt3", "Aldoa", "Arl5a", "Ascl1",
  "Auts2", "Bach2", "Bcan", "Bcl11a", "Btf3", "Btf3", "Ccn1", "Ccnd1",
  "Cd24a", "Cdca2", "Cdk14", "Cdk6", "Cited2", "Cntnap2", "Crabp2", "Csnk1g2",
  "Csnk2b", "Csnk2b", "Cst3", "Ctnna2", "Cyc1", "Cyc1", "Dab1", "Dach1",
  "Dach2", "Dbi", "Dcc", "Dcun1d5", "Dlk1", "Dll1", "Dlx1", "Dlx2",
  "Dlx5", "Dlx6", "Dner", "Dnmt1", "Dpf3", "E2f7", "Ednrb", "Efna2",
  "Efna2", "Efna5", "Enc1", "Epha5", "Ephb1", "Etv1", "Ezr", "Fabp7",
  "Fat3", "Fgf15", "Fgfr2", "Fgfr3", "Fhl1", "Foxn3", "Foxn3", "Fst",
  "Fstl1", "Fzd5", "Fzd8", "Glcci1", "Gli1", "Gli2", "Gli3", "Gng5",
  "Gria4", "Grik2", "Grik3", "Gsk3b", "Hes1", "Hes5", "Hes6", "Hmgn1",
  "Hspd1", "Id4", "Igf1r", "Igfbp5", "Insm1", "Islr2", "Jmjd1c", "Jmjd1c",
  "Jund", "Kif13a", "Kif15", "Ldha", "Lgals1", "Lhx2", "Lhx6", "Lhx8",
  "Lingo2", "Lmo1", "Lsamp", "Lxn", "Maml2", "Mapk10", "Mdh2", "Mdh2",
  "Meis2", "Mfge8", "Mif", "Mta1", "Mybbp1a", "Nav2", "Nav2", "Ncam1",
  "Ncor2", "Ndn", "Ndrg2", "Neil3", "Nes", "Nfia", "Nfib", "Nfix",
  "Nkain2", "Nkx2-1", "Nkx2-2", "Nkx6-2", "Nmral1", "Nmral1", "Nnat", "Notch1",
  "Npas1", "Npas3", "Npm1", "Nr2f1", "Nr2f2", "Nrxn3", "Ntn1", "Olig1",
  "Olig2", "Otx2", "Pbrm1", "Pbx1", "Pcdh10", "Pcdh7", "Pcgf6", "Pdgfra",
  "Pdpn", "Peg10", "Peg3", "Pgam1", "Pgk1", "Pkm", "Plce1", "Plk2",
  "Plxna2", "Pnn", "Polr1d", "Polr2h", "Polr2m", "Polr3k", "Pou3f1", "Pou3f3",
  "Pou3f4", "Prc1", "Prox1", "Prox1", "Ptch1", "Ptn", "Ptprm", "Ptprz1",
  "Rbp1", "Rere", "Robo1", "Rpl26", "Rps4x", "Rps8", "Runx1t1", "Sall3",
  "Satb1", "Sema5b", "Sema6a", "Sfrp1", "Shh", "Six3", "Slc1a3", "Slit1",
  "Slit2", "Snf8", "Sox10", "Sox21", "Sox4", "Sox6", "Spag5", "Sparc",
  "St18", "Sulf1", "Syt1", "Taf7", "Tal2", "Tcf12", "Tcf4", "Tcf7l2",
  "Tead1", "Tfam", "Tgfb3", "Tle1", "Top1", "Tox3", "Tpi1", "Trps1",
  "Ttyh1", "Uchl1", "Vegfb", "Wnt7b", "Xbp1", "Zac1", "Plagl1", "Zbtb20",
  "Zeb1", "Zeb2", "Zfp462", "Zic1", "Zswim5"
)

proj_neu_genes <- c(
  "Ascl1", "Asic2", "Calb1", "Chl1", "Cntnap3", "Cxcl14", "Dlx1", "Dmrtb1",
  "Erbb4", "Etv1", "Foxp1", "Foxp2", "Gbx1", "Gbx2", "Gria1", "Grik1",
  "Id4", "Isl1", "Kctd12", "Lhx8", "Maf", "Mafb", "Mef2c", "Nfia",
  "Nfib", "Nkx2-1", "Npas1", "Npy", "Nr2f1", "Nr2f2", "Ntrk1", "Nxph1",
  "Pam", "Pbx1", "Peg10", "Pnoc", "Pou3f1", "Ptprm", "Reln", "Shh",
  "Six3", "Slit2", "Sox6", "Sphkap", "Sst", "Tac1", "Tcf4", "Tgfb3",
  "Th", "Unc5c", "Zeb2", "Zic1"
)
progenitor_genes <- unique(progenitor_genes)
proj_neu_genes <- unique(proj_neu_genes)
