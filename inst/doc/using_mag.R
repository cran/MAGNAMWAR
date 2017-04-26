## ----echo = F------------------------------------------------------------
library(MAGNAMWAR)

## ------------------------------------------------------------------------
lfrc_fasta <- system.file('extdata', 'RASTtogbk//lfrc.fasta', package='MAGNAMWAR')
lfrc_reference <- system.file('extdata', 'RASTtogbk//lfrc_lookup.csv', package='MAGNAMWAR')
lfrc_path <- system.file('extdata', 'RASTtogbk//lfrc_out.fasta', package='MAGNAMWAR')

RASTtogbk(lfrc_fasta,lfrc_reference,lfrc_path)

## ---- echo=FALSE---------------------------------------------------------
file_groups <- system.file('extdata', 'groups_example_r.txt', package='MAGNAMWAR')

## ---- results="hide"-----------------------------------------------------
# file_groups is the file path to your output file from OrthoMCL

parsed_data <- format_afterOrtho(file_groups, "groups")

## ---- echo = F-----------------------------------------------------------
parsed_data[[1]][,1:13]

## ---- echo = F-----------------------------------------------------------
parsed_data[[2]][,1:2]

## ---- echo = F-----------------------------------------------------------
head(pheno_data,n = 5)

## ---- results="hide"-----------------------------------------------------
mcl_matrix <- analyze_OrthoMCL(parsed_data, 
                               pheno_data, 
                               model = 'lmeR2nest', 
                               species_name = 'Treatment',  
                               resp = 'RespVar', 
                               rndm1 = 'Experiment', 
                               rndm2 = 'Vial')

## ---- echo=F-------------------------------------------------------------
mcl_matrix[,1:3]

## ------------------------------------------------------------------------
dir <- system.file('extdata', 'fasta_dir', package='MAGNAMWAR')
dir <- paste(dir,'/',sep='')
joined_matrix <- join_repseq(after_ortho_format_grps, dir, mcl_mtrx_grps, fastaformat = 'old')

## ---- echo=F-------------------------------------------------------------
joined_matrix[1,8:10]

## ---- fig.width=8, fig.height=4------------------------------------------
pdgplot(pheno_data, mcl_matrix, OG = 'ex_r00002', species_colname = 'Treatment', data_colname = 'RespVar', ylab = "TAG Content")

## ---- fig.width=8, fig.height=4------------------------------------------
tree <- system.file('extdata', 'muscle_tree2.dnd', package='MAGNAMWAR')

pdgplot(pheno_data, mcl_matrix, 'ex_r00002', 'Treatment', 'RespVar', ylab = "TAG Content", tree = tree)

## ---- fig.width=6, fig.height=6.5, fig.align='center'--------------------
phydataerror(tree, pheno_data, mcl_matrix, species_colname = 'Treatment', data_colname = 'RespVar', OG='ex_r00002', xlabel='TAG Content')

## ------------------------------------------------------------------------
pdg_v_OG(parsed_data,0)

## ---- fig.width=4, fig.height=4, fig.align='center'----------------------
qqplotter(mcl_matrix)

## ---- fig.width=6, fig.height=3, fig.align='center'----------------------
manhat_grp(parsed_data, mcl_matrix)

