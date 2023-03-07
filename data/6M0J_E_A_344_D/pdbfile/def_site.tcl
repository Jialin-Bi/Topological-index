mol new {pdbfile/6M0J_mut.pdb} type {pdb} first 0 last 0 step 1 waitfor 1
set prot [atomselect top "within 10.0 of chain E"]
$prot writepdb 6M0J_mut_b.pdb
set prot2 [atomselect top "within 10.0 of resid 344 and chain E"]
$prot2 writepdb 6M0J_mut_m.pdb
exit