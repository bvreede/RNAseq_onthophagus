library(DESeq2)

load("/Volumes/HD1v2/dsxRNAseq/barbara/onthophagus/dds_sv_20170515.rda")
load("/Volumes/HD1v2/dsxRNAseq/barbara/onthophagus/dds_svps_20170515.rda")

# dsx mediated - whole
# headhorns
lmh.sv <- results(dds.sv, contrast=c("group","largemalectrlHH","largemaledsxHH"))
smh.sv <- results(dds.sv, contrast=c("group","smallmalectrlHH","smallmaledsxHH"))
lfh.sv <- results(dds.sv, contrast=c("group","largefemalectrlHH","largefemaledsxHH"))
sfh.sv <- results(dds.sv, contrast=c("group","smallfemalectrlHH","smallfemaledsxHH"))

# thoracic horns
lmt.sv <- results(dds.sv, contrast=c("group","largemalectrlTH","largemaledsxTH"))
smt.sv <- results(dds.sv, contrast=c("group","smallmalectrlTH","smallmaledsxTH"))
lft.sv <- results(dds.sv, contrast=c("group","largefemalectrlTH","largefemaledsxTH"))
sft.sv <- results(dds.sv, contrast=c("group","smallfemalectrlTH","smallfemaledsxTH"))

# genitalia
lmg.sv <- results(dds.sv, contrast=c("group","largemalectrlGEN","largemaledsxGEN"))
smg.sv <- results(dds.sv, contrast=c("group","smallmalectrlGEN","smallmaledsxGEN"))
lfg.sv <- results(dds.sv, contrast=c("group","largefemalectrlGEN","largefemaledsxGEN"))
sfg.sv <- results(dds.sv, contrast=c("group","smallfemalectrlGEN","smallfemaledsxGEN"))

# brains
lmb.sv <- results(dds.sv, contrast=c("group","largemalectrlBR","largemaledsxBR"))
smb.sv <- results(dds.sv, contrast=c("group","smallmalectrlBR","smallmaledsxBR"))
lfb.sv <- results(dds.sv, contrast=c("group","largefemalectrlBR","largefemaledsxBR"))
sfb.sv <- results(dds.sv, contrast=c("group","smallfemalectrlBR","smallfemaledsxBR"))
save.image(file="/Volumes/HD1v2/dsxRNAseq/barbara/onthophagus/sv_v_pseudo.RData")



# sex bias
# headhorns
lch.sv <- results(dds.sv, contrast=c("group","largefemalectrlHH","largemalectrlHH"))
sch.sv <- results(dds.sv, contrast=c("group","smallfemalectrlHH","smallmalectrlHH"))
ldh.sv <- results(dds.sv, contrast=c("group","largefemaledsxHH","largemaledsxHH"))
sdh.sv <- results(dds.sv, contrast=c("group","smallfemaledsxHH","smallmaledsxHH"))

# thoracic horns
lct.sv <- results(dds.sv, contrast=c("group","largefemalectrlTH","largemalectrlTH"))
sct.sv <- results(dds.sv, contrast=c("group","smallfemalectrlTH","smallmalectrlTH"))
ldt.sv <- results(dds.sv, contrast=c("group","largefemaledsxTH","largemaledsxTH"))
sdt.sv <- results(dds.sv, contrast=c("group","smallfemaledsxTH","smallmaledsxTH"))

# genitalia
lcg.sv <- results(dds.sv, contrast=c("group","largefemalectrlGEN","largemalectrlGEN"))
scg.sv <- results(dds.sv, contrast=c("group","smallfemalectrlGEN","smallmalectrlGEN"))
ldg.sv <- results(dds.sv, contrast=c("group","largefemaledsxGEN","largemaledsxGEN"))
sdg.sv <- results(dds.sv, contrast=c("group","smallfemaledsxGEN","smallmaledsxGEN"))

# brains
lcb.sv <- results(dds.sv, contrast=c("group","largefemalectrlBR","largemalectrlBR"))
scb.sv <- results(dds.sv, contrast=c("group","smallfemalectrlBR","smallmalectrlBR"))
ldb.sv <- results(dds.sv, contrast=c("group","largefemaledsxBR","largemaledsxBR"))
sdb.sv <- results(dds.sv, contrast=c("group","smallfemaledsxBR","smallmaledsxBR"))
save.image(file="/Volumes/HD1v2/dsxRNAseq/barbara/onthophagus/sv_v_pseudo.RData")


# dsx mediated - pseudo
#headhorns
lmh.ps <- results(dds.sv.ps, contrast=c("group","largemalectrlHH","largemaledsxHH"))
smh.ps <- results(dds.sv.ps, contrast=c("group","smallmalectrlHH","smallmaledsxHH"))
lfh.ps <- results(dds.sv.ps, contrast=c("group","largefemalectrlHH","largefemaledsxHH"))
sfh.ps <- results(dds.sv.ps, contrast=c("group","smallfemalectrlHH","smallfemaledsxHH"))

#thoracic horns
lmt.ps <- results(dds.sv.ps, contrast=c("group","largemalectrlTH","largemaledsxTH"))
smt.ps <- results(dds.sv.ps, contrast=c("group","smallmalectrlTH","smallmaledsxTH"))
lft.ps <- results(dds.sv.ps, contrast=c("group","largefemalectrlTH","largefemaledsxTH"))
sft.ps <- results(dds.sv.ps, contrast=c("group","smallfemalectrlTH","smallfemaledsxTH"))

#genitalia
lmg.ps <- results(dds.sv.ps, contrast=c("group","largemalectrlGEN","largemaledsxGEN"))
smg.ps <- results(dds.sv.ps, contrast=c("group","smallmalectrlGEN","smallmaledsxGEN"))
lfg.ps <- results(dds.sv.ps, contrast=c("group","largefemalectrlGEN","largefemaledsxGEN"))
sfg.ps <- results(dds.sv.ps, contrast=c("group","smallfemalectrlGEN","smallfemaledsxGEN"))

#brains
lmb.ps <- results(dds.sv.ps, contrast=c("group","largemalectrlBR","largemaledsxBR"))
smb.ps <- results(dds.sv.ps, contrast=c("group","smallmalectrlBR","smallmaledsxBR"))
lfb.ps <- results(dds.sv.ps, contrast=c("group","largefemalectrlBR","largefemaledsxBR"))
sfb.ps <- results(dds.sv.ps, contrast=c("group","smallfemalectrlBR","smallfemaledsxBR"))
save.image(file="/Volumes/HD1v2/dsxRNAseq/barbara/onthophagus/sv_v_pseudo.RData")

# sex bias - pseudo
#headhorns
lch.ps <- results(dds.sv.ps, contrast=c("group","largefemalectrlHH","largemalectrlHH"))
sch.ps <- results(dds.sv.ps, contrast=c("group","smallfemalectrlHH","smallmalectrlHH"))
ldh.ps <- results(dds.sv.ps, contrast=c("group","largefemaledsxHH","largemaledsxHH"))
sdh.ps <- results(dds.sv.ps, contrast=c("group","smallfemaledsxHH","smallmaledsxHH"))

#thoracic horns
lct.ps <- results(dds.sv.ps, contrast=c("group","largefemalectrlTH","largemalectrlTH"))
sct.ps <- results(dds.sv.ps, contrast=c("group","smallfemalectrlTH","smallmalectrlTH"))
ldt.ps <- results(dds.sv.ps, contrast=c("group","largefemaledsxTH","largemaledsxTH"))
sdt.ps <- results(dds.sv.ps, contrast=c("group","smallfemaledsxTH","smallmaledsxTH"))

#genitalia
lcg.ps <- results(dds.sv.ps, contrast=c("group","largefemalectrlGEN","largemalectrlGEN"))
scg.ps <- results(dds.sv.ps, contrast=c("group","smallfemalectrlGEN","smallmalectrlGEN"))
ldg.ps <- results(dds.sv.ps, contrast=c("group","largefemaledsxGEN","largemaledsxGEN"))
sdg.ps <- results(dds.sv.ps, contrast=c("group","smallfemaledsxGEN","smallmaledsxGEN"))

#brains
lcb.ps <- results(dds.sv.ps, contrast=c("group","largefemalectrlBR","largemalectrlBR"))
scb.ps <- results(dds.sv.ps, contrast=c("group","smallfemalectrlBR","smallmalectrlBR"))
ldb.ps <- results(dds.sv.ps, contrast=c("group","largefemaledsxBR","largemaledsxBR"))
sdb.ps <- results(dds.sv.ps, contrast=c("group","smallfemaledsxBR","smallmaledsxBR"))
save.image(file="/Volumes/HD1v2/dsxRNAseq/barbara/onthophagus/sv_v_pseudo.RData")

# nutrition
# headhorns
cmh.sv <- results(dds.sv, contrast=c("group","largemalectrlHH","smallmalectrlHH"))
cfh.sv <- results(dds.sv, contrast=c("group","largefemalectrlHH","smallfemalectrlHH"))
dmh.sv <- results(dds.sv, contrast=c("group","largemaledsxHH","smallmaledsxHH"))
dfh.sv <- results(dds.sv, contrast=c("group","largefemaledsxHH","smallfemaledsxHH"))

# thoracic horns
cmt.sv <- results(dds.sv, contrast=c("group","largemalectrlTH","smallmalectrlTH"))
cft.sv <- results(dds.sv, contrast=c("group","largefemalectrlTH","smallfemalectrlTH"))
dmt.sv <- results(dds.sv, contrast=c("group","largemaledsxTH","smallmaledsxTH"))
dft.sv <- results(dds.sv, contrast=c("group","largefemaledsxTH","smallfemaledsxTH"))

# genitalia
cmg.sv <- results(dds.sv, contrast=c("group","largemalectrlGEN","smallmalectrlGEN"))
cfg.sv <- results(dds.sv, contrast=c("group","largefemalectrlGEN","smallfemalectrlGEN"))
dmg.sv <- results(dds.sv, contrast=c("group","largemaledsxGEN","smallmaledsxGEN"))
dfg.sv <- results(dds.sv, contrast=c("group","largefemaledsxGEN","smallfemaledsxGEN"))

# brains
cmb.sv <- results(dds.sv, contrast=c("group","largemalectrlBR","smallmalectrlBR"))
cfb.sv <- results(dds.sv, contrast=c("group","largefemalectrlBR","smallfemalectrlBR"))
dmb.sv <- results(dds.sv, contrast=c("group","largemaledsxBR","smallmaledsxBR"))
dfb.sv <- results(dds.sv, contrast=c("group","largefemaledsxBR","smallfemaledsxBR"))
save.image(file="/Volumes/HD1v2/dsxRNAseq/barbara/onthophagus/sv_v_pseudo.RData")

#headhorns
cmh.ps <- results(dds.sv.ps, contrast=c("group","largemalectrlHH","smallmalectrlHH"))
cfh.ps <- results(dds.sv.ps, contrast=c("group","largefemalectrlHH","smallfemalectrlHH"))
dmh.ps <- results(dds.sv.ps, contrast=c("group","largemaledsxHH","smallmaledsxHH"))
dfh.ps <- results(dds.sv.ps, contrast=c("group","largefemaledsxHH","smallfemaledsxHH"))

#thoracic horns
cmt.ps <- results(dds.sv.ps, contrast=c("group","largemalectrlTH","smallmalectrlTH"))
cft.ps <- results(dds.sv.ps, contrast=c("group","largefemalectrlTH","smallfemalectrlTH"))
dmt.ps <- results(dds.sv.ps, contrast=c("group","largemaledsxTH","smallmaledsxTH"))
dft.ps <- results(dds.sv.ps, contrast=c("group","largefemaledsxTH","smallfemaledsxTH"))

#genitalia
cmg.ps <- results(dds.sv.ps, contrast=c("group","largemalectrlGEN","smallmalectrlGEN"))
cfg.ps <- results(dds.sv.ps, contrast=c("group","largefemalectrlGEN","smallfemalectrlGEN"))
dmg.ps <- results(dds.sv.ps, contrast=c("group","largemaledsxGEN","smallmaledsxGEN"))
dfg.ps <- results(dds.sv.ps, contrast=c("group","largefemaledsxGEN","smallfemaledsxGEN"))

#brains
cmb.ps <- results(dds.sv.ps, contrast=c("group","largemalectrlBR","smallmalectrlBR"))
cfb.ps <- results(dds.sv.ps, contrast=c("group","largefemalectrlBR","smallfemalectrlBR"))
dmb.ps <- results(dds.sv.ps, contrast=c("group","largemaledsxBR","smallmaledsxBR"))
dfb.ps <- results(dds.sv.ps, contrast=c("group","largefemaledsxBR","smallfemaledsxBR"))

save.image(file="/Volumes/HD1v2/dsxRNAseq/barbara/onthophagus/sv_v_pseudo.RData")